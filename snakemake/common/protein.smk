rule orfanage_find_orfs:
    resources:
        threads = 16,
        nodes = 3
    conda:
        'base'
    shell:
        """
        conda activate /gpfs/projects/bsc83/utils/conda_envs/orfanage
        orfanage \
            --cleanq \
            --mode LONGEST_MATCH \
            --reference {input.fa} \
            --query {input.gtf} \
            --output {output.gtf} \
            --threads {resources.threads} \
            {input.annot_gtf}
        """

rule orfanage_filter:
    resources:
        threads = 1,
        nodes = 2
    shell:
        """python {params.scripts_path}/filter_orfanage.py \
            --orfanage_gtf_file_path {input.cds} \
            --output_path_to_be_predicted {output.gtf_pred} \
            --output_path_filtered {output.gtf} \
            --minimum_orf_length {params.min_orf_len}"""


rule gtf_cds_gtf_stop_codon:
    resources:
        threads = 1,
        nodes = 2
    conda:
        'base'
    shell:
        """
        conda activate /gpfs/projects/bsc83/utils/conda_envs/gffread
        gffread -x {output.fa} -g {input.fa} {input.gtf}
        """

rule correct_stop_codon:
    resources:
        threads = 1,
        nodes = 2
    shell:
        """
        python {params.scripts_path}/correct_stop_codon_orfanage.py \
            --orfanage_gtf_file_path {input.gtf} \
            --output_path {output.gtf}
        """

rule cds_for_cpat:
    resources:
        threads = 1,
        nodes = 2
    conda:
        'base'
    shell:
        """
        conda activate /gpfs/projects/bsc83/utils/conda_envs/gffread
        gffread -w {output.fa} -g {input.fa} {input.cds}
        """

rule cpat:
    resources:
        threads = 1,
        nodes = 1
    conda:
        'base'
    shell:
        """
        module load gcc
        module load R/4.3.0
        cpat \
            -x {input.hexamer} \
            -d {input.logit_model} \
            -g {input.fa} \
            --min-orf={params.min_orf} \
            --top-orf={params.top_orf} \
            -o {params.opref} \
            1> {params.opref}_cpat.output \
            2> {params.opref}_cpat.error
        """

rule filt_cpat:
    resources:
        threads = 1,
        nodes = 1
    shell:
        """
        python {params.scripts_path}/filter_cpat.py \
            --input_file_path {input.prob} \
            --orf_input_seq_path {input.fa} \
            --output_path {output.tsv} \
            --first_cutoff {params.first_cutoff} \
            --second_cutoff {params.second_cutoff}
        """

# rule postprocess_check_orf_completeness:
#     input:
#         cpat_seqs="results/{sample}.ORF_seqs.fa",
#         orfanage_seqs="results/{sample}_orfanage_orfs.fa",
#         cpat_info="results/{sample}.ORF_remaining.tsv",
#         orfanage_info="results/{sample}_orfanage_cds_filtered_stop_codon_corrected.gtf",
#     output:
#         "results/{sample}_ORF_completeness.tsv",
#     container:
#         f"docker://condaforge/mambaforge:{config['mambaforge_version']}"
#     conda:
#         "envs/default.yml"
#     shell:
#         """python workflow/scripts/check_orf_completeness.py \
#             --cpat_seqs {input.cpat_seqs} \
#             --orfanage_seqs {input.orfanage_seqs} \
#             --cpat_info {input.cpat_info} \
#             --orfanage_info {input.orfanage_info} \
#             --output_path {output}"""
