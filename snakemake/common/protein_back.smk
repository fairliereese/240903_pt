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

rule find_complete_orfs:
    resources:
        threads = 1,
        nodes = 2
    shell:
        """
        python {params.scripts_path}/check_orf_completeness.py \
            --cpat_seqs {input.cpat_seqs} \
            --orfanage_seqs {input.orfanage_seqs} \
            --cpat_info {input.cpat_info} \
            --orfanage_info {input.orfanage_info} \
            --output_path {output.tsv}
        """

rule get_cpat_cds_coords:
    resources:
        threads = 1,
        nodes = 2
    shell:
        """
        python {params.scripts_path}/create_cpat_CDS_coordinates.py \
                                --name {params.opref} \
                                --sample_gtf {input.sample_gtf} \
                                --called_orfs {input.called_orfs}
        """


rule postprocess_combine_cds_gtf:,
    params:
        threads = 1,
        nodes = 1
    conda:
        'base'
    shell:
        """
        conda activate /gpfs/projects/bsc83/utils/conda_envs/gffread
        cat {input.cpat_cds} {input.orfanage_cds} | gffread -T - | sort -k1,1V -k4,4n -k5,5rn -k3,3r | gffread -T - > {output.gtf}
        """

rule prep_make_query:
    input:
        query="data-raw/{sample}.gtf",
        genome="data-raw/{sample}_genome.fa",
    output:
        "results/{sample}_orfanage_ready_query.gtf",
    container:
        f"docker://condaforge/mambaforge:{config['mambaforge_version']}"
    conda:
        "envs/gffread.yml"
    shell:
        "cat {input.query} | gffread -g {input.genome} -T -o {output} -"

rule postprocess_amend_cds_source:
    input:
        cds=config['lr']['orfanage']['protein_unsourced'],
        source_gtf="results/{sample}_orfanage_ready_query.gtf",
        cpat_cds="results/{sample}_cpat_with_cds.gtf",
        orfanage_cds="results/{sample}_orfanage_cds_filtered_stop_codon_corrected.gtf",
    output:
        "results/protein_prediction/{sample}_protein.gtf",
    shell:
        """python {params.scripts_path}/recover_source.py \
                                --combined_cds_path {input.cds} \
                                --source_gtf_path {input.source_gtf} \
                                --cpat_cds_path {input.cpat_cds} \
                                --orfanage_cds_path {input.orfanage_cds} \
                                --output_path {output}"""


# rule postprocess_extract_orf_fasta:
#     input:
#         protein_gtf="results/protein_prediction/{sample}_protein.gtf",
#         genome="data-raw/{sample}_genome.fa",
#     output:
#         "results/protein_prediction/{sample}_ORF.fa",
#     container:
#         f"docker://condaforge/mambaforge:{config['mambaforge_version']}"
#     conda:
#         "envs/gffread.yml"
#     shell:
#         "gffread -x {output} -g {input.genome} {input.protein_gtf}"
#
#
# rule postprocess_extract_protein_fasta:
#     input:
#         protein_gtf="results/protein_prediction/{sample}_protein.gtf",
#         genome="data-raw/{sample}_genome.fa",
#     output:
#         "results/protein_prediction/{sample}_protein.fa",
#     container:
#         f"docker://condaforge/mambaforge:{config['mambaforge_version']}"
#     conda:
#         "envs/gffread.yml"
#     shell:
#         "gffread -y {output} -g {input.genome} {input.protein_gtf}"
#
#
# rule postprocess_prepare_sqanti_protein_gtf:
#     input:
#         protein_gtf="results/protein_prediction/{sample}_protein.gtf",
#         annotation="data-raw/{sample}_annotation.gtf",
#     output:
#         "results/{sample}_gencode.transcript_exons_only.gtf",
#         "results/{sample}_gencode.cds_renamed_exon.gtf",
#         "results/{sample}.transcript_exons_only.gtf",
#         "results/{sample}.cds_renamed_exon.gtf",
#     container:
#         f"docker://condaforge/mambaforge:{config['mambaforge_version']}"
#     conda:
#         "envs/proteogenomics.yml"
#     threads: config["per_rule_threads_multi"]
#     shell:
#         """python workflow/scripts/rename_cds_to_exon.py \
#             --sample_gtf {input.protein_gtf} \
#             --sample_name results/{wildcards.sample} \
#             --reference_gtf {input.annotation} \
#             --num_cores {threads}"""
#
#
# rule postprocess_prepare_sqanti_protein_tsv:
#     input:
#         transcript_only_exons="{sample}.transcript_exons_only.gtf",
#         cds_renamed="{sample}.cds_renamed_exon.gtf",
#     output:
#         "{sample}_best_orf.tsv",
#     container:
#         f"docker://condaforge/mambaforge:{config['mambaforge_version']}"
#     conda:
#         "envs/default.yml"
#     shell:
#         """python workflow/scripts/create_orf_table_for_sqanti_protein.py \
#                     --transcript_exons_path {input.transcript_only_exons} \
#                     --cds_only_path {input.cds_renamed} \
#                     --output_prefix {wildcards.sample}"""
#
#
# rule postprocess_run_sqanti_protein:
#     input:
#         best_orfs="results/{sample}_best_orf.tsv",
#         renamed_exons="results/{sample}.transcript_exons_only.gtf",
#         cds_only="results/{sample}.cds_renamed_exon.gtf",
#         gencode_renamed_exons="results/{sample}_gencode.transcript_exons_only.gtf",
#         gencode_cds_only="results/{sample}_gencode.cds_renamed_exon.gtf",
#     output:
#         "results/{sample}.sqanti_protein_classification.tsv",
#     container:
#         f"docker://condaforge/mambaforge:{config['mambaforge_version']}"
#     conda:
#         "envs/proteogenomics.yml"
#     shell:
#         """python workflow/scripts/sqanti3_protein.py \
#                     {input.renamed_exons} \
#                     {input.cds_only} \
#                     {input.best_orfs} \
#                     {input.gencode_renamed_exons} \
#                     {input.gencode_cds_only} \
#                     -d ./ \
#                     -p results/{wildcards.sample}"""
#
#
# rule postprocess_summarize_all:
#     input:
#         best_orf="results/{sample}_best_orf.tsv",
#         protein_classification="results/{sample}.sqanti_protein_classification.tsv",
#         orf_completeness="results/{sample}_ORF_completeness.tsv",
#         original_gtf="data-raw/{sample}.gtf",
#         gtf_predicted="results/protein_prediction/{sample}_protein.gtf",
#         protein_fasta="results/protein_prediction/{sample}_protein.fa",
#         blastp="results/protein_prediction/{sample}_{gencode_vers}_blastp.out",
#     output:
#         "results/protein_prediction/{sample}_protein_annotation_{gencode_vers}.tsv",
#     container:
#         f"docker://condaforge/mambaforge:{config['mambaforge_version']}"
#     conda:
#         "envs/default.yml"
#     shell:
#         """python workflow/scripts/create_protein_overview_table.py \
#             --best_orf_path {input.best_orf} \
#             --sqanti_protein_path {input.protein_classification} \
#             --orf_completeness_path {input.orf_completeness} \
#             --output_name {output} \
#             --gtf_original_path {input.original_gtf} \
#             --gtf_predicted_path {input.gtf_predicted} \
#             --protein_fasta_path {input.protein_fasta} \
#             --blastp_path {input.blastp}
#             """
#
#
# rule postprocess_prepare_protein_fasta_for_blast:
#     input:
#         protein_fasta="data-raw/gencode.{gencode_vers}.pc_translations.fa",
#     output:
#         "results/gencode.{gencode_vers}.pc_translations_renamed.fa",
#     container:
#         f"docker://condaforge/mambaforge:{config['mambaforge_version']}"
#     conda:
#         "envs/blast.yml"
#     shell:
#         "sed -r 's/\|[^\|]*//2g' {input.protein_fasta} > {output}"
#
#
# rule postprocess_create_blast_db:
#     input:
#         protein_fasta="results/gencode.{gencode_vers}.pc_translations_renamed.fa",
#     output:
#         "results/gencode.{gencode_vers}.pc_translations_renamed.pog",
#     container:
#         f"docker://condaforge/mambaforge:{config['mambaforge_version']}"
#     conda:
#         "envs/blast.yml"
#     shell:
#         """makeblastdb \
#                 -in {input.protein_fasta} \
#                 -dbtype prot \
#                 -parse_seqids \
#                 -out results/gencode.{wildcards.gencode_vers}.pc_translations_renamed"""
#
#
# rule postprocess_run_blast:
#     input:
#         protein_fasta="results/protein_prediction/{sample}_protein.fa",
#         protein_reference="results/gencode.{gencode_vers}.pc_translations_renamed.fa",
#         dbs="results/gencode.{gencode_vers}.pc_translations_renamed.pog",
#     output:
#         "results/protein_prediction/{sample}_{gencode_vers}_blastp.out",
#     container:
#         f"docker://condaforge/mambaforge:{config['mambaforge_version']}"
#     conda:
#         "envs/blast.yml"
#     params:
#         blast_evalue=config["blast_e_value"],
#     threads: config["per_rule_threads_multi"]
#     shell:
#         """blastp \
#             -evalue {params.blast_evalue} \
#             -num_threads {threads} \
#             -outfmt 6 \
#             -db results/gencode.{wildcards.gencode_vers}.pc_translations_renamed \
#             -query {input.protein_fasta} > \
#             {output}"""
