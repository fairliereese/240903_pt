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
    conda:
        'base'
    shell:
        """
        conda activate /gpfs/projects/bsc83/utils/conda_envs/cpat
        cpat.py \
                -x {input.hexamer} \
                -d {input.logit_model} \
                -g {input.query} \
                --min-orf={params.min_orf} \
                --top-orf={params.top_orf} \
                -o {params.opref} \
                1> {params.opref}_cpat.output \
                2> {params.opref}_cpat.error
        """
