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
