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
