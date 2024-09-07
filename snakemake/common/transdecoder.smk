rule find_orfs:
    resources:
        threads = 1,
        nodes = 2
    shell:
        """
        conda activate /gpfs/projects/bsc83/utils/conda_envs/transdecoder
        TransDecoder.LongOrfs \
            -t {input.fa} \
            -m 100 \
            -O {params.odir}
        """
