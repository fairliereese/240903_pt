rule compute_gc:
    conda:
        'deeptools'
    resources:
        threads = 2,
        nodes = 1
    shell:
        """
        computeGCContent -r {input.fa} -w 5 -o {output.bw}
        """
