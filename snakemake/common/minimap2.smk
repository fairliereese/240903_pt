rule minimap2:
    resources:
        threads = 64,
        nodes = 4
    shell:
        """
        module load minimap2
        minimap2 \
            -ax splice \
            -t {resources.threads} \
            --MD \
            --secondary=no \
            -L \
            -o {output.sam} \
            {input.fq} \
            {input.fa}
        """
