rule minimap2:
    resources:
        threads = 8,
        nodes = 16
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

rule minimap2_with_index:
    resources:
        threads = 8,
        nodes = 32
    shell:
        """
        module load minimap2
        minimap2 \
            -ax splice \
            -t {resources.threads} \
            -a {input.ind} \
            --MD \
            --secondary=no \
            -L \
            -o {output.sam} \
            {input.fq} \

        """

rule minimap2_index:
    resources:
        threads = 8,
        nodes = 16
    shell:
        """
        module load minimap2
        minimap2 \
            -d {output.ind} \
            {input.fa}
        """
