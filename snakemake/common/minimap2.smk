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
            --MD \
            -o {output.sam} \
            -a {input.ind} \
            {input.fq}
        """

rule minimap2_index:
    resources:
        threads = 8,
        nodes = 16
    shell:
        """
        module load minimap2
        minimap2 -x map-ont -t 112 -d {output.ind} {input.fa}
        """
