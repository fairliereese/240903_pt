rule whatshap_polyphase:
    resources:
        threads = 8,
        nodes = 4
    conda:
        'whatshap'
    shell:
        """
        whatshap polyphase \
            {input.vcf} \
            {input.bam} \
            --ploidy {params.ploidy} \
            --reference {input.fa} \
            -o {output.vcf}
        """

rule whatshap_stats:
    resources:
        threads = 1,
        nodes = 1
    conda:
        'whatshap'
    shell:
        """
        whatshap stats \
            {input.vcf} \
            --chr-lengths {input.chr_lens} \
            --tsv {output.tsv} > {output.log}
        """
