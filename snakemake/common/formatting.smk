rule gzip:
    resources:
        threads = 1,
        nodes = 1
    shell:
        """
        gzip -c {input.gz} > {output.ofile}
        """

rule gunzip:
    resources:
        threads = 1,
        nodes = 1
    shell:
        """
        gunzip -c {input.ifile} > {output.ofile}
        """

rule concat_vcf:
    resources:
        threads = 1,
        nodes = 1
    conda:
        "base"
    shell:
        """
        bcftools concat {input.ifiles} > {output.ofile}
        """

rule bed_to_vcf:
    resources:
        threads = 1,
        nodes = 1
    run:
        bed_to_vcf(input.bed,
                   output.vcf,
                   input.template)

# get chromosome lengths using faidx
rule fa_get_chr_lens:
    resources:
        threads = 1,
        nodes = 1
    shell:
        """
        faidx {input.fa} \
            -i chromsizes \
            > {output.chr_lens}
        """
