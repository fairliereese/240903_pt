rule gzip:
    resources:
        threads = 1,
        nodes = 1
    shell:
        """
        gzip -c {input.ifile} > {output.ofile}
        """

rule gunzip:
    resources:
        threads = 1,
        nodes = 1
    shell:
        """
        gunzip -c {input.gz} > {output.ofile}
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

rule bam_to_fastq:
    resources:
        threads = 1,
        nodes = 1,
        time = "1:00:00"
    shell:
        """
        module load bedtools
        bedtools bamtofastq \
            -i {input.bam} \
            -fq {output.fq}
        """

rule spliced_bam2gff:
    resources:
        threads = 8,
        nodes = 1
    conda:
        'base'
    shell:
        """
        conda activate /gpfs/projects/bsc83/utils/conda_envs/spliced_bam2gff
        spliced_bam2gff \
            -M {input.align} \
            -g \
            -d 20 \
            -t {resources.threads} > {output.gff}
        """
