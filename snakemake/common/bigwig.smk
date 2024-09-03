def get_strand_flag(wc):
    """
    Get strandedness argument for bigwig creation
    """
    if wc.strand == 'fwd':
        flag = '--filterRNAstrand reverse'
    elif wc.strand == 'rev':
        flag = '--filterRNAstrand forward'
    return flag

rule bam_to_bw_strand:
    resources:
        threads = 4,
        nodes = 1
    params:
        strand_flag = lambda wc:get_strand_flag(wc)
    shell:
        """
        bamCoverage -b {input.bam} -o {output.bw} {params.strand_flag}
        """

rule bam_to_bw:
    resources:
        threads = 4,
        nodes = 1
    conda:
        'deeptools'
    shell:
        """
        bamCoverage -b {input.bam} -o {output.bw}
        """
