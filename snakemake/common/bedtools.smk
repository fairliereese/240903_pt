
rule fa_to_bed_windows:
    resources:
        threads = 1,
        nodes = 1
    shell:
        """
        module load bedtools
        bedtools makewindows -g {input.chrom_sizes} \
            -w {params.window_size} \
            > {output.bed}
        """

rule bed_to_fasta:
    resources:
        threads = 1,
        nodes = 1
    shell:
        """
        module load bedtools
        bedtools getfasta -fi {input.fa} -bed {input.bed} > {output.fa}
        """

# rule get_gc:
#     conda:
#         'ucsctools'
#     resources:
#         threads = 1,
#         nodes = 1
#     shell:
#         """
#         bedtools makewindows \
#             -g {input.chrom_sizes} \
#             -w {params.window_size} | \
#         bedtools nuc \
#             -fi {input.fa} \
#             -bed - | \
#         awk 'NR>1 {print $1"\t"$2"\t"$3"\t"$5}' > {output.bedgraph}
#         bedGraphToBigWig \
#             {output.bedgraph} \
#             {input.chrom_sizes} \
#             {output.bw}
#         """
rule samtobed:
    resources:
        threads = 1,
        nodes = 1
    shell:
        """
        module load samtools
        module load bedtools
        samtools view -b {input.sam} | bedtools bamtobed -i - > {output.bed}
        """
