rule fa_to_bam:
    resources:
        threads = 1,
        nodes = 1
    shell:
        """
        module load samtools
        samtools view -b -T {input.fa} {input.fa} > {input.bam}
        """

rule sam_to_bam:
    resources:
        threads = 16,
        nodes = 1,
        time = "1:00:00"
    shell:
        """
        module load samtools
        samtools view -hSb {input.sam} > {output.bam}
        """

rule add_rg:
    resources:
        threads = 4,
        nodes = 1
    shell:
        """
        module load samtools
        samtools addreplacerg \
            -r "SM:{wildcards.dataset}" \
            -r "ID:{wildcards.dataset}" \
            -m overwrite_all \
            -@ {resources.threads} \
            -o {output.align} \
            {input.align}
        """

rule sort_bam:
    resources:
        threads = 16,
        nodes = 1,
        time = "1:00:00"
    shell:
        """
        module load samtools
        samtools sort \
            --threads {resources.threads} \
            -O bam {input.bam} > {output.bam}
        """

rule index_bam:
    resources:
        threads = 16,
        nodes = 1,
        time = "1:00:00"
    shell:
        """
        module load samtools
        samtools index -@ {resources.threads} {input.bam}
        """

rule alignment_stats:
    resources:
        threads = 2,
        nodes = 1
    shell:
        """
        module load samtools
        samtools stats {input.alignment} | grep ^SN | cut -f 2- | grep -e 'reads map
ped' -e 'reads unmapped' -e 'average length' -e 'maximum length' | sed '/reads mapped and paired/d' > {output.stats}
        """

rule merge_alignment:
    resources:
        threads = 32,
        nodes = 1,
        time = "2:00:00"
    shell:
        """
        module load samtools
        samtools merge \
            -@ {resources.threads} \
            -o {output.bam} {input.files}
        """

rule count_primary_mappings:
    resources:
        threads = 8,
        nodes = 1
    shell:
        """
        module load samtools
        samtools view -F 256 -c {input.align} > {output.out}
        """

rule map_query_cov:
    resources:
        nodes = 1,
        threads = 16
    run:
        df = compute_query_coverage(input.align, resources.threads)
        df.to_csv(output.out, sep='\t', index=False)

rule map_query_cov_summary:
    resources:
        nodes = 2,
        threads = 1
    run:
        temp = pd.DataFrame()
        for d, f in zip(params.datasets, input.tsvs):
            temp2 = pd.read_csv(f, sep='\t')
            temp2['dataset'] = d
            temp = pd.concat([temp, temp2], axis=0)
        temp.to_csv(output.summ, sep='\t', index=False)

rule count_lines_summary:
    resources:
        threads = 1,
        nodes = 2
    run:
        temp = pd.DataFrame()
        for d, f in zip(params.datasets, input.tsvs):
            temp2 = pd.read_csv(f, sep='\t')
            temp3 = pd.DataFrame()
            temp3['dataset'] = [d]
            temp3['n_dupe_reads'] = [len(temp2.index)]
            temp = pd.concat([temp, temp3], axis=0)
        temp.to_csv(output.summ, sep='\t', index=False)



rule primary_mappings_filt:
    resources:
        threads = 8,
        nodes = 1,
        time = "1:00:00"
    shell:
        """
        module load samtools
        samtools view -h -F 256 {input.sam} > {output.out}
        """

# filter out non-primary, unmapped, and supp. alignments
rule filt_non_prim_unmap_supp:
    resources:
        nodes = 4,
        threads = 64
    shell:
        """
        module load samtools
        samtools view \
            -@ {resources.threads} \
            -h \
            -F 256 \
            -F 4 \
            -F 2048 \
            {input.align} > {output.align}
        """

# filter out non-primary, unmapped, and supp. alignments
rule filt_non_prim_unmap_supp_contigs:
    resources:
        nodes = 4,
        threads = 64
    shell:
        """
        module load samtools
        samtools view \
            -@ {resources.threads} \
            -h \
            -F 256 \
            -F 4 \
            -F 2048 \
            -L {input.bed} \
            {input.align} > {output.align}
        """

rule cov_filt_read_ids_min_max:
    resources:
        nodes = 2,
        threads = 1
    run:
        df = pd.read_csv(input.query_cov, sep='\t')
        df = df.loc[(df.query_cov>=params.min_cov)&(df.query_cov<params.max_cov)]
        df = df[['read_id']]
        df.to_csv(output.out, sep='\t', index=False, header=None)


rule cov_filt_read_ids:
    resources:
        nodes = 2,
        threads = 1
    run:
        df = pd.read_csv(input.query_cov, sep='\t')
        df = df.loc[df.query_cov>=params.min_cov]
        df = df[['read_id']]
        df.to_csv(output.out, sep='\t', index=False, header=None)

rule cov_filt:
    resources:
        nodes = 2,
        threads = 1
    shell:
        """
        module load samtools
        samtools view -h --qname-file {input.txt} {input.sam} > {output.out}
        """
#         """
#         module load samtools
#         temp_fname={wildcards.dataset}_query_cov_filt
#         awk -v dataset="{wildcards.dataset}" '$2 > 0.9 && $3 == dataset {{print $1}}' {input.query_cov} > ${{temp_fname}}
#         samtools view -h --qname-file ${{temp_fname}} {input.sam} > {output.out}
#         """

rule bam_reads_per_chr:
    resources:
        threads = 4,
        nodes = 1,
    shell:
        """
        module load samtools
        samtools idxstats {input.align} | cut -f 1,3 > {output.out}
        """

rule reads_per_chr_summary:
    resources:
        nodes = 1,
        threads = 1
    run:
        df = pd.DataFrame()
        for f, d in zip(input.tsvs, params.datasets):
            temp = pd.read_csv(f, sep='\t')
            temp.columns = ['chrom', 'n_reads']
            temp['dataset'] = d
            df = pd.concat([df, temp], axis=0)
        # df.columns = ['chr', 'n_reads', 'dataset']
        df.to_csv(output.summ, sep='\t', index=False)

# NM:i:count Number of differences (mismatches plus inserted and deleted bases) between the sequence and
# reference, counting only (case-insensitive) A, C, G and T bases in sequence and reference as potential
# matches, with everything else being a mismatch. Note this means that ambiguity codes in both
# sequence and reference that match each other, such as ‘N’ in both, or compatible codes such as ‘A’ and
# ‘R’, are still counted as mismatches. The special sequence base ‘=’ will always be considered to be a
# 3
# match, even if the reference is ambiguous at that point. Alignment reference skips, padding, soft and
# hard clipping (‘N’, ‘P’, ‘S’ and ‘H’ CIGAR operations) do not count as mismatches, but insertions and
# deletions count as one mismatch per base.
# Note that historically this has been ill-defined and both data and tools exist that disagree with this
# definition
# https://samtools.github.io/hts-specs/SAMtags.pdf
rule mm_per_read:
    resources:
        threads = 4,
        nodes = 2,
    run:
        if input.align.endswith('.bam'):
            in_mode = 'rb'
        else:
            in_mode = 'r'
        input =  pysam.AlignmentFile(input.align, in_mode, threads=resources.threads)
        read_mm = []
        read_ids = []
        for read in input:
            read_mm.append(read.get_tag('NM'))
            read_ids.append(read.query_name)
        df = pd.DataFrame()
        df['read_id'] = read_ids
        df['n_mm'] = read_mm
        df.to_csv(output.tsv, sep='\t', index=False)

rule mm_per_read_summary:
    resources:
        nodes = 1,
        threads = 1
    run:
        df = pd.DataFrame()
        for f, d in zip(input.tsvs, params.datasets):
            temp = pd.read_csv(f, sep='\t')
            temp.columns = ['read_id', 'n_mm']
            temp['dataset'] = d
            df = pd.concat([df, temp], axis=0)
        # df.columns = ['chr', 'n_reads', 'dataset']
        df.to_csv(output.summ, sep='\t', index=False)

rule supp_filt:
    resources:
        nodes = 1,
        threads = 8
    shell:
        """
        module load samtools
        samtools view -hF 2048 {input.align} > {output.align}
        """

rule downsample_prop:
    resources:
        nodes = 1,
        threads = 8
    shell:
        """
        module load samtools
        samtools view {input.align} \
            -h \
            -@ {resources.threads} \
            --subsample {params.prop} \
            --subsample-seed 42 > {output.align}
        """

rule dedupe_reads:
    resources:
        nodes = 2,
        threads = 16
    run:
        tiebreak_supp_reads(input.align, resources.threads, output.align)

rule get_dupe_reads:
    resources:
        nodes = 2,
        threads = 4
    run:
        get_dupe_read_names(input.align, resources.threads, output.txt)

# get flag

# get mapq values
rule get_mapq:
    resources:
        nodes = 2,
        threads = 4
    run:
        df = get_mapq(input.align, resources.threads)
        df.to_csv(output.fname, sep='\t', index=False)

rule filt_mapq:
    resources:
        nodes = 2,
        threads = 8
    shell:
        """
        module load samtools
        samtools view -q {params.mapq} -b {input.align} -o {output.align}
        """

rule bedcov:
    resources:
        threads = 8,
        nodes = 1
    shell:
        """
        module load samtools
        samtools bedcov {input.bed} {input.align} > {output.txt}
        """

# get the reads that map into regions from a bed file
# output in bam format
rule get_bed_reads:
    resources:
        threads = 16,
        nodes = 2
    shell:
        """
        module load samtools
        samtools view \
            -b \
            --region-file {input.bed} \
            {input.bam} > {output.bam}
        """

# get the reads matching specific read ids
rule get_bam_from_read_ids:
    resources:
        threads = 8,
        nodes = 1
    shell:
        """
        module load samtools
        samtools view \
            -h \
            -b {input.align} \
            --qname-file {input.read_ids} \
            -o {output.align}
        """

# fwd strand orient all reads
rule flip_reads:
    resources:
        threads = 8,
        nodes = 2
    run:
        import pysam
        reverse_strand = {0: 16, 16: 0}
        with pysam.AlignmentFile(input.bam, "rb",
                threads=resources.threads) as input_bam:
            with pysam.AlignmentFile(output.bam, "wb",
                                     template=input_bam,
                                     threads=resources.threads) as output_bam:
                for read in input_bam:
                    if read.has_tag('ts') and read.flag in reverse_strand:
                        if read.get_tag('ts') == '-':
                            read.flag = reverse_strand[read.flag]
                            read.set_tag('ts', '+')
                        output_bam.write(read)

rule bam_to_tss_bed:
    resources:
        threads = 1,
        nodes = 1
    shell:
        """
        module load bedtools
        bedtools bamtobed -i {input.bam} | awk '{{
            if ($6 == "+") {{
                print $1"\\t"$2"\\t"($2+1)"\\t"$4"\\t"$5"\\t"$6
            }} else {{
                print $1"\\t"($3-1)"\\t"$3"\\t"$4"\\t"$5"\\t"$6
            }}
        }}' > {output.bed}
        """
