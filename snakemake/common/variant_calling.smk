rule ref_make_gatk_dict:
    resources:
        threads = 1,
        nodes = 1
    shell:
        """
        module load gatk
        gatk-launch CreateSequenceDictionary -R {input.fa}
        """

rule ref_make_fa_ind:
    resources:
        threads = 1,
        nodes = 1
    shell:
        """
        module load samtools
        samtools faidx {input.fa}
        """

rule call_variants_gatk:
    resources:
        threads = 112,
        nodes = 6
    shell:
        """
        module purge
        module load java-openjdk/17.0.11+9
        module load gatk/4.5.0.0
        export OMP_NUM_THREADS={resources.threads}
        gatk \
            --java-options "-Xmx36g -Xms32g" \
            HaplotypeCaller \
            -I {input.bam} \
            -R {input.fa} \
            -O {output.vcf} \
            -ERC BP_RESOLUTION \
            --sample-ploidy {params.ploidy} \
            --max-reads-per-alignment-start 0 \
            --dont-use-soft-clipped-bases true \
            --native-pair-hmm-threads {resources.threads} \
            --max-alternate-alleles 6 \
            --max-genotype-count 1024 \
            --bamout {output.bam} \
            --create-output-bam-index \
            --sample-name {wildcards.dataset}
        """

# rule filt_variants:
#     params:
#         qual_thresh = 1000
#     resources:
#         threads = 1,
#         nodes = 1
#     run:
#         df = read_vcf(input.vcf)
#
#         # remove non-ref alleles
#         l1 = len(df.index)
#         df = df.loc[df.ALT!='<NON_REF>']
#         l2 = len(df.index)
#         # assert l1 != l2 # don't need this anymore as we discovered a setting that turns it off
#
#         # remove reads w/ qual under thresh
#         df.QUAL = df.QUAL.astype('float')
#         df = df.loc[df.QUAL>=params.qual_thresh]
#         write_vcf(df, output.vcf, input.vcf)

rule intersect_variants_with_bed:
    resources:
        threads = 1,
        nodes = 1
    shell:
        """
        module load bedtools
        bedtools intersect -a {input.vcf} -b {input.bed} > {output.bed}
        """

rule rev_intersect_variants_with_bed:
    resources:
        threads = 1,
        nodes = 1
    shell:
        """
        module load bedtools
        bedtools intersect -v -a {input.vcf} -b {input.bed} > {output.bed}
        """

rule merge_variants:
    resources:
        threads = 1,
        nodes = 1
    conda:
        'base'
    shell:
        """
        bcftools merge \
            --threads {resources.threads} \
            --use-header {input.header_vcf} \
            {input.vcfs} > {output.vcf}
        """

rule bgzip:
    resources:
        threads = 1,
        nodes = 1
    conda:
        'htslib'
    shell:
        """
        bgzip -c {input.ifile} > {output.gz}
        """

rule bcftools_vcf_index:
    resources:
        threads = 16,
        nodes = 2
    conda:
        'base'
    shell:
        """
        bcftools index \
            {input.vcf} \
            --threads {resources.threads}
        """


rule vcf_index:
    resources:
        threads = 1,
        nodes = 1
    conda:
        'bcftools'
    shell:
        """
        tabix -p vcf {input.vcf}
        """

# remove PL tag because it doesn't work with bcftools
# at higher ploidies
rule vcf_rm_PL:
    resources:
        threads = 1,
        nodes = 1
    conda:
        'base'
    shell:
        """
        if [ -s "{input.vcf}" ]; then
            bcftools annotate \
                -x ^INFO/DP,^FORMAT/GT,^FORMAT/AD,^FORMAT/DP,^FORMAT/GQ \
                {input.vcf} > {output.vcf}
        else
            touch {output.vcf}
        fi
        """

rule vcf_norm:
    resources:
        threads = 1,
        nodes = 2
    conda:
        'base'
    shell:
        """
        if [ -s "{input.vcf}" ]; then
            bcftools norm \
              -a \
              -m -any \
              --fasta-ref {input.fa} \
              --old-rec-tag INFO \
              {input.vcf} > {output.vcf}
        else
            touch {output.vcf}
        fi
        """
# https://stackoverflow.com/questions/78118781/separating-alleles-for-bcftools-merge-when-ploidy-2

rule vcf_rm_tags:
    resources:
        threads = 1,
        nodes = 2
    conda:
        'base'
    shell:
        """
        bcftools annotate  \
            -x INFO \
            -O z \
            -o {output.vcf} {input.vcf}
        """

rule bcftools_subset_on_samples:
    resources:
        threads = 1,
        nodes = 1
    conda:
        'base'
    shell:
        """
        bcftools view \
        --samples {params.samples} \
        -Ov \
        {input.vcf} > {output.vcf}
        """
rule vcftools_calc_af:
    resources:
        threads = 1,
        nodes = 1
    shell:
        """
        vcftools \
            --vcf {input.vcf} \
            --freq \
            --out {params.opref}
        """

rule bcftools_concat:
    resources:
        threads = 8,
        nodes = 3
    conda:
        'base'
    shell:
        """
        # output type = v --> compressed vcf
        bcftools concat \
            --output-type v \
            --threads {resources.threads} \
            -o {output.vcf} \
            {params.cli_vcfs}
            # chr1.vcf chr2.vcf chr3.vcf ... chrX.vcf
        """

rule bcftools_filter_on_regions:
    resources:
        threads = 8,
        nodes = 2
    conda:
        'base'
    shell:
        """
        bcftools view \
            {input.vcf} \
            --threads {resources.threads} \
            --regions-file {input.bed} \
            --output-type v \
            --output {output.vcf}
        """

rule vcftools_012:
    resources:
        threads = 1,
        nodes = 2
    shell:
        """
        vcftools \
            --gzvcf {input.vcf} \
            --012 \
            --out {params.opref}
        """
