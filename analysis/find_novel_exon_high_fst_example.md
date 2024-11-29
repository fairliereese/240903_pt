
```bash
gid='ENSG00000136146.15'
sample='YRI3'
sample_1kg='NA18906'

# get ref coords of gene from gtf
gtf='../data/transcripts_filt_with_gene.gtf'
grep "$gid" "$gtf" | awk '$3 == "gene" { print $1, $4 - 1, $5, $9, ".", $7 }' OFS="\t" > $gid.bed

module load vcftools
vcf='/gpfs/projects/bsc83/Data/1000G/x30_vcf/original_data/1kGP_high_coverage_Illumina.all_chr.filtered.SNV_INDEL_SV_phased_panel_biallelic.vcf.gz'
vcftools \
    --gzvcf $vcf \
    --bed $gid.bed \
    --indv $sample_1kg \
    --out ${gid}_1000g_var_in_exons.vcf \
    --recode \
    --recode-INFO-all


gid='ENSG00000205220.12'
sample='YRI6'
sample_1kg='NA19129'

# get ref coords of gene from gtf
gtf='../data/transcripts_filt_with_gene.gtf'
grep "$gid" "$gtf" | awk '$3 == "gene" { print $1, $4 - 1, $5, $9, ".", $7 }' OFS="\t" > $gid.bed

module load vcftools
vcf='/gpfs/projects/bsc83/Data/1000G/x30_vcf/original_data/1kGP_high_coverage_Illumina.all_chr.filtered.SNV_INDEL_SV_phased_panel_biallelic.vcf.gz'
vcftools \
    --gzvcf $vcf \
    --bed $gid.bed \
    --indv $sample_1kg \
    --out ${gid}_1000g_var_in_exons.vcf \
    --recode \
    --recode-INFO-all
```

# replace line thing w/ this `##INFO=<ID=END2,Number=1,Type=Integer,Description="Position of breakpoint on CHR2">`
