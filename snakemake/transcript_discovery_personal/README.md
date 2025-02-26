## too lazy to do this w/ snakemake
```bash
head -20000 /gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/09_other_analyses/02_prepare_vep_inputs/1000G_variants_in_exons_headed.vcf |grep ^#  > ~/mele_lab/projec
ts/240903_pt/data/1000g/1000G_vcf_header.vcf
scp ~/mele_lab/projects/240903_pt/data/1000g/1000G_vcf_header.vcf .

cd ~/mele_lab/projects/240903_pt/data/td_personal/exonic_snps/

fname=HG00621_exon_ss_vars_int_loj.bed
grep 60792131 $fname > ${fname}_pos

fname=HG00621_exon_vars_int_loj.bed
grep 60792131 $fname > ${fname}_pos

fname=HG00621_ss_vars_int_loj.bed
grep 60792131 $fname > ${fname}_pos
```


I have no idea what's going on with these dumb intersections:
```bash
# vcf for one sampel encompassing region that seems to be missing an SNP
1000g_vcf=/gpfs/projects/bsc83/Data/1000G/x30_vcf/original_data/1kGP_high_coverage_Illumina.all_chr.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
vcf=~/mele_lab/projects/240903_pt/data/td_personal/exonic_snps/HG00621_wha_happen.vcf

module load bedtools
bcftools view \
    --samples HG00621 \
    --threads {resources.threads} \
    -Ov \
    --min-ac=1 \
    --regions chr10:60791471-60792700 \
    $1000g_vcf > $vcf

sj_10nt_bed=~/mele_lab/projects/240903_pt/data/td_personal/sqanti/sj_10nt.bed
sj_12nt_bed=~/mele_lab/projects/240903_pt/data/td_personal/sqanti/sj_ss_12nt.bed
ss_2nt_bed=~/mele_lab/projects/240903_pt/data/td_personal/sqanti/ss_2nt.bed

sj_10nt_bed_int=~/mele_lab/projects/240903_pt/data/td_personal/sqanti/sj_10nt_int.bed
sj_12nt_bed_int=~/mele_lab/projects/240903_pt/data/td_personal/sqanti/sj_ss_12nt_int.bed
ss_2nt_bed_int=~/mele_lab/projects/240903_pt/data/td_personal/sqanti/ss_2nt_int.bed

# # this workds fine; I get results
# bedtools intersect \
#          -a $sj_10nt_bed \
#          -b $vcf \
#          -wa -wb > $sj_10nt_bed_int


# this works fine too; I get results
module load bedtools
module load bcftools
bcftools view \
   --samples HG00621 \
   -Ov \
   --min-ac=1 \
   $vcf | \
bedtools intersect \
    -a $sj_10nt_bed \
    -b stdin \
    -wa \
    -u > $sj_10nt_bed_int




# this works fine too; I get results
module load bedtools
module load bcftools
bcftools view \
   --samples HG00621 \
   -Ov \
   --min-ac=1 \
   $vcf | \
bedtools intersect \
    -a $sj_10nt_bed \
    -b stdin \
    -wa \
    -u > $sj_10nt_bed_int

# try without header
module load bedtools
module load bcftools
bedtools intersect \
    -a $sj_10nt_bed \
    -b stdin \
    -wa \
    -u > $sj_10nt_bed_int



# this is not empty but does not contain the 2 variants that it should be overlapping that I get above
bcftools view \
    --samples HG00621 \
    -Ov \
    --min-ac=1 \
    $1000g_vcf | \
bedtools intersect \
     -a $sj_10nt_bed \
     -b stdin \
     -wa \
     -u > /gpfs/projects/bsc83/Projects/pantranscriptome/fairlie/240903_pt/data/td_personal/exonic_snps/HG00621_exon_vars_int.bed

```
