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
vcf=/Users/fairliereese/Documents/programming/mele_lab/projects/240903_pt/data/td_personal/exonic_snps/HG00621_wha_happen.vcf

sj_10nt_bed=/Users/fairliereese/Documents/programming/mele_lab/projects/240903_pt/data/td_personal/sqanti/sqanti/sj_10nt.bed
sj_12nt_bed=/Users/fairliereese/Documents/programming/mele_lab/projects/240903_pt/data/td_personal/sqanti/sj_ss_12nt.bed
ss_2nt_bed=/Users/fairliereese/Documents/programming/mele_lab/projects/240903_pt/data/td_personal/sqanti/ss_2nt.bed

sj_10nt_bed_int=/Users/fairliereese/Documents/programming/mele_lab/projects/240903_pt/data/td_personal/sqanti/sqanti/sj_10nt_int.txt
sj_12nt_bed_int=/Users/fairliereese/Documents/programming/mele_lab/projects/240903_pt/data/td_personal/sqanti/sj_ss_12nt_int.txt
ss_2nt_bed_int=/Users/fairliereese/Documents/programming/mele_lab/projects/240903_pt/data/td_personal/sqanti/ss_2nt_int.txt

bedtools intersect \
         -a $sj_10nt_bed \
         -b $vcf \
         -wa -wb > $sj_10nt_bed_int

```
