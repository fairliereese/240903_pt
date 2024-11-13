Run whats in the ipynb
```bash
# gffread
gffread \
  -M \
  -T \
  test_gtf.gtf > test_merge.gtf

# gffread
gffread \
  -MTQ \
  test_gtf.gtf > test_merge_qd.gtf

# buildLoci
module load bedtools
awk '$3 != "gene"' test_gtf.gtf > test_gtf_no_gene.gtf # no genes # only exons

bedtools intersect \
    -s \
    -wao \
    -a test_gtf_no_gene.gtf \
    -b test_gtf_no_gene.gtf | \
    /gpfs/home/bsc/bsc083001/mele_lab/bin/buildLoci/buildLoci.pl - \
        locPrefix test_build_loci \
        > test_build_loci.gtf
```
