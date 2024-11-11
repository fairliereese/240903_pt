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
bedtools intersect \
    -s \
    -wao \
    -a test_gtf.gtf \
    -b test_gtf.gtf | \
    buildLoci.pl - \
        locPrefix test_build_loci \
        > test_build_loci.gtf
```
