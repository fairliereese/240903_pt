```bash
grep "^>" ../../ref/hg38.fa | sed 's/>//g' | awk '{print $1}' | sort > chromosomes.txt

awk '{print $1}' ../../ref/annot.gtf | sort | uniq > annot_chromosomes.txt

awk '{print $1}' ../../ref/annot.gtf | sort | uniq > enc_chromosomes.txt

```
