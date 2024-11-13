```bash
grep "^>" ../../ref/hg38.fa | sed 's/>//g' | awk '{print $1}' | sort > chromosomes.txt

awk '{print $1}' ../../ref/annot.gtf | sort | uniq > annot_chromosomes.txt

awk '{print $1}' ../../data/enc/enc.gtf | sort | uniq > enc_chromosomes.txt

```

```python
import pyranges as pr
import pandas as pd

chr_df = pd.read_csv('rename_chrom.tsv', sep='\t')
chr_map = dict(zip(chr_df['chr'], chr_df['new_chr']))


df = pr.read_gtf('../../data/enc/enc.gtf').df
df['Chromosome'] = df['Chromosome'].map(chr_map)
df = pr.PyRanges(df)
df.to_gtf('enc_renamed_chr.gtf')
```
