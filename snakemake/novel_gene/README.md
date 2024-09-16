
With the file mapping transcript to gene id, add gene ids
```bash
conda activate cerberus
```
```python
import pyranges as pr
import pandas as pd
import cerberus

gtf = '../../data/transcripts_novel_gene_loci.gtf'
tsv = '/gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/04_transcriptome_assembly/04_evaluation/02_sqanti/data/240909merge/240909merge_associatedgene2isoform_noambigousISM_FSM_genic.tsv'

df = pd.read_csv(tsv, sep='\t', header=None)
df.columns = ['tid', 'gid']
gtf_df = pr.read_gtf(gtf).df

# novel genes separate
nov_gene = gtf_df.loc[gtf_df.gene_id.notnull()].copy(deep=True)
nov_gene['gene_name'] = nov_gene['gene_id']
gtf_df = gtf_df.loc[gtf_df.gene_id.isnull()].copy(deep=True)

# filter based on whether transcript is even in the dataset
gtf_df = gtf_df.loc[gtf_df.transcript_id.isin(df.tid.tolist())]

# get dict mapping tid:gid
tg_dict = dict([(entry['tid'], entry['gid']) for ind, entry in df.iterrows()])
tg_dict[list(tg_dict.keys())[0]]

# add gene id just for known genes
gtf_df['gene_id'] = gtf_df.transcript_id.map(tg_dict)
gtf_df['gene_name'] = gtf_df['gene_id']
assert len(gtf_df.loc[gtf_df.gene_id.isnull()].index) == 0

# cat the novel and known thing toegehter
gtf_df = pd.concat([gtf_df, nov_gene], axis=0)
assert len(gtf_df.loc[gtf_df.gene_id.isnull()].index) == 0

# add gene entries
df = gtf_df.copy(deep=True)
l1 = len(df.gene_id.unique().tolist())
g_df = make_hier_entry(df, how='g')
g_df['Source'] = 'ChatGPT'
g_df['Frame'] = '.'
g_df['Score'] = '.'
l2 = len(g_df.loc[g_df.Feature=='gene'].index)
assert l1 == l2
df = pd.concat([df, g_df], axis=0)
df = cerberus.sort_gtf(df)



# save
out_gtf = '../../data/novel_gene/transcripts_novel_gene_loci_filt_gene_name.gtf'
gtf_df = pr.PyRanges(gtf_df)
gtf_df.to_gtf(out_gtf)
