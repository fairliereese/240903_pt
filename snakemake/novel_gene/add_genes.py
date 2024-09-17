import pyranges as pr
import pandas as pd
import cerberus

def make_hier_entry(df, how='t'):
    """
    kind {'g','t'}
    """
    agg_dict = {'min_coord': 'min', 'max_coord': 'max'}
    t_df = df.copy(deep=True)
    t_df['min_coord'] = t_df[['Start', 'End']].min(axis=1)
    t_df['max_coord'] = t_df[['Start', 'End']].max(axis=1)
    if how == 't':
        gb_cols = ['Chromosome', 'Strand', 'gene_name',
                   'gene_id', 'transcript_id', 'transcript_name',
                   'tss_id', 'tes_id',
                   'new_transcript_id', 'original_transcript_id',
                   'original_transcript_name', 'ag1', 'ag2']
    elif how == 'g':
        gb_cols = ['Chromosome', 'Strand', 'gene_name',
                   'gene_id']
    gb_cols = list(set(gb_cols)&(set(t_df.columns)))

    cols = gb_cols + ['min_coord', 'max_coord']
    t_df = t_df[cols]
    t_df = t_df.groupby(gb_cols, observed=True).agg(agg_dict).reset_index()
    t_df.rename({'min_coord': 'Start', 'max_coord': 'End'}, axis=1, inplace=True)
    if how == 't':
        t_df['Feature'] = 'transcript'
    elif how == 'g':
        t_df['Feature'] = 'gene'

    return t_df

nov_gtf = '../../data/transcripts_novel_gene_loci.gtf'
tsv = '../../data/240909merge_associatedgene2isoform_noambigousISM_FSM_genic.tsv'
known_gtf = '/gpfs/projects/bsc83/Projects/pantranscriptome/novelannotations/merged/240909_merged_withoutchrEBV.gtf'

df = pd.read_csv(tsv, sep='\t', header=None)
df.columns = ['tid', 'gid']
# remove "novelGenes" from sqanti mappings
df = pd.read_csv(tsv, sep='\t', header=None)
df.columns = ['tid', 'gid']
l1 = len(df.index)
df = df.loc[~df.gid.str.startswith('novelGene')]
l2 = len(df.index)
assert l1 != l2

nov_gtf_df = pr.read_gtf(nov_gtf).df

# novel genes separate, remove gene entries
nov_gene = nov_gtf_df.loc[nov_gtf_df.gene_id.notnull()].copy(deep=True)
nov_gene['gene_name'] = nov_gene['gene_id']
nov_gene['Source'] = 'ChatGPT'
nov_gene = nov_gene.loc[nov_gene.Feature!='gene']

# known genes
gtf_df = pr.read_gtf(known_gtf).df
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

# cat the novel and known thing together
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

print(len(gtf_df.index))
df = pd.concat([gtf_df, g_df], axis=0)
print(len(df.index))

df = cerberus.sort_gtf(df)


assert len(df.loc[df.Feature=='gene']) == len(df.gene_id.unique())


# save
out_gtf = '../../data/novel_gene/transcripts_novel_gene_loci_filt_gene_name.gtf'
df = pr.PyRanges(df)
df.to_gtf(out_gtf)
