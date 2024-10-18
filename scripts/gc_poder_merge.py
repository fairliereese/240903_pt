import sys
import pyranges as pr
import cerberus
import pandas as pd

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

def merge_gtfs(gc_gtf_file,
               poder_gtf_file,
               output_file):

    gc_df = pr.read_gtf(gc_gtf_file)
    p_df = pr.read_gtf(poder_gtf_file)

    # get list of ics from each
    gc_ics = cerberus.get_ic(gc_df)
    p_ics = cerberus.get_ic(p_df)

    # get uniq ICs (gene+chrom+ic)
    # to poder that are not in gc
    gc_ics['id'] = gc_ics['Chromosome'].astype(str)+'_'+\
                   gc_ics['gene_id'].astype(str)+'_'+\
                   gc_ics['ic'].astype(str)
    p_ics['id'] = p_ics['Chromosome'].astype(str)+'_'+\
                   p_ics['gene_id'].astype(str)+'_'+\
                   p_ics['ic'].astype(str)
    p_ics['in_gc'] = p_ics['id'].isin(gc_ics['id'].tolist())
    p_ics[['in_gc', 'transcript_id']].groupby('in_gc').count()
    nov_tids = p_ics.loc[~p_ics['id'].isin(gc_ics['id'])].transcript_id.tolist()

    # filter poder gtf based on just novel transcripts
    p_df = p_df.df
    p_df = p_df.loc[p_df.transcript_id.isin(nov_tids)]
    assert len(p_df.loc[p_df.transcript_id.notnull()].transcript_id.unique())==len(nov_tids)

    # restrict to just transcript + exon entries
    gc_df = gc_df.df
    gc_df = gc_df.loc[gc_df.Feature.isin(['transcript', 'exon'])]
    p_df = p_df.loc[p_df.Feature.isin(['transcript', 'exon'])]

    # concatenate different dfs
    df = pd.concat([gc_df, p_df], axis=0)

    # drop gene name entries cause all they do is screw everything up
    df.drop('gene_name', axis=1, inplace=True)

        # # make sure that we don't have duplicates from gene name / gene id combos
    # l1 = len(df[['gene_name', 'gene_id']].drop_duplicates())
    # l2 = len(df.gene_id.unique())
    # print(l1)
    # print(l2)

        # temp = df[['gene_name', 'gene_id']].drop_duplicates()
    # temp.loc[temp.gene_id.duplicated(keep=False)].sort_values(by='gene_id').head()

    # make new gene entries for everything
    l1 = len(df.gene_id.unique().tolist())
    # make gene entry
    g_df = make_hier_entry(df, how='g')

    g_df['Source'] = 'v47_poder'
    g_df['Frame'] = '.'
    g_df['Score'] = '.'
    l2 = len(g_df.loc[g_df.Feature=='gene'].index)
    assert l1 == l2

    # concat them and then sort gtf
    df = pd.concat([df, g_df], axis=0)
    df = cerberus.sort_gtf(df)

    # make sure we've added the same number of transcripts
    n_gc_t = len(gc_df.loc[gc_df.transcript_id.notnull()].transcript_id.unique())
    n_nov_t = len(nov_tids)
    n_gc_p_t = len(df.loc[df.transcript_id.notnull()].transcript_id.unique())

    print(n_gc_t)
    print(n_nov_t)
    print(n_gc_p_t)
    assert n_gc_t+n_nov_t == n_gc_p_t

    df = pr.PyRanges(df)
    df.to_gtf(output_file)

def main():
    if len(sys.argv) != 4:
        print("Usage: python gc_poder_merge.py <gencode gtf> <poder gtf> <output_gtf_file>")
        sys.exit(1)

    gc_gtf_file = sys.argv[1]
    poder_gtf_file = sys.argv[2]
    output_file = sys.argv[3]

    merge_gtfs(gc_gtf_file,
               poder_gtf_file,
               output_file)

if __name__ == "__main__":
    main()
