import cerberus
import pandas as pd
import pyranges as pr
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Process GTF and SJ files.')
parser.add_argument('--cl_file', required=True, help='Path to classification file')
parser.add_argument('--gtf_file', required=True, help='Path to GTF file')
parser.add_argument('--sjs_file', required=True, help='Path to splice junction file')
parser.add_argument('--sjs_summary', required=True, help='Path to output splice junction summary')
parser.add_argument('--ic_summary', required=True, help='Path to output IC summary')
parser.add_argument('--cell_line_id', required=True, help='Cell line ID')
parser.add_argument('--map_genome', required=True, help='Mapping genome')
parser.add_argument('--sqanti_genome', required=True, help='SQANTI genome')
args = parser.parse_args()

######################
## IC DF
#######################
# get ics + tid from cerberus
df = pr.read_gtf(args.gtf_file, duplicate_attr=True)
df = df.df
df = cerberus.add_stable_gid(df)
df = pr.PyRanges(df)
df = cerberus.get_ic(df)

# remove monoexonic
df = df.loc[df.ic != '-']

# make sure no monoexonic
assert len(df.loc[df.ic == '-']) == 0

# remove gene id; we're not considering this here
df.drop('gene_id', axis=1, inplace=True)

# get ic_id from chr strand and ic; then drop constituent cols
df['ic_id'] = df.Chromosome.astype(str) + '_' +\
              df.Strand.astype(str) + '_' +\
              df.ic.astype(str)
df.drop(['Chromosome', 'Strand', 'ic'], axis=1, inplace=True)

# add the classification from the sqanti file
sq_df = pd.read_csv(args.cl_file, sep='\t')
sq_df = sq_df[['isoform', 'structural_category']].rename({'isoform': 'transcript_id'}, axis=1)
df = df.merge(sq_df, how='left', on='transcript_id')

# gb ic and structural cat; record tids from each
# make sure we're only annotating 1 structural category per ic
df = df.groupby(['ic_id', 'structural_category']).agg(','.join).reset_index()
assert len(df.loc[df.ic_id.duplicated(keep=False)].index) == 0

# add sample, map/td genome, sqanti genome metadata and concat
df['cell_line_id'] = args.cell_line_id
df['map_genome'] = args.map_genome
df['sqanti_genome'] = args.sqanti_genome

df.to_csv(args.ic_summary, sep='\t', index=False)

######################
## SJ DF
#######################
df = pd.read_csv(args.sjs_file, sep='\t')
df = df[['chrom', 'strand', 'genomic_start_coord', 'genomic_end_coord', 'splice_site', 'canonical']]
df['sj_id'] = df['chrom'].astype(str) + '_' +\
              df['strand'].astype(str) + '_' +\
              df['genomic_start_coord'].astype(str) + '_' +\
              df['genomic_end_coord'].astype(str)
df = df[['sj_id', 'splice_site', 'canonical']].drop_duplicates()
df = df.rename({'splice_site': 'splice_motif'}, axis=1)

# add sample, map/td genome, sqanti genome metadata and concat
df['cell_line_id'] = args.cell_line_id
df['map_genome'] = args.map_genome
df['sqanti_genome'] = args.sqanti_genome

df.to_csv(args.sjs_summary, sep='\t', index=False)
