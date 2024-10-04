Genomes all linked here: https://github.com/human-pangenomics/HPP_Year1_Assemblies/blob/main/assembly_index/Year1_assemblies_v2_genbank.index


```python
# make an ez parse table for fabien
import pandas as pd
import os
import sys

p = os.path.dirname(os.path.dirname(os.getcwd()))+'/scripts/'
sys.path.append(p)

from sm_utils import *
from utils import *

meta_file = '../config.tsv'
meta_file_2 = 'config.tsv'
genomes_file = 'genomes_config.tsv'

df = parse_config(meta_file)
df2 = pd.read_csv(meta_file_2, sep='\t')
df2['tech_rep'] = df2.cell_line_id+'_1'

# TODO test
# df2 = df2.loc[df2.cell_line_id == 'GM24385']

# get the genomes to download
g_df = pd.read_csv(genomes_file, sep='\t')

# maternal haplotypes
g_df['aws_mat_link'] = g_df['hap2_aws_fasta']
g_df = g_df.loc[g_df['aws_mat_link'].notnull()]
assert len(g_df.loc[g_df['aws_mat_link'].str.contains('maternal')].index) == len(g_df.index)

# paternal haplotypes
g_df['aws_pat_link'] = g_df['hap1_aws_fasta']
g_df = g_df.loc[g_df['aws_pat_link'].notnull()]
assert len(g_df.loc[g_df['aws_pat_link'].str.contains('paternal')].index) == len(g_df.index)


genome_cols = ['same_population_sample', 'european_sample',	'afr_sample']
g_df = g_df.loc[(g_df['sample'].isin(df2[genome_cols[0]]))|
                (g_df['sample'].isin(df2[genome_cols[1]]))|
                (g_df['sample'].isin(df2[genome_cols[2]]))]

# a little more df2 formatting
df2 = df2[['tech_rep', 'same_population_sample',
           'european_sample', 'afr_sample']].melt(id_vars='tech_rep')
df2 = df2.reset_index()
df2 = df2.rename({'variable':'assembly_status',
                  'value': 'assembly_sample'},
                  axis=1)

# limit just to the samples where we'll do this
df = df.loc[df.tech_rep.isin(df2.tech_rep.tolist())]

# get a key for assembly status, assembly sample, and actual sample
df2['dataset_key'] = df2.assembly_status+'_'+\
                     df2.assembly_sample+'_'+\
                     df2.tech_rep

```
