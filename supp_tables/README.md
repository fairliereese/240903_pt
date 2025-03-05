
## Code to fix metadata table

```python
import pandas as pd
df = pd.read_csv('00_sample_metadata.tsv', sep='\t')
df = df.loc[df.merged_run_mode==True]
df = df.loc[df.mixed_samples == False]
assert len(df.index) == 43
df.to_csv('00_sample_metadata.tsv', sep='\t')

df = pd.read_csv('00_sample_metadata.tsv', sep='\t')
df.drop(['lab_number_sample',
         'wgs_sr_ill',
         'wgs_lr_pb',
         'wgs_lr_ont',
         'mixed_samples'], axis=1,
       inplace=True)
df.to_csv('00_sample_metadata.tsv', sep='\t', index=False)


```

## Add tau values to the master table
(Code stolen from `tau_vs_pop_specific.ipynb`)
```python
import pandas as pd
min_cpm = 1
min_samples = 1
f = f'../analysis/241108_{min_cpm}_{min_samples}_mean_tau.tsv'
print(f)
tau_df = pd.read_csv(f, sep='\t')

mt_df = pd.read_csv('04_poder_mt.tsv', sep='\t')
```

## Get tau values for Enh. GENCODE MAGE
From `tau_mage_enh_gencode.ipynb`.
```python
import pandas as pd
min_cpm = 1
min_samples = 1
tau_df = pd.read_csv(f'../analysis/250127_{min_cpm}_{min_samples}_mage_population_mean_tau.tsv', sep='\t')
tau_df.rename({'tid':'isoform'}, axis=1, inplace=True)
tau_df.to_csv('13_mage_tau.tsv', sep='\t', index=False)
```
