```bash
conda activate pt_snakemake
snakemake \
  -s Snakefile \
  -j 100 \
  --latency-wait 120 \
  --use-conda \
  --cluster \
    "sbatch \
    --nodes {resources.nodes} \
    -q gp_bscls \
    -A bsc83 \
    -c {resources.threads}  \
    --mail-user=freese@bsc.es \
    --mail-type=START,END,FAIL \
    --time=12:00:00" \
    -n

  snakemake \
    -s Snakefile \
    -j 100 \
    --latency-wait 120 \
    --use-conda \
    --cluster \
      "sbatch \
      --nodes {resources.nodes} \
      -q highmem \
      -A bsc83 \
      -c {resources.threads}  \
      --mail-user=freese@bsc.es \
      --mail-type=START,END,FAIL \
      --time=12:00:00" \
      -n

snakemake \
  -s Snakefile \
  -j 100 \
  --latency-wait 120 \
  --use-conda \
  -n

```
