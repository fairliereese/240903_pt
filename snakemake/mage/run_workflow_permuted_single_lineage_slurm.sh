#!/bin/bash

#SBATCH --job-name=T.CD4_sQTLseeker_raw
#SBATCH --output=slurm_logs/T.CD4_raw_%j.out
#SBATCH --error=slurm_logs/T.CD4_raw_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --qos=gp_bscls
#SBATCH --time="09:00:00"
#SBATCH --account=bsc83

module load nextflow/21.04.0
module load singularity/3.11.5
module load R/4.3.2

lin="T.CD4"

image="/gpfs/projects/bsc83/Projects/scRNAseq/imestres/sQTLseeker/runs/sqtlseeker2-nf.sif"
prog="/gpfs/projects/bsc83/Projects/scRNAseq/imestres/sQTLseeker/sqtlseeker2-nf/sqtlseeker2.nf"
data_dir="/gpfs/projects/bsc83/Projects/scRNAseq/imestres/sQTLseeker/sqtlseeker2-nf/new_data/raw/lineages"
out_dir="/gpfs/projects/bsc83/Projects/scRNAseq/imestres/sQTLseeker/results/run_20241209/${lin}/"

f_genes="${data_dir}/lineage_genes-EUB.AFB.bed"
f_tcc="${data_dir}/${lin}_tcc-exp-EUB.AFB.tsv.gz"
f_metadata="${data_dir}/lineage_metadata-EUB.AFB.tsv"
f_genotype="${data_dir}/lineage_genotype-EUB.AFB.vcf.gz"

mkdir -p $out_dir
cd $out_dir

nextflow run $prog \
	--genotype $f_genotype \
	--trexp $f_tcc \
	--metadata $f_metadata \
	--genes $f_genes \
	--dir $out_dir \
	--mode "permuted" \
	--covariates "true" \
	--max_perm 1000 \
  --svqtl "true" \
  --fdr 0.01 \
  --fdr_svqtl 0.01 \
  --min_gene_exp 5 \
  --min_transcript_exp 1 \
  --min_dispersion 0.01 \
  --min_md 0.05 \
  --min_proportion 0.4 \
  --with-singularity $image
