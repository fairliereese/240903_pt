#!/bin/bash

module load R/4.3.2
module load vcftools/0.1.16
module load bcftools/1.19

script_data="/gpfs/projects/bsc83/Projects/scRNAseq/imestres/sQTLseeker/sqtlseeker2-nf/scripts/prepare_lineage_data.R"
script_genotype="/gpfs/projects/bsc83/Projects/scRNAseq/imestres/sQTLseeker/sqtlseeker2-nf/scripts/prepare_lineage_genotype.R"

pca_dir="/gpfs/projects/bsc83/Projects/scRNAseq/imestres/sQTLseeker/analyses/pca/all_lineages/lineage_EUB-AFB.eigenvec"
metadata_tsv="/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/00_MetaData/sample_metadata_expanded.tsv"

tcc_exp_dir="/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/02_Processing/04_pseudobulk/03_split-group-by-Condition/01_split_Condition-group_POP/01_tcc/01_lineage/"
tcc_suffix="-kb-tcc-pseudobulk-avg-reps-split_Condition-group_POP-filt.rds"

tcc_gene_ids_dir="/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/03_Differential/02_DTU/02_DTU/02_popDTU/01_popDTU/00_Archive-never-rm/01_Run-Age-POP/"
ids_suffix="-tcc-avg-reps-popDTU.rds"

annotation_gtf="/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/01_Randolph2021/01_Quantification/00_IndexData/01_GENCODEv42/gencode.v42.annotation.gtf"

out_dir="/gpfs/projects/bsc83/Projects/scRNAseq/imestres/sQTLseeker/sqtlseeker2-nf/data/lineages/"

geno_dir="${out_dir}genotype_012"
geno_vcf="/gpfs/projects/bsc83/Projects/scRNAseq/imestres/sQTLseeker/runs/data/merged_genotypes.vcf.gz"

lineage="MONO,B,NK,T.CD4,T.CD8"


echo "Creating input files...."
Rscript $script_data \
	--gtf $annotation_gtf \
	--metadata $metadata_tsv \
	--tcc-dir $tcc_exp_dir \
	--suffix=$tcc_suffix \
	--id-dir $tcc_gene_ids_dir \
	--id-rds=$ids_suffix \
	--pca-dir $pca_dir \
	--out-dir $out_dir \
    	--lineage=$lineage \
	--verbose

echo "Processing merged genotype...."
mkdir -p "$geno_dir"

vcftools --gzvcf "$geno_vcf" --012 --out "$geno_dir/matrix"

matrix_012="$geno_dir/matrix.012"
pos_012="$geno_dir/matrix.012.pos"
indv_012="$geno_dir/matrix.012.indv"

echo "Building up genotype from matrices...."
Rscript $script_genotype \
    --genotype $geno_vcf \
    --metadata "${out_dir}lineage_metadata-EUB.AFB.tsv" \
    --in-012 $matrix_012 \
    --indv-012 $indv_012 \
    --pos-012 $pos_012 \
    --out-dir $out_dir \
    --verbose

echo "Done!"
meta=/gpfs/projects/bsc83/Projects/scRNAseq/imestres/sQTLseeker/sqtlseeker2-nf/data/lineages/lineage_metadata-EUB.AFB.tsv
