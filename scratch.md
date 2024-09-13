/gpfs/projects/bsc83/utils/conda_envs/transdecoder/opt/transdecoder/PerlLib/Pipeliner.pm

```bash
cd /gpfs/projects/bsc83/Projects/pantranscriptome/fairlie/240903_pt/snakemake/protein

orfanage             --cleanq             --mode LONGEST_MATCH             --reference /gpfs/projects/bsc83/Data/assemblies/GRCh38/GRCh38.primary_assembly.genome.fa             --query ../../data/transcripts_no_spike.gtf             --output ../../data/orfanage/transcripts.cds.gtf             --threads 2             /gpfs/projects/bsc83/Data/gene_annotations/gencode/v47/gencode.v47.primary_assembly.annotation.gtf

module load gcc
module load R/4.3.0
cpat             -x ../../ref/hexamer.tsv             -d ../../ref/cpat/cpat_model.RData             -g ../../data/orfanage/transcripts.cds_for_cpat.fa             --min-orf=100             --top-orf=50             -o ../../data/cpat/transcripts.no             1> ../../data/cpat/transcripts.no_cpat.output             2> ../../data/cpat/transcripts.no_cpat.error
  # \
  # 1> ../../data/cpat/transcripts.no_cpat.output \
  # 2> ../../data/cpat/transcripts.no_cpat.error

  python /gpfs/projects/bsc83/Projects/pantranscriptome/fairlie/240903_pt/scripts//refmt_gtf.py             ../../data/novel_gene/novel_genes_built.gtf             ../../data/novel_gene/novel_genes_fmt.gtf

  python /gpfs/projects/bsc83/Projects/pantranscriptome/fairlie/240903_pt/scripts//create_cpat_CDS_coordinates.py                                 --name ../../data/protein/transcripts                                 --sample_gtf ../../data/protein/transcripts_cpat_cds_to_be_predicted.gtf --called_orfs ../../data/protein/transcripts.ORF_remaining.tsv

  python /gpfs/projects/bsc83/Projects/pantranscriptome/fairlie/240903_pt/scripts//create_orf_table_for_sqanti_protein.py                     --trans
cript_exons_path ../../data/protein/transcripts.transcript_exons_only.gtf                     --cds_only_path ../../data/protein/transcripts.cds_renamed_exon.gtf                     --output_prefix transcripts
