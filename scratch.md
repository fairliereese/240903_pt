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

  python /gpfs/projects/bsc83/Projects/pantranscriptome/fairlie/240903_pt/scripts//create_orf_table_for_sqanti_protein.py                     --transcript_exons_path ../../data/protein/transcripts.transcript_exons_only.gtf                     --cds_only_path ../../data/protein/transcripts.cds_renamed_exon.gtf                     --output_prefix transcripts

conda activate /gpfs/projects/bsc83/utils/conda_envs/SQANTI3-5.2.1
  python /gpfs/projects/bsc83/Projects/pantranscriptome/fairlie/240903_pt/scripts//sqanti3_protein.py                     ../../data/protein/transcripts.transcript_exons_only.gtf                     ../../data/protein/transcripts.cds_renamed_exon.gtf                     ../../data/protein/transcripts_best_orf.tsv                     ../../data/protein/transcripts_gencode.transcript_exons_only.gtf                     ../../data/protein/transcripts_gencode.cds_renamed_exon.gtf                     -d ./                     -p ../../data/protein/transcripts.sqanti


  python /gpfs/projects/bsc83/Projects/pantranscriptome/fairlie/240903_pt/scripts//create_protein_overview_table.py             --best_orf_path ../../data/protein/transcripts_best_orf.tsv             --sqanti_protein_path ../../data/protein/transcripts.sqanti_protein_classification.tsv             --orf_completeness_path ../../data/protein/transcripts_ORF_completeness.tsv             --output_name ../../data/protein/transcripts_protein_annotation.tsv --gtf_original_path ../../data/transcripts_no_spike.gtf             --gtf_predicted_path ../../data/protein/transcripts_protein.gtf             --protein_fasta_path ../../data/protein/transcripts_protein.fa             --blastp_path ../../data/protein/transcripts_blastp.out


python /gpfs/projects/bsc83/Projects/pantranscriptome/fairlie/240903_pt/scripts/filter_cpat.py \
          --input_file_path ../../data/protein/transcripts.ORF_prob.tsv \
          --orf_input_seq_path ../../data/protein/transcripts.ORF_seqs.fa \
          --output_path ../../data/protein/transcripts.ORF_remaining.tsv \
          --first_cutoff 0.725 \
          --second_cutoff 0.364

python /gpfs/projects/bsc83/Projects/pantranscriptome/fairlie/240903_pt/scripts//recover_source.py \
  --combined_cds_path ../../data/protein/transcripts_protein_unsourced.gtf \
  --source_gtf_path ../../data/transcripts_no_spike_no_ebv.gtf \
  --cpat_cds_path ../../data/protein/transcripts_cpat_with_cds.gtf \
  --orfanage_cds_path ../../data/protein/transcripts_orfanage_cds_filtered_stop_codon_corrected.gtf \
  --output_path ../../data/protein/transcripts_protein.gtf

  python /gpfs/projects/bsc83/Projects/pantranscriptome/fairlie/240903_pt/scripts//create_protein_overview_table.py \
               --best_orf_path ../../data/protein/transcripts_best_orf.tsv \
               --sqanti_protein_path ../../data/protein/transcripts.sqanti_protein_classification.tsv \
               --orf_completeness_path ../../data/protein/transcripts_ORF_completeness.tsv \
               --output_name ../../data/protein/transcripts_protein_annotation.tsv \
               --gtf_original_path ../../data/transcripts_no_spike_no_ebv.gtf \
               --gtf_predicted_path ../../data/protein/transcripts_protein.gtf \
               --protein_fasta_path ../../data/protein/transcripts_protein.fa \
               --blastp_path ../../data/protein/transcripts_blastp.out

# to run locally
 python /Users/fairliereese/Documents/programming/mele_lab/projects/240903_pt/scripts/create_protein_overview_table.py \
              --best_orf_path ../../data/protein/transcripts_best_orf.tsv \
              --sqanti_protein_path ../../data/protein/transcripts.sqanti_protein_classification.tsv \
              --orf_completeness_path test \
              --output_name ../../data/protein/transcripts_protein_annotation.tsv \
              --gtf_original_path ../../data/transcripts_no_spike_no_ebv.gtf \
              --gtf_predicted_path ../../data/protein/transcripts_protein.gtf \
              --protein_fasta_path ../../data/protein/transcripts_protein.fa \
              --blastp_path ../../data/protein/transcripts_blastp.out

```bash
cd ~/mele_lab/bin/kallisto/build
rm -r *
cmake .. -DBUILD_FUNCTESTING=ON -DENABLE_AVX2=OFF -DCOMPLIATION_ARCH=OFF -DMAX_KMER_SIZE=64
make
make install
```

```bash
python /Users/fairliereese/Documents/programming/mele_lab/projects/240903_pt/scripts/check_orf_completeness.py \
  --cpat_seqs ../../data/protein/transcripts.ORF_seqs.fa \
  --orfanage_seqs ../../data/protein/transcripts_orfanage_orfs.fa \
  --cpat_info ../../data/protein/transcripts.ORF_remaining.tsv \
  --orfanage_info ../../data/protein/transcripts_orfanage_cds_filtered_stop_codon_corrected.gtf \
  --output_path test
```

```python
!conda activate /gpfs/home/bsc/bsc083001/miniconda3/envs/lr-kallisto
import ngs_tools as ngs
gtf_file = '/gpfs/projects/bsc83/Projects/pantranscriptome/novelannotations/merged/240926_filtered_with_genes.gtf'
gtf = ngs.gtf.genes_and_transcripts_from_gtf(gtf_file)

```

<!-- ```python
r'''
 29         ^(?P<chromosome>.+?)\s+ # chromosome

 30         .*?\t                   # source
 31         (?P<feature>.+?)\s+     # feature: transcript, exon, etc.
 32         (?P<start>[0-9]+?)\s+   # start position (1-indexed)
 33         (?P<end>[0-9]+?)\s+     # end position (1-indexed, inclusive)
 34         .*?\s+                  # score
 35         (?P<strand>\+|-|\.)\s+  # +, -, . indicating strand
 36         .*?\s+                  # frame
 37         (?P<attributes>.*)      # attributes
``` -->

```bash
module load intel/2023.0
module load samtools
fa=/gpfs/projects/bsc83/Projects/pantranscriptome/fairlie/240903_pt/ref/HG002/HG002.fa.gz
samtools faidx $fa chrX > test
```
