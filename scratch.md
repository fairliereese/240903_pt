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

```bash
module load samtools
bam=/gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/03_mapping/data/pantrx_general_mapping/genomic12_NI7_GM19240.bam
var=C
chrom=chr6
pos=29942581
samtools view -h ${bam} ${chrom}:${pos}-${pos} | \
        awk 'BEGIN {{pos = ${pos}; allele = "${var}"}} \
         $0 ~ /^@/ || (substr($10, pos - $4 + 1, 1) == allele) {{print $0}}' > ~/test.sam
```

```bash
conda activate /gpfs/home/bsc/bsc083001/miniconda3/envs/kallisto
/gpfs/home/bsc/bsc083001/miniconda3/envs/kallisto/bin/kallisto bus \
     -t 32 \
     -x bulk \
     -i ../../ref/v47_kallisto_short/hg38_v47_k-31.idx \
     -g ../../ref/v47_kallisto_short/hg38_v47.t2g \
     -o ../../data/mage/v47_kallisto/NA19704_batch11_rep1/ \
     --paired \
     ../../data/mage/raw/NA19704_batch11_rep1_r1.fq.gz \
     ../../data/mage/raw/NA19704_batch11_rep1_r2.fq.gz

conda activate /gpfs/home/bsc/bsc083001/miniconda3/envs/kallisto
kb count \
  -x BULK \
  -o ../../data/mage/v47_kallisto/NA19704_batch11_rep1/ \
              -i ../../ref/v47_kallisto_short/hg38_v47_k-31.idx             -g ../../ref/v47_kallisto_short/hg38_v47.
t2g             --parity=paired             --strand=unstranded             --tcc             --matrix-to-directories             ../../data/mage/raw/NA19704_batch11_rep1_r1.fq.gz ../../data/mage/raw/NA19704_batch11_re

p1_r2.fq.gz



         bus_sort: ../../data/mage/v47_kallisto/{sample}/output_sorted.bus

    count_pref: ../../data/mage/v47_kallisto/{sample}/count
    count_mtx: ../../data/mage/v47_kallisto/{sample}/count.mtx
    count_ec: ../../data/mage/v47_kallisto/{sample}/count.ec.txt
    matrix_tsv: ../../data/mage/v47_kallisto/{sample}/matrix.abundance.tsv

    quant:
      odir: ../../data/mage/v47_kallisto_quant/{sample}
      matrix: ../../data/mage/v47_kallisto_quant/{sample}/matrix.abundance.mtx
      matrix_tpm: ../../data/mage/v47_kallisto_quant/{sample}/matrix.abundance.tpm.mtx
      matrix_tsv: ../../data/mage/v47_kallisto_quant/{sample}/matrix.abundance.tsv
      matrix_tpm_tsv: ../../data/mage/v47_kallisto_quant/{sample}/matrix.abundance.tpm.tsv
      merge_matrix_tsv: ../../data/mage/v47_kallisto_quant/matrix.abundance.tsv
      merge_matrix_tpm_tsv: ../../data/mage/v47_kallisto_quant/matrix.abundance.tpm.tsv

  enh_v47_kallisto:
    odir: ../../data/mage/enh_v47_kallisto/{sample}/
    bus: ../../data/mage/enh_v47_kallisto/{sample}/output.bus
    transcripts: ../../data/mage/enh_v47_kallisto/{sample}/transcripts.txt
    matrix: ../../data/mage/enh_v47_kallisto/{sample}/matrix.ec
    flens: ../../data/mage/enh_v47_kallisto/{sample}/flens.txt
    bus_sort: ../../data/mage/enh_v47_kallisto/{sample}/output_sorted.bus

    count_pref: ../../data/mage/enh_v47_kallisto/{sample}/count
    count_mtx: ../../data/mage/enh_v47_kallisto/{sample}/count.mtx
    count_ec: ../../data/mage/enh_v47_kallisto/{sample}/count.ec.txt
    matrix_tsv: ../../data/mage/enh_v47_kallisto/{sample}/matrix.abundance.tsv

    quant:
      odir: ../../data/mage/enh_v47_kallisto_quant/{sample}
      matrix: ../../data/mage/enh_v47_kallisto_quant/{sample}/matrix.abundance.mtx
      matrix_tpm: ../../data/mage/enh_v47_kallisto_quant/{sample}/matrix.abundance.tpm.mtx
      matrix_tsv: ../../data/mage/enh_v47_kallisto_quant/{sample}/matrix.abundance.tsv
      matrix_tpm_tsv: ../../data/mage/enh_v47_kallisto_quant/{sample}/matrix.abundance.tpm.tsv
      merge_matrix_tsv: ../../data/mage/enh_v47_kallisto_quant/matrix.abundance.tsv
      merge_matrix_tpm_tsv: ../../data/mage/enh_v47_kallisto_quant/matrix.abundance.tpm.tsv

module load bedtools
bedtools intersect             -s             -wao             -a ../../data/personal_genome/gtf/HG002.fq_HG002.fa_transcript_models.gtf  -b ../../data/personal_genome/gtf/HG002.fq_HG002.fa_transcript_models.gtf
```

Trying to make the weird input files:
```bash
Rscript /gpfs/projects/bsc83/Projects/pantranscriptome/fairlie/240903_pt/scripts//prepare_lineage_genotype.R \
  --genotype ../../data/mage/poder_kallisto/sqtl/1kg_mage_cis_norm_snp_biallelic.vcf.gz  \
  --metadata ../../data/mage/poder_kallisto/sqtl/metadata.tsv \
  --in-012 ../../data/mage/poder_kallisto/sqtl/1kg_mage_cis.012 \
  --indv-012 ../../data/mage/poder_kallisto/sqtl/1kg_mage_cis.012.indv \
  --pos-012 ../../data/mage/poder_kallisto/sqtl/1kg_mage_cis.012.pos \
  --out-file ../../data/mage/poder_kallisto/sqtl/genotype.tsv.gz \
  --verbose
```

```bash
d=/gpfs/projects/bsc83/Projects/pantranscriptome/fairlie/240903_pt/snakemake/mage/work/8c/e6d8a42a96a4f5a7161f82dbc7a17b
d2=/gpfs/projects/bsc83/Projects/scRNAseq/imestres/sQTLseeker/sqtlseeker2-nf/bin/
Rscript $d2/prepare_trexp.R \
  --transcript_expr $d/matrix.abundance.tpm.tsv \
  --metadata $d/metadata.tsv \
  --gene_location $d/gene_bodies.tsv  \
  --group "1" \
  --covariates \
  --min_gene_expr 1  \
  --min_proportion 0.4 \
  --min_transcript_expr 0.1 \
  --min_dispersion 0.1  \
  --output_tre tre.df.RData \
  --output_gene genes.ss.bed
```
```R
library(dplyr)
library("data.table")
genes.bed.f="/gpfs/projects/bsc83/Projects/pantranscriptome/fairlie/240903_pt/snakemake/mage/work/8c/e6d8a42a96a4f5a7161f82dbc7a17b/gene_bodies.tsv"
metadata.f='/gpfs/projects/bsc83/Projects/pantranscriptome/fairlie/240903_pt/snakemake/mage/work/8c/e6d8a42a96a4f5a7161f82dbc7a17b/metadata.tsv'
trans.expr.f="/gpfs/projects/bsc83/Projects/pantranscriptome/fairlie/240903_pt/snakemake/mage/work/8c/e6d8a42a96a4f5a7161f82dbc7a17b/matrix.abundance.tpm.tsv"
sel.group = "1"

# this stuff is just right from prepare_trexp.R
genes.bed <- read.table(genes.bed.f, header = TRUE, as.is = TRUE, sep = "\t")

## 3. Getting the IDs of the samples in the group of interest

metadata <- read.table(metadata.f, header = TRUE,
                       as.is = TRUE, sep = "\t")

subset.df <- subset(metadata, group == sel.group)                                   # Select samples of interest
subset.samples <- subset.df$sampleId

## 4. Prepare transcript expression

if( grepl("\\.gz$", trans.expr.f) ){
    trans.expr.f <- paste0("zcat < '", trans.expr.f, "'")
}

te.df <- as.data.frame(fread(input=trans.expr.f, header = TRUE, sep = "\t"))
colnames(te.df)[1:2] <- c("trId", "geneId")                                         # Proper 1,2 colnames
subset.samples <- subset.samples[subset.samples%in%colnames(te.df)]                 # Get samples that have quantifications

te.df <- te.df[, c("trId", "geneId", subset.samples)]                               # Select subgroup of samples = the ones from the group of interest

# try before starting to subset everything

# hellow i'm the problem
te.df <- subset(te.df, geneId %in% genes.bed$geneId)         
```

```bash
module load nextflow/21.04.0
module load singularity/3.11.5
module load R/4.3.2

nextflow run /gpfs/projects/bsc83/Projects/scRNAseq/imestres/sQTLseeker/sqtlseeker2-nf/sqtlseeker2.nf \
--genotype ../../data/mage/v47_kallisto/sqtl/genotype.tsv.gz \
--trexp ../../data/mage/v47_kallisto/sqtl/matrix.abundance.tpm.tsv \
--metadata ../../data/mage/v47_kallisto/sqtl/metadata.tsv \
--genes ../../data/mage/v47_kallisto/sqtl/gene_bodies.tsv \
--dir ../../data/mage/v47_kallisto/sqtlseeker/ \
--mode "nominal" \
--covariates "true" \
  --svqtl "true" \
  --fdr 0.01 \
  --fdr_svqtl 0.01 \
  --min_gene_exp 1 \
  --min_transcript_exp 0.1 \
  --min_md 0.05 \
  --min_proportion 0.8 \
  --with-singularity /gpfs/projects/bsc83/Projects/scRNAseq/imestres/sQTLseeker/runs/sqtlseeker2-nf.sif
```

```bash
bedtools intersect -loj -wa  -a ../../data/hprc/kinnex/map/merge/HG03704_q10_tss.bed -b ../../ref/ccre/pls.bed.gz > ../../data/hprc/kinnex/map/merge/HG03704_tss_pls_int.tsv
```


240611 testing flip reads + call tss bed locations code
```bash
bam=/gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/03_mapping/data/assembly_mapping/genomic/25_IS2_GM22299.bam
module load samtools
samtools view -c $bam # 13,693,804
samtools view $bam | grep "ts:" |  wc -l # like 12 million but forgot to copy
# so these have ts tags

bam=~/mele_lab/projects/240903_pt/data/hprc/kinnex/map/merge/HG03704_q10.bam
samtools view -c $bam # 40,195,562
samtools view $bam | grep "ts:" |  wc -l # 31,349,064




```
