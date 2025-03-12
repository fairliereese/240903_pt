## Data processing

All computationally-intense data processing was done using Snakemake. Each somewhat-distinct task is in its own folder to make it possible to run several Snakemakes in parallel with one another.

## Organization

In the parent directory, there are a few shared resources
* [`config.yml`](https://github.com/fairliereese/240903_pt/blob/main/snakemake/config.yml): Defines all the locations of the files used throughout this repo (including, notable, in the [analysis](https://github.com/fairliereese/240903_pt/tree/main/analysis) folder as well).
* [common](https://github.com/fairliereese/240903_pt/tree/main/snakemake/common): Snakemake rule definitions for rules that are repeatedly used throughout the subfolder tasks.

In each subdirectory, the important files which are usually there are as follows:
* `Snakefile`: Used to run the data processing / analysis workflow
* `snakefile_dev.ipynb` (and other `*dev.ipynb`): Jupyter notebooks used to debug input / output Snakemake files or other tasks that are run during `Snakefile` execution.
* Oftentimes there are additional `*.txt`, `*.md`, or `*.tsv` files that help outline additional input / output information (especially related to external dataset use) or other information / code used to run the `Snakefile`.

## Subfolder descriptions

Subfolders that are not listed here to not contains analyses / processing that was ultimately used to write this manuscript.

<!-- * [1000g](https://github.com/fairliereese/240903_pt/tree/main/snakemake/1000g):  -->
<!-- * [lapa](https://github.com/fairliereese/240903_pt/tree/main/snakemake/lapa): Quantification / identification of TSS usage from LR-RNA-seq -->
<!-- Also attempt at running sQTLseeker and suppa -->
<!-- * [map](https://github.com/fairliereese/240903_pt/tree/main/snakemake/map): Run mapping and compute mapping statistics using T2T, GRCh38, and [African CAAPA contigs](https://www.biorxiv.org/content/10.1101/2023.11.04.564839v1). -->
<!-- * [personal_genome](https://github.com/fairliereese/240903_pt/tree/main/snakemake/personal_genome): What is this -->
<!-- * [personal_lr-kallisto](https://github.com/fairliereese/240903_pt/tree/main/snakemake/personal_lr-kallisto): Run lr-kallisto on personal assembly-mapped 6 samples using PODER liftOff to personal assembly.. -->
<!-- * [suppa](https://github.com/fairliereese/240903_pt/blob/main/snakemake/suppa): Run SUPPA on our LR-RNA-seq data to quantify alternative splicing events. -->
<!-- * [unmerged_lr-kallisto](https://github.com/fairliereese/240903_pt/tree/main/snakemake/unmerged_lr-kallisto): Running lr-kallisto on FASTQs before they were merged across sequencing runs for the same sample... for some reason. -->
<!-- * [unmerged_v47_lr-kallisto](https://github.com/fairliereese/240903_pt/tree/main/snakemake/unmerged_v47_lr-kallisto): Running lr-kallisto on FASTQs before they were merged across sequencing runs for the same sample... for some reason. Using GENCODE v47 as reference annotation. -->
<!-- * [v47_lr-kallisto](https://github.com/fairliereese/240903_pt/tree/main/snakemake/v47_lr-kallisto): Quantification of our LR-RNA-seq dataset with lr-kallisto and GENCODE v47.  -->
<!-- * [v47_personal_lr-kallisto](https://github.com/fairliereese/240903_pt/blob/main/snakemake/v47_personal_lr-kallisto): Run lr-kallisto on personal assembly-mapped 6 samples using GENCODE v47 liftOff to personal assembly. -->
* [ics_inter_catalog](https://github.com/fairliereese/240903_pt/tree/main/snakemake/ics_inter_catalog): Extract unique intron chains from each catalog (CHESS3, ENCODE4, GTEx) and run SQANTI on them.
* [astu_example](https://github.com/fairliereese/240903_pt/tree/main/snakemake/astu_example): Systematically output LR-RNA-seq BAM files split by allele to facilitate visualization of allele-specific transcript usage (ASTU).
* [encode](https://github.com/fairliereese/240903_pt/tree/main/snakemake/encode): Dertermine genetic ancestry (if annotated) of ENCODE LR-RNA-seq catalog. (Not Snakemake).
* [gtex_lr-kallisto](https://github.com/fairliereese/240903_pt/tree/main/snakemake/gtex_lr-kallisto): Re-quantify the GTEx LR-RNA-seq dataset from FASTQ using lr-kallisto and PODER.
* [lr-kallisto](https://github.com/fairliereese/240903_pt/tree/main/snakemake/lr-kallisto): Quantification of our LR-RNA-seq dataset with lr-kallisto and PODER.
* [mage](https://github.com/fairliereese/240903_pt/tree/main/snakemake/mage): Quantification of [MAGE](https://github.com/mccoy-lab/MAGE) RNA-seq dataset using kallisto and multiple annotations (PODER, GENCODE, Enhanced GENCODE).
* [map_personal](https://github.com/fairliereese/240903_pt/tree/main/snakemake/map_personal): Download of personal assemblies from pangenome matching our LR-RNA-seq samples. <!-- Also mapping but we ended up using Fabien's -->
* [merge_espresso](https://github.com/fairliereese/240903_pt/tree/main/snakemake/merge_espresso): Merging ESPRESSO transcripts from downsampling experiment when using GENCODE v47 as the annotation by intron chain.
* [merge_espresso_refseq](https://github.com/fairliereese/240903_pt/tree/main/snakemake/merge_espresso_refseq): Merging ESPRESSO transcripts from downsampling experiment when using RefSeq v110 as the annotation by intron chain.
* [merge_v47_poder](https://github.com/fairliereese/240903_pt/blob/main/snakemake/merge_v47_poder): Merge novel PODER transcript GTF with GENCODE v47 GTF to create the Enhanced GENCODE annotation.
* [novel_annotation_add_gene](https://github.com/fairliereese/240903_pt/blob/main/snakemake/novel_annotation_add_gene): Add gene entries to PODER GTF which only had transcript and exon entries, and a gene ID tag.
* [novel_gene](https://github.com/fairliereese/240903_pt/tree/main/snakemake/novel_gene): Add novel gene IDs to intergenic transcripts found in PODER using [buildLoci](https://github.com/julienlag/buildLoci).
* [pfam](https://github.com/fairliereese/240903_pt/blob/main/snakemake/pfam): Run PFAM protein domain finder on GENCODE v47 reference annotation and PODER predicted proteins.
* [poder_protein](https://github.com/fairliereese/240903_pt/tree/main/snakemake/poder_protein): Run protein prediction pipeline on PODER.
* [pop_div_exon_fsts](https://github.com/fairliereese/240903_pt/tree/main/snakemake/pop_div_exon_fsts): Find known exons and novel transcribed regions of novel exons. Compute the Fst values pairwise between each set of populations in the 1000G that overlap ours for the SNP variants that fall within these regions.
* [protein](https://github.com/fairliereese/240903_pt/tree/main/snakemake/protein): Run protein prediction pipeline on UMA.
* [transcript_discovery_personal](https://github.com/fairliereese/240903_pt/tree/main/snakemake/transcript_discovery_personal): Map our LR-RNA-seq data from the 30 1000G-overlapping samples to their corresponding 2 personalized-GRCh38 haplotypes (ie GRCh38 with the SNPs from each sample incorporated). Run SQANTI on the ESPRESSO results to get splice junctions from each mapping (personalized-GRCh38s and GRCh38 alone). Intersect biallelic SNPs with splice-junction proximal exonic SNPs and splice site intronic SNPs to determine how splice junctions' discovery is affected by genetics.



## Snakemake calls

```bash
# normal
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
    --time=48:00:00" \
    -n


# gpu version
snakemake \
  -s Snakefile \
  -j 100 \
  --latency-wait 120 \
  --use-conda \
  --cluster \
    "sbatch \
    --nodes {resources.nodes} \
    -q acc_bscls \
    -A bsc83 \
    -c {resources.threads}  \
    --mail-user=freese@bsc.es \
    --mail-type=START,END,FAIL \
    --time=2:00:00" \
    -n

# debug version
conda activate pt_snakemake
snakemake \
  -s Snakefile \
  -j 100 \
  --latency-wait 120 \
  --use-conda \
  --cluster \
    "sbatch \
    --nodes {resources.nodes} \
    -q gp_debug \
    -A bsc83 \
    -c {resources.threads}  \
    --mail-user=freese@bsc.es \
    --mail-type=START,END,FAIL \
    --time=2:00:00" \
    -n

# non-cluster version
snakemake \
  -s Snakefile \
  -j 100 \
  --latency-wait 120 -n

# plotting version
snakemake --forceall --dag | dot -Tpdf > dag.pdf
```
