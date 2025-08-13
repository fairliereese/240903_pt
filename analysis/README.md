## This document outlines the locations for analyses done for the manuscript

#### Notes:
* Most non-UpSet plot figures made by me were re-plotted in R for better aesthetics, and thus the visualizations in the following notebooks, though correct, look different. Similarly, sometimes the output is just a saved table used later to plot in R.
* Similarly, many statistical tests in these notebooks were re-run in R for compatibility with plotting.
* Any notebooks not explicitly described here do not contain work used to write the manuscript.

#### *A population-diverse long-read RNA-seq dataset enables the discovery of many high-confidence novel transcripts*; (ie generation and description of PODER).

* [`spliced_sirv_sqanti_reads_filter.ipynb`](https://github.com/fairliereese/240903_pt/blob/main/analysis/spliced_sirv_sqanti_reads_filter.ipynb): Analysis of spliced SIRVs to validate our experimental and bioinformatic pipelines.
* [`long_sirv_sqanti_reads.ipynb`](https://github.com/fairliereese/240903_pt/blob/main/analysis/long_sirv_sqanti_reads.ipynb): Analysis of long SIRVs to explore potential read length / coverage biases.s
* [`compare_annotation_to_external.ipynb`](https://github.com/fairliereese/240903_pt/blob/main/analysis/figure_characterize_poder/compare_annotation_to_external.ipynb): UpSet plots for overlap with external catalogs, including by novelty category.
Plots showing # and % of novel transcripts from PODER, CHESS, GTEx, and ENCODE that are unsupported between one another.
* [`external_sqanti_cats.ipynb`](https://github.com/fairliereese/240903_pt/blob/main/analysis/external_sqanti_cats.ipynb): SQANTI categories for the different catalogs (ENCODE, CHESS, GTEx).
* [`characterize_novel_exons.ipynb`](https://github.com/fairliereese/240903_pt/blob/main/analysis/characterize_novel_exons.ipynb): Finding, characterization, and counting of novel exons in PODER.
* [`pop_div_compare_novel_exon_parts_to_known.ipynb`](https://github.com/fairliereese/240903_pt/blob/main/analysis/pop_div_compare_novel_exon_parts_to_known.ipynb): Comparing CEU <-> all other 1000G populations for biallelic SNPs present in known and novel exon parts.
* [`explore_orfs_aa_3.ipynb`](https://github.com/fairliereese/240903_pt/blob/main/analysis/explore_orfs_aa_3.ipynb): Predicted protein-coding sequence results.
* [`ggtranscript_hlab.ipynb`](https://github.com/fairliereese/240903_pt/blob/main/analysis/ggtranscript_hlab.ipynb): Transcript browser plot for `HLA-B`.
* [`compare_total_mage_counts.ipynb`](https://github.com/fairliereese/240903_pt/blob/main/analysis/compare_total_mage_counts.ipynb): Computing number of reads used to quantify MAGE data using Enhanced GENCODE versus GENCODE.
* [`gene_len.ipynb`](https://github.com/fairliereese/240903_pt/blob/main/analysis/gene_len.ipynb): Explore relationship between gene length and different transcript or gene features in PODER.
* [`transcript_len.ipynb`](https://github.com/fairliereese/240903_pt/blob/main/analysis/transcript_len.ipynb): Explore relationship between transcript length and different transcript or gene features in PODER.
* [`tool_sharing_upset.ipynb`](https://github.com/fairliereese/240903_pt/blob/main/analysis/tool_sharing_upset.ipynb): Plot the intersection for tool sharing for transcripts in PODER (ie, post-filtering).




#### *Current gene annotations are less representative of transcriptomes from non-European descent individuals*; (description and characterization of non-European bias).

* [`tau.ipynb`](https://github.com/fairliereese/240903_pt/blob/main/analysis/tau.ipynb): Compute Tau values using PODER lr-kallisto quantification values.
* [`tau_vs_pop_specific.ipynb`](https://github.com/fairliereese/240903_pt/blob/main/analysis/tau_vs_pop_specific.ipynb): Statistical tests for Tau population-specific transcripts versus population-specific discovered transcripts.
* [`tau_mage_enh_gencode.ipynb`](https://github.com/fairliereese/240903_pt/blob/main/analysis/tau_mage_enh_gencode.ipynb): Compute Tau values of population specificity using the kallisto + Enhanced GENCODE quantification of the MAGE RNA-seq dataset. Compare population specificity across our dataset and MAGE.
* [`pop_shared_t.ipynb`](https://github.com/fairliereese/240903_pt/blob/main/analysis/pop_shared_t.ipynb): Examine attributes of sample / population-shared transcripts and their distributions of detection across samples / populations.
- [`perc_nov_pop_spec.ipynb`](https://github.com/fairliereese/240903_pt/blob/main/analysis/perc_nov_pop_spec.ipynb): Compute % of novel transcripts per category that are population-specific.

#### *A population diverse-annotation enhances discovery of allele-specific transcript usage in non-European populations*; description and characterization of the differential findings for allele-specific analyses using different reference gene annotations.
* [`n_transcripts_per_annot.ipynb`](https://github.com/fairliereese/240903_pt/blob/main/analysis/n_transcripts_per_annot): Compute number of transcripts per reference gene annotation.

#### *Personalized genome reference assemblies enhance novel splice junction discovery*.
* [`td_personal_ic_n_novel_not_hg38_over_novel_hg38.ipynb`](https://github.com/fairliereese/240903_pt/blob/main/analysis/td_personal_ic_n_novel_not_hg38_over_novel_hg38.ipynb): Number of novel ICs we can find using personalized-GRCh38s vs. GRCh38 for 30 1000G samples. UpSet intersections of detection of novel transcripts between personalized-GRCh38s and GRCh38 alone.
* [`td_personal_relative_gain_p9.ipynb`](https://github.com/fairliereese/240903_pt/blob/main/analysis/td_personal_relative_gain_p9.ipynb): Look at saturation of transcript discovery with increasing numbers of samples, using either the GRCh38 mapping / transcript discovery strategy or the personalized-GRCh38 strategy.
* [`td_personal_gain_28_sample_experiment_p9.ipynb.ipynb`](https://github.com/fairliereese/240903_pt/blob/main/analysis/td_personal_gain_28_sample_experiment_p9.ipynb.ipynb): Examine the gain in unique transcripts starting from 28 merge_espresso_downsample_refseq discovered based on performing GRCh38 mapping / transcript discovery strategy on 2 samples (2x lrRNA-seq) or the personalized-GRCh38 strategy on 1 sample (1x lrRNA-seq + whole genome sequencing).
