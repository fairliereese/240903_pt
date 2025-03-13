## This document outlines the locations for analyses done pertaining to the section "A population-diverse long-read RNA-seq dataset enables the discovery of many high-confidence novel transcripts" (ie generation and description of PODER).

#### Notes:
* Most non-UpSet plot figures made by me were re-plotted in R for better aesthetics, and thus the visualizations in the following notebooks, though correct, look different. Similarly, sometimes the output is just a saved table used later to plot in R.
* Similarly, many statistical tests in these notebooks were re-run in R for compatibility with plotting.
* Any notebooks not explicitly described here do not contain work used to write the manuscript.


* [`compare_annotation_to_external.ipynb`](https://github.com/fairliereese/240903_pt/blob/main/analysis/figure_characterize_poder/compare_annotation_to_external.ipynb): UpSet plots for overlap with external catalogs, including by novelty category.
Plots showing # and % of novel transcripts from PODER, CHESS, GTEx, and ENCODE that are unsupported between one another.
* [`external_sqanti_cats.ipynb`](https://github.com/fairliereese/240903_pt/blob/main/analysis/external_sqanti_cats.ipynb): SQANTI categories for the different catalogs (ENCODE, CHESS, GTEx).
* [`characterize_novel_exons.ipynb`](https://github.com/fairliereese/240903_pt/blob/main/analysis/characterize_novel_exons.ipynb): Finding, characterization, and counting of novel exons in PODER.
* [`explore_orfs_aa_3.ipynb`](https://github.com/fairliereese/240903_pt/blob/main/analysis/explore_orfs_aa_3.ipynb): Predicted protein-coding sequence results.
* [`ggtranscript_hlab.ipynb`](https://github.com/fairliereese/240903_pt/blob/main/analysis/ggtranscript_hlab.ipynb): Transcript browser plot for `HLA-B`.
* [`tau.ipynb`](https://github.com/fairliereese/240903_pt/blob/main/analysis/tau.ipynb): Compute Tau values using PODER lr-kallisto quantification values.
* [`tau_vs_pop_specific.ipynb`](https://github.com/fairliereese/240903_pt/blob/main/analysis/tau_vs_pop_specific.ipynb): Statistical tests for Tau population-specific transcripts versus population-specific discovered transcripts.
* [`td_personal_ic_n_novel_not_hg38_over_novel_hg38.ipynb`](https://github.com/fairliereese/240903_pt/blob/main/analysis/https://github.com/fairliereese/240903_pt/blob/main/analysis/td_personal_ic_n_novel_not_hg38_over_novel_hg38.ipynb): Number of novel ICs we can find using personalized-GRCh38s vs. GRCh38 for 30 1000G samples. UpSet intersections of detection of novel transcripts between personalized-GRCh38s and GRCh38 alone.
