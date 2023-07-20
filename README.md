# COVID-19 Severity-Related Gene Expression Analysis

This repository contains the R scripts and data files used in the study "Leveraging RNA-seq Data to Understand the Molecular Mechanisms Associated with Varying Levels of COVID-19 Severity". The research identifies dynamic gene expression patterns correlated with disease progression using both bulk RNA-seq and single-cell RNA-seq (scRNA-seq) data from COVID-19 patients and healthy individuals.

## Repository Structure

The repository is organized as follows:

- `Identification and visualization of gene expression patterns.R`: This R script is used to identify and visualize the gene expression patterns associated with COVID-19 severity based on RNA-seq data. Four distinct expression trends have been identified.

- `Gene_enrichment.R`: This R script is used to conduct pathway enrichment analyses on the genes of interest. The analysis revealed a range of impacted pathways, from immune response to circadian rhythm.

- `scRNA-seq analysis.R`: This R script is used to perform in-depth single-cell analysis. The results provide further insights into the heterogeneity and specific cellular localization of the different expression patterns.

- `model_gene_plot.pdf`: This is the output visualization of the gene expression model, showing the distinct gene expression trends associated with COVID-19 severity.

- `Required files`: This folder contains necessary data files for the scripts:
  - `all_results.Rdata`: This file stores the results from the analyses.
  - `GSE152418_filtered_metadata.csv`: This file contains the filtered metadata for the GSE152418 dataset used in the analyses.
  - `GSE152418_GeneLevel_Raw_data.csv`: This file contains the raw gene-level data for the GSE152418 dataset.

## Usage

To reproduce the analyses, first, clone this repository. Make sure you have the latest version of R and necessary R packages installed. Then, simply run the R scripts in the order provided above. Note that you might need to set your working directory to the location where the scripts are saved. 

## Related Resources

All findings have been consolidated in the COVID-19 Severity Related Database (http://covid.haoran.pub/), a web-based platform that aims to facilitate further research and potential therapeutic developments.

## Contact

For any questions or concerns, please open an issue in this repository, or contact the maintainers directly.

## Acknowledgments

We would like to thank all those who have contributed to this project, and the brave individuals who participated in the studies that provided the data for our analyses.
