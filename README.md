# CRCbiome_virome_2023

This repository contains code for processing of sequencing data and virus predictions, and data analysis used to generate results in Istvan et al. For details on this analysis see the preprint at [medRxiv](https://www.medrxiv.org/content/10.1101/2023.08.24.23294548v1).

## Bioinformatic pipeline
The input files for the pipeline are
- VirSorter output files (from sample-based assembly)
- Contig statistics
- Paired end reads passing QC
- A file with a list of sample names corresponding to the preceding files.

The bioinformatic analyses are written for the Linux OS as a Snakemake pipeline, and are intended to be run on an HPC cluster system. Software dependencies provided in conda environment specification files in `workflow/envs/`. 

A more user-friendly version of this pipeline derived from the current pipeline is available as [VirMake](https://github.com/Rounge-lab/VirMake).

## Scripts for post-processing analysis
R scripts for analysis of data are available in `analyses/`. Required R packages include `tidyverse`, `vegan`, `ape`, `gridExtra`, `broom`, `effectsize`, `ggrepel`, `UpSetR`, and `MaAsLin2`. 

To test the scripts, mock datasets for 1000 samples and a subset of the viral genomes identified in the study are available in the folder `mock_datasets`. To run the a mock analysis, run the script called `analyses/scripts/run_mock_analysis.R`. The working directory for the script should be the repository folder. Running the mock datasets should take less than one hour on a regular desktop computer and will result in generation of statsitical tests, figures and tables reported in the manuscript, albeit with limited associations due to the random generation of input data.


