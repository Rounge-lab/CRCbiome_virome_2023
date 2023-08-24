# CRCbiome_virome_2023

This repository contains code for processing of sequencing data and virus predictions, and data analysis used to generate results in Istvan et al. For details on this analysis see (link to biorxiv).

The input files for the pipeline are
- VirSorter output files (from sample-based assembly)
- Contig statistics
- Paired end reads passing QC
- A file with a list of sample names corresponding to the preceding files.

The pipeline is run using snakemake.

A more user-friendly version of this pipeline derived from the current pipeline is available as [VirMake](https://github.com/Rounge-lab/VirMake).

Scripts for analysis of data are available in `analyses/`. These depend on information on samples and individuals that is unfortunately not possible to share. 

