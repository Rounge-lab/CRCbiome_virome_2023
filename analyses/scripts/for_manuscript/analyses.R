


# env ---------------------------------------------------------------------

library(tidyverse)
library(vegan)
library(ape)
library(gridExtra)
library(broom)
library(effectsize)
library(ggrepel)
library(UpSetR)

script_folder <- "analyses/scripts/"

if ("mock" %in% ls()) script_folder <- paste0("../", script_folder)

# functions ---------------------------------------------------------------

source(paste0(script_folder, "utils/utils.R"))

# preprocessing -----------------------------------------------------------

# source("analyses/scripts/preprocess_datasets.R")

# load data ---------------------------------------------------------------

sample_data <- read_rds("data/participant_data/sample_meta.Rds")

meta_dat <- read_rds("data/participant_data/metadata_selected_variables.Rds")
meta_variables <- read_rds("data/participant_data/metadata_variables.Rds")

virus_abundance <- read_rds("data/abundance_tables/viral_abundance.Rds")

# fraction_viral <- read_tsv("data/abundance_tables/fraction_viral_contigs.tsv")

vOTU_stats <- read_tsv("data/vOTU_stats/vOTU_stats.tsv")
# vOTU_lifecycle_cat <- read_tsv("data/vOTU_stats/vOTU_lyfecycle_cat.tsv")
amgs_cat_per_vOTU <- read_tsv("data/amgs/amg_cat_by_vOTU_mid_min75.tsv")

checkv_res_all <- read_rds("data/abundance_tables/checkV_results_all.Rds")
checkV_res_all_by_sample <- read_rds("data/abundance_tables/checkV_results_by_sample.Rds")


# Data generation tasks ---------------------------------------------------

## Summarize vOTU stats
source(paste0(script_folder, "for_manuscript/vOTU_based_stats.R"))
create_vOTU_based_stats()

## Perform diversity and composition analyses
## Do PERMANOVA analyses
source(paste0(script_folder, "auxiliary/alpha_beta_div_for_lifestyle_demography.R"))
alpha_diversity_calculations()
beta_diversity_calculations()

# presentations -----------------------------------------------------------

## Create contig quality plot
source(paste0(script_folder, "for_manuscript/contig_qual_plot.R"))
plot_contig_qual()

## Create technical plot
source(paste0(script_folder, "for_manuscript/technical_plot.R"))
create_technical_fig()

## Create annotation plot
source(paste0(script_folder, "for_manuscript/vOTU_descriptive_plot.R"))
make_vOTU_descriptive_plots()
create_amg_plots()

## Create host associations plot
source(paste0(script_folder, "for_manuscript/host_associations_plot.R"))
create_host_associations_figures()

