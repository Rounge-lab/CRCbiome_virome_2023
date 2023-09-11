
library(tidyverse)

norm_with_min <- function(n, m, s, min_val) {
  rnorm(n = n, mean = m, sd = s) %>% enframe() %>% mutate(value = case_when(value < min_val ~ value+abs(min(value))+min_val, TRUE ~ value)) %>% pull(value)
}

create_mock_dataset <- function(dataset, n, seed = 1, subject_ids = NA, sample_ids = NA) {
  set.seed(seed)
  if (any(is.na(subject_ids))) subject_ids <- sample(seq(1e5, 6e5), size = n, replace = FALSE)
  if (any(is.na(sample_ids))) sample_ids <- paste("S", sample(seq(5.8e6, 1.9e7), size = n, replace = FALSE), sep = "_")
  
  if (dataset %in% "sample_info") {
    return(tibble("sample_id" = sample_ids,
                  "rekv_nr" = str_remove(sample_id, "S_") %>% as.integer(),
                  "Plate_Nr" = paste0("Plate_", sample(seq(1, 23), size = n, replace = TRUE)),
                  "Well" = paste(sample(LETTERS[1:8], size = n, replace = TRUE), sample(seq(12), size = n, replace = TRUE), sep = "_"),
                  "run" = sample(c(1,2,3,4,"multiple"), size = n, prob = c(rep(0.245,4),0.02), replace = TRUE),
                  "deltaker_id" = subject_ids,
                  "invitasjon_id" = sample(seq(7e5, 9.9e5), size = n, replace = FALSE),
                  "ffq_ref_nr" = sample(seq(1e4, 2e5), size = n, replace = FALSE),
                  "Prøvetype" = "Baseline",
                  "investigation_date" = sample(as.Date("2017-10-01"):as.Date("2021-04-01"), size = n, replace = TRUE) %>% as.Date(origin = "1970-01-01"),
                  "collected_date" = investigation_date - sample(0:10, size = n, replace = TRUE),
                  "FIT_value" = rexp(n*1000, rate = 0.1) %>% enframe() %>% mutate(value = value/0.01) %>% filter(value > 76) %>% slice_sample(n = n) %>% pull(value),
                  "ekstraksjonsbatch" = sample(1:100, size = n, replace = TRUE),
                  "ekstraksjons_dato" = investigation_date + sample(30:300, size = n, replace = TRUE),
                  "qubit_total_dna" = norm_with_min(n = n, m = 5.5, s = 3.5, min_val = 0.7),
                  "nanodrop_total_dna" = runif(n = n, min = 0.4, max = 1.1)*qubit_total_dna,
                  "sampling_to_extraction" = as.integer(investigation_date - collected_date),
                  "Total_Bases_QC_ATLAS" = norm_with_min(n = n, m = 3e9, s = 1e9, min_val = 1e9),
                  "Total_Reads_QC_ATLAS" = Total_Bases_QC_ATLAS/(2*151),
                  "N50" = norm_with_min(n = n, m = 3000, s = 1500, min_val = 300),
                  "n_contigs" = norm_with_min(n = n, m = 50000, s = 15000, min_val = 15000)))
    
  }
  
  if (dataset %in% "metadata") {
    return(tibble("deltaker_id" = subject_ids,
                  "kjonn" = sample(c("Female", "Male"), size = n, replace = TRUE) %>% factor(levels = c("Female", "Male")),
                  "age_cat" = sample(c("50-60", "60-70", "70-80"), size = n, replace = TRUE) %>% factor(levels = c("50-60", "60-70", "70-80")),
                  "senter" = sample(c("Bærum", "Moss"), size = n, replace = TRUE) %>% factor(levels = c("Bærum", "Moss")),
                  "BMI" = rnorm(n = n, mean = 27, sd = 4.1),
                  "Energi_kcal" = rnorm(n = n, mean = 2222, sd = 650),
                  "Prot_energi" = rnorm(n = n, mean = 16.5, sd = 2.5),
                  "Karboh_energi" = rnorm(n = n, mean = 41.5, sd = 6.67),
                  "Sukker_energi" = norm_with_min(n = n, m = 5.1, s = 3.9, min_val = 0),
                  "Fiber" = norm_with_min(n = n, m = 29.1, s = 10.3, min_val = 0),
                  "Fett_energi" = norm_with_min(n = n, m = 34.7, s = 6, min_val = 0),
                  "Mettet_energi" = norm_with_min(n = n, m = 12, s = 2.7, min_val = 0),
                  "C_enum_energi" = norm_with_min(n = n, m = 13.1, s = 2.9, min_val = 0),
                  "C_flerum_energi" = norm_with_min(n = n, m = 6.4, s = 1.6, min_val = 0),
                  "Trans_u_energi" = norm_with_min(n = n, m = 0.3, s = 0.13, min_val = 0),
                  "Alko" = norm_with_min(n = n, m = 13.4, s = 16.5, min_val = 0),
                  "wcrf_index_main" = norm_with_min(n = n, m = 3.5, s = 1, min_val = 0),
                  "PhysAct_Score" = norm_with_min(n = n, m = 199, s = 223, min_val = 0),
                  "Smoking" = sample(c("Non-smoker", "Smoker", "Missing"), 
                                     size = n, replace = TRUE, prob = c(rep(0.95/2, 2), 0.05)) %>% 
                    factor(levels = c("Non-smoker", "Smoker", "Missing")),
                  "Snus" = sample(c("Non-snuser", "Snuser", "Missing"), 
                                  size = n, replace = TRUE, prob = c(rep(0.95/2, 2), 0.05)) %>% 
                    factor(levels = c("Non-snuser", "Snuser", "Missing")),
                  "Utdanning" = sample(c("Primary school", "University/college", "High school", "Missing"), 
                                       size = n, replace = TRUE, prob = c(rep(0.95/3, 3), 0.05)) %>% 
                    factor(levels = c("Primary school", "High school", "University/college", "Missing")),
                  "Sivilstatus_cat2" = sample(c("Not married/non-cohabiting", "Married/cohabiting", "Missing"), 
                                              size = n, replace = TRUE, prob = c(0.2, 0.75, 0.05)) %>% 
                    factor(levels = c("Not married/non-cohabiting", "Married/cohabiting", "Missing")),
                  "Arbeid_lump" = sample(c("Employed", "Retired", "Other"), 
                                         size = n, replace = TRUE, prob = c(0.3, 0.55, 0.15)) %>% 
                    factor(levels = c("Employed", "Retired", "Other")),
                  "Nasj_cat2" = sample(c("Native", "Non-native", "Missing"), 
                                       size = n, replace = TRUE, prob = c(0.9, 0.075, 0.025)) %>% 
                    factor(levels = c("Native", "Non-native", "Missing")),
                  "Antibiotics" = sample(c("No", "Yes", "Missing"), 
                                         size = n, replace = TRUE, prob = c(0.8, 0.15, 0.5)) %>% 
                    factor(levels = c("No", "Yes", "Missing")),
                  "Antacids" = sample(c("No", "Yes", "Missing"), 
                                      size = n, replace = TRUE, prob = c(0.7, 0.25, 0.5)) %>% 
                    factor(levels = c("No", "Yes", "Missing"))
    ))
  }
  if (dataset %in% "screening") {
    return(tibble("deltaker_id" = subject_ids,
                  "final_result" = sample(c("not cancer", "6. Cancer"), size = n, replace = TRUE, prob = c(0.95, 0.05))))
  }
  
  if (dataset %in% "abundance") {
    
    
    return(lapply(seq(ceiling(n/nrow(virus_abundance))), function(x) {
      virus_abundance
    }) %>% 
      bind_rows() %>% 
      pivot_longer(-sample_id) %>% 
      mutate(value = sample(value)) %>% 
      group_by(name) %>% 
      slice_sample(n = n) %>% 
      mutate(sample_id = sample(sample_ids)) %>% 
      ungroup() %>% 
      pivot_wider(names_from = name, values_from = value))
      
  }
  
  if (dataset %in% "checkv_all") {
    return(checkv_res_all %>% 
             select(contig_length, checkv_quality) %>% 
             sample_frac(size = 1, replace = FALSE))
  }
  
  if (dataset %in% "checkv_sample_summary") {
    return(tibble("sample_id" = sample_ids,
                  "avg_contig_len" = norm_with_min(n = n, m = 6000, s = 800, min_val = 0),
                  "putative_viral_scaffolds" = norm_with_min(n = n, m = 1675, s = 460, min_val = 0) %>% round(),
                  "n_complete" = norm_with_min(n = n, m = 3.5, s = 2.2, min_val = 0) %>% round(),
                  "n_high_qual" = norm_with_min(n = n, m = 15, s = 6, min_val = 0) %>% round(),
                  "n_med_qual" = norm_with_min(n = n, m = 30, s = 10, min_val = 0) %>% round(),
                  "n_low_qual" = norm_with_min(n = n, m = 860, s = 250, min_val = 0) %>% round(),
                  "n_non_deter" = norm_with_min(n = n, m = 770, s = 211, min_val = 0) %>% round(),
                  "n_contigs" = norm_with_min(n = n, m = 5e4, s = 1.5e5, min_val = 0) %>% round()))
  }
}

dir.create("mock/")
setwd("mock/")

mock <- TRUE

## Mock data based on real
vOTU_stats <- read_tsv("from_real_data/vOTU_stats.tsv") %>% 
  slice_sample(n = nrow(.)) %>% 
  filter(cluster_contigs > 5)

virus_abundance <- read_rds("from_real_data/viral_abundance.Rds") %>% 
  ## Mock dataset only for subset of vOTUs
  select(sample_id, any_of(vOTU_stats %>% pull(vOTU)))

checkv_res_all <- read_rds("from_real_data/checkV_results_all.Rds")

meta_variables <- read_rds("from_real_data/metadata_variables.Rds") 

amgs_cat_per_vOTU <- read_tsv("from_real_data/amg_cat_by_vOTU_mid_min75.tsv") %>% 
  filter(vOTU %in% vOTU_stats$vOTU)

## Mock individual data
sample_data <- create_mock_dataset(dataset = "sample_info", n = 1000)
meta_dat <- create_mock_dataset("metadata", n = 1000, subject_ids = sample_data$deltaker_id)
screening_data <- create_mock_dataset("screening", n = 1000, subject_ids = sample_data$deltaker_id)
virus_abundance <- create_mock_dataset("abundance", n = 1000, sample_ids = sample_data$sample_id)

checkv_res_all <- create_mock_dataset(dataset = "checkv_all")
checkV_res_all_by_sample <- create_mock_dataset(dataset = "checkv_sample_summary", n = 1000, sample_ids = sample_data$sample_id)

dir.create("data/participant_data/", recursive = TRUE)
dir.create("data/vOTU_stats/", recursive = TRUE)
dir.create("data/amgs/", recursive = TRUE)
dir.create("data/abundance_tables/", recursive = TRUE)

## write files
sample_data %>% 
  write_rds("data/participant_data/sample_meta.Rds")
meta_dat %>% 
  write_rds("data/participant_data/metadata_selected_variables.Rds")
meta_variables %>% 
  write_rds("data/participant_data/metadata_variables.Rds")
screening_data %>% 
  write_rds("data/participant_data/screening_data.Rds")
virus_abundance %>% 
  write_rds("data/abundance_tables/viral_abundance.Rds")

vOTU_stats %>% 
  write_tsv("data/vOTU_stats/vOTU_stats.tsv")
amgs_cat_per_vOTU %>% 
  write_tsv("data/amgs/amg_cat_by_vOTU_mid_min75.tsv")

checkv_res_all %>% 
  write_rds("data/abundance_tables/checkV_results_all.Rds")
checkV_res_all_by_sample %>% 
  write_rds("data/abundance_tables/checkV_results_by_sample.Rds")

## Make output directories
dir.create("data/diff_abund/", recursive = TRUE)
dir.create("data/diff_abund/no_crc/", recursive = TRUE)
dir.create("tables/vOTU_based_stats/", recursive = TRUE)
dir.create("data/diversity/alpha_div/", recursive = TRUE)
dir.create("data/diversity/beta_div/", recursive = TRUE)
dir.create("figures/contig_stats/", recursive = TRUE)
dir.create("figures/misc/", recursive = TRUE)
dir.create("figures/vOTU_summary/", recursive = TRUE)
dir.create("figures/diversity_composition/", recursive = TRUE)

## Requires MaAsLin2
source("../analyses/scripts/Maaslin2_analysis.R")

source("../analyses/scripts/for_manuscript/analyses.R")