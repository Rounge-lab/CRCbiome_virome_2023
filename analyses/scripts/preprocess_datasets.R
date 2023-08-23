
# environment -------------------------------------------------------------

setwd("/ess/p1068/data/durable/007-f_smei/001-trro/CRCbiome/papers/Virome/")

library(tidyverse)

# load data ---------------------------------------------------------------

sample_meta <- read_rds("data/participant_data/unprocessed/data_by_sample.Rds")


screening_data <- 
  read_rds("data/participant_data/unprocessed/151222_screening_proc.rds") %>% 
  tibble() %>% 
  mutate(age_cat = cut(age_invitation, breaks = c(50,60,70,80), labels = c("50-60", "60-70", "70-80")))

lifestyle_data <- 
  read_rds("data/participant_data/unprocessed/151222_lifestyle_proc.Rds") %>% 
  tibble()

diet_data <- 
  read_rds("data/participant_data/unprocessed/151222_diet_proc.rds") %>% 
  tibble()


## WCRF - Already excluded those with low quality questionnaires. Also excluded stage IV cancers (n = 2)
diet_wcrf <- read_csv2("data/participant_data/unprocessed/151222_wcrf_proc_CRCbiome.csv")

## energy-adjusted foodgroups
diet_ea <- read_csv2("data/participant_data/unprocessed/151222_nutrients_inc_suppl_residuals.csv")



## Virus abundance - mid, min75cov
virus_abundance <- 
  read_tsv("data/mid/abundance/coverage_table_min75.tsv") %>% 
  pivot_longer(-ID, names_to = "sample_id") %>% 
  mutate(sample_id = gsub("-", "_", sample_id)) %>% 
  inner_join(sample_meta %>% filter(Total_Bases_QC_ATLAS >= 1e9, Prøvetype %in% "Baseline") %>% select(sample_id, Total_Reads_QC_ATLAS)) %>% 
  group_by(ID) %>% 
  mutate(any_obs = any(value > 0)) %>% 
  filter(any_obs) %>% 
  select(-any_obs) %>%
  ungroup() %>% 
  mutate(value = value/Total_Reads_QC_ATLAS*1e6) %>% 
  select(-Total_Reads_QC_ATLAS) %>% 
  pivot_wider(names_from = "ID", values_from = value)


## Virus taxonomy
virus_taxonomy <- read_csv("data/mid/vcontact/results_vcontact2_graph.csv") %>% 
  select(-1) %>% 
  dplyr::rename(vOTU = Scaffold)


## Virus_contigs

read_tsv("data/mid/stats/viral_contigs_by_sample.tsv") %>% 
  mutate(sample_id = gsub("-", "_", sample_id)) %>% 
  mutate(fraction_viral_contigs = n_contigs_viral/n_contigs_over_threshold,
         fraction_viral_bases = n_bases_viral/n_bases_over_threshold) %>% 
  write_tsv("data/abundance_tables/fraction_viral_contigs.tsv")

## Virus AMGs
## To construct these tables, run this script
source("analyses/scripts/auxiliary/amg_by_vOTU_and_sample.R")
process_amg_data()
rm(process_amg_data)



amg_cat_by_sample <- read_tsv("data/amgs/amg_cat_by_sample_mid_min75.tsv") %>% mutate(sample_id = gsub("-", "_", sample_id))
amg_subhead_by_sample <- read_tsv("data/amgs/amg_subhead_by_sample_mid_min75.tsv") %>% mutate(sample_id = gsub("-", "_", sample_id))
amg_module_by_sample <- read_tsv("data/amgs/amg_module_by_sample_mid_min75.tsv") %>% mutate(sample_id = gsub("-", "_", sample_id))



# checkv results ----------------------------------------------------------

cvnames <- c("contig_id", "contig_length", "provirus", "proviral_length", 
             "gene_count", "viral_genes", "host_genes", "checkv_quality", 
             "miuvig_quality", "completeness", "completeness_method", 
             "contamination", "kmer_freq", "warnings")
checkV_results <- read_tsv("data/mid/dereplication/checkV_summary.tsv", col_names = cvnames, skip = 1) %>% 
  mutate(sample_id = str_extract(contig_id, "^S-[:digit:]*"))

checkv_res_all <- read_tsv("data/abundance_tables/checkV_results_all.tsv") %>% 
  filter(sample_id %in% gsub("_", "-", sample_meta$sample_id))



# vOTU stats --------------------------------------------------------------

old_to_new_ids <- read_tsv("data/mid/dereplication/old_to_new_ids.tsv")
derep_clusters <- read_tsv("data/mid/dereplication/galah_clusters.tsv", col_names = FALSE)

number_of_samples_per_vOTU <- 
  derep_clusters %>% 
  rename(old_id = X1) %>% 
  mutate(old_id = str_replace(old_id, "data/mid/dereplication/split/", ""),
         old_id = str_replace(old_id, ".fna", "")) %>% 
  mutate(sample_name = str_extract(X2, "S-[:digit:]*")) %>% 
  left_join(old_to_new_ids %>% rename(vOTU = new_id)) %>% 
  group_by(vOTU) %>% 
  summarize(n_samples = length(unique(sample_name)))


## Get some information from dramv
dramv_summary <- 
  read_tsv("data/mid/dramv/distillate/vMAG_stats.tsv") %>% 
  rename(scaffold = 1) %>% 
  mutate(vOTU = str_extract(scaffold, "^vOTU[:digit:]*")) %>% 
  select(vOTU, everything())

## Get some information from virsorter second round
virsorter_round2_scores <- 
  read_tsv("data/mid/dereplication/virsorter/final-viral-score.tsv") %>% 
  mutate(vOTU = str_extract(seqname, "^vOTU[:digit:]*")) %>% 
  select(vOTU, everything())

## Get actual fasta seq lengths
# fasta_seq_lengths <- read_delim("fasta_len.txt", delim = " ", col_names = FALSE) %>% 
#   rename(vOTU = 1, seq_length = 2) %>% 
#   mutate(vOTU = str_replace(vOTU, ">", ""))

## Get number of proteins
n_proteins <- 
  read_csv("data/mid/vcontact/merged_df.csv") %>% 
  filter(str_detect(contig_id, "vOTU[:digit:]*")) %>% 
  rename(vOTU = contig_id, genes = proteins) %>% 
  select(vOTU, genes)

## There is some difference between genes used for dramv annotation and vcontact
# dramv_summary %>% 
#   left_join(n_proteins) %>% 
#   select(vOTU, `Gene count`, genes) %>% 
#   filter(`Gene count` != genes) # %>% 
  # View()
  # summarize(total_diff = sum(abs(`Gene count`-genes)))

# checkV_results %>% 
#   inner_join(old_to_new_ids %>% rename(contig_id = old_id, vOTU = new_id)) %>% 
#   left_join(virsorter_round2_scores %>% select(vOTU, length)) %>% 
#   select(vOTU, contig_length, provirus, proviral_length, length, gene_count, viral_genes) %>% 
#   left_join(fasta_seq_lengths) %>% 
#   left_join(n_proteins) %>% 
#   View()


# clusters_votu <-
vOTU_stats <-
  checkV_results %>% 
  inner_join(old_to_new_ids %>% rename(contig_id = old_id, vOTU = new_id)) %>% 
  select(-c(contig_id, sample_id)) %>% 
  mutate(genome_length = case_when(provirus %in% "Yes" ~ proviral_length,
                                   provirus %in% "No" ~ contig_length)) %>% 
  left_join(n_proteins) %>% 
  select(vOTU, genome_length, genes, everything(), -c(provirus, contig_length, 
                                               proviral_length, gene_count, 
                                               viral_genes, host_genes, contamination,
                                               kmer_freq, warnings)) %>% 
  left_join(derep_clusters %>% 
              mutate(old_id = str_replace(X1, "data/mid/dereplication/split/", ""),
                     old_id = str_replace(old_id, ".fna", ""),
                     contig_id = str_replace(X2, "data/mid/dereplication/split/", ""),
                     contig_id = str_replace(contig_id, ".fna", "")) %>% 
              left_join(old_to_new_ids %>% rename(vOTU = new_id)) %>% 
              left_join(checkV_results %>% select(contig_id, provirus)) %>% 
              count(vOTU, provirus) %>% 
              pivot_wider(names_from = provirus, values_from = n, values_fill = 0) %>% 
              mutate(cluster_contigs = Yes + No) %>% 
              rename(cluster_proviruses = Yes,
                     cluster_not_proviruses = No) %>% 
              select(vOTU, cluster_contigs, cluster_proviruses, cluster_not_proviruses)) %>% 
  left_join(number_of_samples_per_vOTU) %>% 
  left_join(virus_abundance %>% 
              pivot_longer(-sample_id, names_to = "vOTU", values_to = "abundance") %>% 
              group_by(sample_id) %>% 
              mutate(rel_abundance = abundance/sum(abundance)) %>% 
              ungroup() %>% 
              group_by(vOTU) %>% 
              summarize(prevalence = sum(abundance > 0)/n(),
                        mean_prevalent_abundance = mean(abundance[abundance > 0]),
                        mean_prevalent_rel_abund = mean(rel_abundance[abundance > 0])) %>% 
              ungroup())

# vOTU_stats %>% 
#   filter(cluster_contigs > 2) %>% 
#   ggplot(aes(y = mean_prevalent_abundance, x = prevalence)) +
#   geom_point() +
#   scale_x_continuous(trans = "log10") +
#   scale_y_continuous(trans = "log10")
#   
# vOTU_stats %>% 
#   filter(cluster_contigs > 2) %>% 
#   ggplot(aes(y = mean_prevalent_abundance, x = mean_prevalent_rel_abund)) +
#   geom_point() +
#   scale_x_continuous(trans = "log10") +
#   scale_y_continuous(trans = "log10")


# summarize lifecycle -----------------------------------------------------

ref_prop_provirus <- 
  vOTU_stats %>%
  summarize(prop = sum(cluster_proviruses)/(sum(cluster_proviruses)+sum(cluster_not_proviruses))) %>% 
  pull(prop)

## Script with function to test preponderance of lifecycle state of contigs per vOTU. 
## Takes background rate as argument. It seems reasonable to use the overall observed fraction as reference here.

source("analyses/scripts/auxiliary/categorize_lifecycle.R")
categorize_lifecycle(ref_prop_provirus)


# sample_meta -------------------------------------------------------------

sample_meta %>% 
  filter(Total_Bases_QC_ATLAS >= 1e9,
         Prøvetype %in% "Baseline") %>% 
  # filter(train_test %in% "train") %>% 
  # select(-train_test) %>% 
  write_rds("data/participant_data/sample_meta.Rds")


# screening_data ----------------------------------------------------------

screening_data %>% 
  inner_join(sample_meta %>% filter(Total_Bases_QC_ATLAS >= 1e9, Prøvetype %in% "Baseline") %>% select(deltaker_id)) %>%
  write_rds("data/participant_data/screening_data.Rds")

# lifestyle ---------------------------------------------------------------

lifestyle_data %>%
  left_join(screening_data %>% select(id, deltaker_id)) %>% 
  inner_join(sample_meta %>% filter(Total_Bases_QC_ATLAS >= 1e9, Prøvetype %in% "Baseline") %>% select(deltaker_id)) %>% 
  write_rds("data/participant_data/lifestyle_data.Rds")

# diet --------------------------------------------------------------------

## Summarize exclusion
diet_data %>% 
  left_join(screening_data %>% select(id, deltaker_id)) %>% 
  inner_join(sample_meta %>% filter(Total_Bases_QC_ATLAS >= 1e9, Prøvetype %in% "Baseline") %>% select(deltaker_id)) %>% 
  filter(exclude_nocolo %in% 0,
         exclude_reserved %in% 0) %>% ## FFQ for those meeting general inclusion criteria = 1562 participants
  summarize(FFQ_qual_low = sum(exclude_FFQquality_low %in% 1),
            energy_low = sum(exclude_energylow %in% 1),
            energy_high = sum(exclude_energyhigh %in% 1),
            remaining = sum(exclude_FFQquality_low %in% 0 & exclude_energylow %in% 0 & exclude_energyhigh %in% 0)) ## Number excluded by reason for participants with FFQs

## Exclude those not meeting criteria:
diet_data_proc <- 
  diet_data %>%
  left_join(screening_data %>% select(id, deltaker_id)) %>% 
  inner_join(sample_meta %>% filter(Total_Bases_QC_ATLAS >= 1e9, Prøvetype %in% "Baseline") %>% select(deltaker_id)) %>% 
  filter(exclude_nocolo %in% 0,
         exclude_reserved %in% 0,
         exclude_FFQquality_low %in% 0,
         exclude_energyhigh %in% 0,
         exclude_energylow %in% 0) %>% 
  # left_join(screening_data %>% select(deltaker_id, id)) %>%
  mutate(BMI_cat = case_when(!is.na(BMI_cat3) ~ as.character(BMI_cat3),
                             is.na(BMI_cat3) ~ "Missing")) %>% 
  mutate(BMI_cat = factor(BMI_cat, levels = c("Normal weight", "Overweight", "Obese", "Missing")))

diet_data_proc %>% 
  write_rds("data/participant_data/diet_data.Rds")

# diet wcrf ---------------------------------------------------------------

diet_wcrf %>% 
  left_join(screening_data %>% select(id, deltaker_id)) %>% 
  inner_join(sample_meta %>% filter(Total_Bases_QC_ATLAS >= 1e9, Prøvetype %in% "Baseline") %>% select(deltaker_id)) %>% 
  write_rds("data/participant_data/diet_wcrf.Rds")

# diet energy adjusted ----------------------------------------------------

diet_data_proc %>% 
  select(id) %>% 
  left_join(diet_ea) %>% 
  write_rds("data/participant_data/diet_ea.Rds")

# Virus abundance ---------------------------------------------------------

## Viruses

virus_abundance %>% 
  write_rds("data/abundance_tables/viral_abundance.Rds")

n_vOTUs <- virus_abundance %>% pivot_longer(-sample_id) %>% group_by(sample_id) %>% summarize(n_vOTU = sum(value > 0)) %>% ungroup()



# vOTU stats --------------------------------------------------------------

dramv_annotations <- read_tsv("data/mid/dramv/annotations.tsv") %>% 
  mutate(vOTU = str_extract(scaffold, "^vOTU[:digit:]*")) %>% 
  rename(gene = 1)

## Check whether dramv annotations file contains all genes, or only a subset
dramv_annotations %>% 
  group_by(vOTU) %>% 
  summarize(max_gene_pos = max(gene_position),
            n_genes = n()) %>% 
  ungroup() %>% 
  mutate(same_n = max_gene_pos == n_genes) %>% 
  count(same_n)

fraction_annotated <-
  dramv_annotations %>% 
  mutate(any_annotation = if_any(contains("hit"), .fns = function(x) !is.na(x))) %>% 
  group_by(vOTU) %>% 
  summarize(n_genes = n(),
            n_annotated = sum(any_annotation),
            fraction_annotated = sum(any_annotation)/n()) %>% 
  ungroup()

vOTU_stats %>% 
  left_join(virus_taxonomy %>% select(vOTU, Class, Order, Family, Subfamily, Genus, BaltimoreGroup, Host)) %>% 
  left_join(fraction_annotated %>% select(vOTU, annotated_genes = n_annotated, fraction_annotated)) %>% 
  #left_join(integrase_cat) %>% 
  write_tsv("data/vOTU_stats/vOTU_stats.tsv")

# amgs --------------------------------------------------------------------

amg_cat_frac_by_sample <- 
  amg_cat_by_sample %>% 
  pivot_longer(-sample_id, names_to = "amg_cat", values_to = "n_amg_OTUs") %>% 
  left_join(n_vOTUs) %>% 
  mutate(frac_amg_vOTUs = n_amg_OTUs/n_vOTU) %>% 
  select(-c(n_amg_OTUs, n_vOTU)) %>% 
  pivot_wider(names_from = amg_cat, values_from = frac_amg_vOTUs)

amg_subhead_frac_by_sample <- 
  amg_subhead_by_sample %>% 
  pivot_longer(-sample_id, names_to = "amg_cat", values_to = "n_amg_OTUs") %>% 
  left_join(n_vOTUs) %>% 
  mutate(frac_amg_vOTUs = n_amg_OTUs/n_vOTU) %>% 
  select(-c(n_amg_OTUs, n_vOTU)) %>% 
  pivot_wider(names_from = amg_cat, values_from = frac_amg_vOTUs)

amg_module_frac_by_sample <- 
  amg_module_by_sample %>% 
  pivot_longer(-sample_id, names_to = "amg_cat", values_to = "n_amg_OTUs") %>% 
  left_join(n_vOTUs) %>% 
  mutate(frac_amg_vOTUs = n_amg_OTUs/n_vOTU) %>% 
  select(-c(n_amg_OTUs, n_vOTU)) %>% 
  pivot_wider(names_from = amg_cat, values_from = frac_amg_vOTUs)

amg_cat_frac_by_sample %>% 
  write_tsv("data/amgs/amg_cat_frac_by_sample_mid_min75.tsv")
amg_subhead_frac_by_sample %>% 
  write_tsv("data/amgs/amg_subhead_frac_by_sample_mid_min75.tsv")
amg_module_frac_by_sample %>% 
  write_tsv("data/amgs/amg_module_frac_by_sample_mid_min75.tsv")


# lifestyle and demography variables --------------------------------------

variables_tmp <- bind_rows(
  c("Sex" = "kjonn", 
    "Age" = "age_cat", 
    "Region" = "senter",
    "BMI" = "BMI",
    "Physical activity" = "PhysAct_Score",
    "Smoking" = "Smoking", 
    "Snus" = "Snus",
    "Education" = "Utdanning", 
    "Marital status" = "Sivilstatus_cat2",
    "Working status" = "Arbeid_lump",
    "Nationality" = "Nasj_cat2") %>% 
    enframe(name = "var_name", value = "var_id") %>% 
    mutate(dataset = case_when(var_name %in% c("BMI", "Physical activity", "Smoking", "Snus") ~ "lifestyle", 
                               TRUE ~ "demography")),
  
  c("Energy (kcal/day)" = "Energi_kcal", 
    "Protein (E%)" = "Prot_energi", 
    "Carbohydrates (E%)" = "Karboh_energi", 
    "Added sugar (E%)" = "Sukker_energi", 
    "Fiber (g/day)" = "Fiber", 
    "Fat (E%)" = "Fett_energi", 
    "SFA (E%)" = "Mettet_energi", 
    "MUFA (E%)" = "C_enum_energi",
    "PUFA (E%)" = "C_flerum_energi", 
    "TFA (E%)" = "Trans_u_energi", 
    "Alcohol (g/day)" = "Alko",
    "HLI" = "wcrf_index_main") %>% 
    enframe(name = "var_name", value = "var_id") %>% 
    mutate(dataset = case_when(var_name %in% "WCRF" ~ "diet and lifestyle",
                               TRUE~"diet")),
  
  c("Antibiotics use" = "Antibiotics",
    "Antacid use" = "Antacids") %>% 
    enframe(name = "var_name", value = "var_id") %>% 
    mutate(dataset = "medication")) #%>%
  
meta_dat <- 
  read_rds("data/participant_data/sample_meta.Rds") %>% 
  select(deltaker_id) %>% 
  left_join(read_rds("data/participant_data/screening_data.Rds") %>% select(any_of(variables_tmp %>% pull(var_id)), deltaker_id, id), by = "deltaker_id") %>% 
  left_join(read_rds("data/participant_data/diet_data.Rds") %>% select(any_of(variables_tmp %>% pull(var_id)), id), by = "id") %>% 
  left_join(read_rds("data/participant_data/diet_wcrf.Rds") %>% select(any_of(variables_tmp %>% pull(var_id)), id), by = "id") %>% 
  left_join(read_rds("data/participant_data/lifestyle_data.Rds") %>% select(any_of(variables_tmp %>% pull(var_id)), id), by = "id") %>% 
  select(-id)

variables <- 
  variables_tmp %>% 
  mutate(ref = case_when(var_id %in% "kjonn" ~ "Female",
                         var_id %in% "age_cat" ~ "50-60",
                         var_id %in% "senter" ~ "Moss",
                         var_id %in% "Antibiotics" ~ "No",
                         var_id %in% "Antacids" ~ "No",
                         var_id %in% "Smoking" ~ "Non smoker",
                         var_id %in% "Snus" ~ "Non snuser",
                         var_id %in% "Utdanning" ~ "Primary school",
                         var_id %in% "Sivilstatus_cat2" ~ "Married/cohabiting",
                         var_id %in% "Arbeid_lump" ~ "Employed",
                         var_id %in% "Nasj_cat2" ~ "Native")) %>% 
  mutate(var_id = factor(var_id, levels = names(meta_dat)[ names(meta_dat) %in% var_id]))


variables %>% 
  write_rds("data/participant_data/metadata_variables.Rds")
meta_dat %>% 
  write_rds("data/participant_data/metadata_selected_variables.Rds")
meta_dat %>% 
  mutate(across(where(is.double), .fns = function(x) cat_func(x) %>% str_replace("negative", "low") %>% str_replace("positive", "high") %>% factor(levels = c("low", "mid", "high")))) %>% 
  write_rds("data/participant_data/metadata_selected_variables_cat.Rds")
