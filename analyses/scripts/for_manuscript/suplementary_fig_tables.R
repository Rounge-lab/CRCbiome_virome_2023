
# environment -------------------------------------------------------------

setwd("/ess/p1068/data/durable/007-f_smei/001-trro/CRCbiome/papers/Virome/")
source("analyses/scripts/utils/utils.R")

library(tidyverse)
library(vegan)
library(ape)


# Load data ---------------------------------------------------------------

vOTU_stats <- read_tsv("data/vOTU_stats/vOTU_stats.tsv")
vOTU_taxonomy <- read_csv("data/mid/vcontact/results_vcontact2_graph.csv") %>% 
  select(-1) %>% 
  dplyr:: rename(vOTU = Scaffold)
amgs_cat_per_vOTU <- read_tsv("data/amgs/amg_cat_by_vOTU_mid_min75.tsv")
vOTU_by_taxonomy_lifecycle <- 
  read_rds("tables/vOTU_based_stats/vOTU_by_taxonomy_lifecycle.Rds")
fasta_file <-readDNAStringSet("data/mid/dereplication/repr_viral_seqs.fasta")
old_new_names <- read_tsv("data/mid/dereplication/old_to_new_ids.tsv")
clusters <- read_tsv("data/mid/dereplication/galah_clusters.tsv", col_names = FALSE)  
virus_abundance <-read_rds("data/abundance_tables/viral_abundance.Rds")
dramv_amg <- read_tsv("data/mid/dramv/distillate/amg_summary.tsv")%>%
  mutate(vOTU = str_remove(scaffold, "__full-cat.*"))
diff_abund_table <- read_tsv("data/diff_abund/maaslin2_vOTUs_amgs_vs_life_demo_clin.tsv")
meta_dat <- read_rds("data/participant_data/metadata_selected_variables.Rds")
meta_variables <- read_rds("data/participant_data/metadata_variables.Rds")

#Table_1 with tbl_summary function ---------------------------------------------------------------


Virome_Variables <-meta_dat  %>% select(-deltaker_id)

Virome_table1 <- tbl_summary(Virome_Variables,
                             missing = "no",
                             digits = list(
                               gtsummary::all_continuous()  ~ 1, 
                               gtsummary::all_categorical() ~ 1)) #%>% 
#add_overall() %>% 
#add_p(test  = list(
#gtsummary::all_continuous()  ~ "kruskal.test", 
#gtsummary::all_categorical() ~ "chisq.test"))

Virome_table1 %>%
  as_flex_table() %>% 
  flextable::save_as_docx(path = "/tsd/p1068/data/durable/007-f_smei/001-trro/CRCbiome/papers/Virome/tables/table1_gtsummaaary.docx")


#vOTU_Subset >5samples---------------------------------------------------------------

number_of_samples_per_vOTU <- 
  clusters %>% 
  dplyr::rename(old_id = X1) %>% 
  mutate(old_id = str_replace(old_id, "data/mid/dereplication/split/", ""),
         old_id = str_replace(old_id, ".fna", "")) %>% 
  mutate(sample_name = str_extract(X2, "S-[:digit:]*")) %>% 
  left_join(old_new_names) %>% 
  group_by(new_id) %>% 
  summarize(n_samples = length(unique(sample_name)))


vOTU_stats_subset <- 
  number_of_samples_per_vOTU %>% 
  filter(n_samples > 5) %>% 
  dplyr::rename(vOTU = new_id) %>% 
  left_join(vOTU_stats %>% select(genome_length, cluster_contigs,completeness_method,n_samples, 
                                  cluster_proviruses, cluster_not_proviruses,
                                  checkv_quality, miuvig_quality,completeness,
                                  genes, annotated_genes, fraction_annotated, Class,Family, 
                                  BaltimoreGroup, Host, vOTU)) %>%
  left_join(virus_ab_prev, by = "vOTU") 

write_tsv(vOTU_stats_subset, "Export/vOTU_stats_subset.tsv")

 ## AMG

amg_vOTU<-
  number_of_samples_per_vOTU %>%
  #filter(n_samples > 5) %>% 
  dplyr::rename(vOTU = new_id)%>%
  left_join (dramv_amg%>%  select (gene,gene_id, gene_description, 
                                   category, header, subheader, module,vOTU)) %>%
  select(!n_samples) %>%
  filter(!if_all(gene_id, is.na)) %>%
  dplyr:: rename(gene_name = gene)


amg_vOTU%>%
  write_tsv("Export/amg_vOTU.tsv")


#vOTUs_diversity_index---------------------------------------------------------------

beta_dist <- 
  virus_abundance %>% 
  as.data.frame() %>% 
  column_to_rownames(var = "sample_id") %>% 
  as.matrix() %>% 
  vegan::vegdist()

print("mean and sd bray curtis dissimilarity")

vOTUs_diversity_index <- 
  virus_abundance %>% 
  pivot_longer(-sample_id, names_to = "vOTU", values_to = "abundance") %>% 
  group_by(sample_id) %>% 
  summarize(observed = sum(abundance>0),
            shannon_diversity = vegan::diversity(abundance, index = "shannon"),
            invsimpson_diversity = vegan::diversity(abundance, index = "invsimpson")) %>%
  mutate_at(vars(-sample_id), list(~signif(., 3)))  %>% 
  ungroup() %>% 
  pivot_longer(-sample_id, names_to = "index", values_to = "diversity") %>% 
  group_by(index) %>% 
  summarize(mean = mean(diversity),
            sd = sd(diversity),
            median = median(diversity),
            min = min(diversity),
            max = max(diversity)) %>% 
  ungroup() %>% 
  bind_rows(beta_dist %>% 
              as.matrix() %>% 
              as.data.frame() %>% 
              rownames_to_column("sample_id_1") %>% 
              pivot_longer(-sample_id_1, names_to = "sample_id_2", values_to = "dissimilarity") %>% 
              ## upper tri
              filter(as.integer(as.factor(sample_id_1)) > as.integer(as.factor(sample_id_2))) %>%
              summarize(mean = mean(dissimilarity),
                        sd = sd(dissimilarity),
                        min = min(dissimilarity),
                        max = max(dissimilarity),
                        median = median(dissimilarity)) %>% 
              mutate(index = "Bray-Curtis dissimilarity")) %>% 
  mutate(index = case_when(index %in% "invsimpson_diversity" ~ "Inverse Simpson",
                           index %in% "observed" ~ "Observed",
                           index %in% "shannon_diversity" ~ "Shannon",
                           TRUE ~ index) %>% 
           factor(levels = c("Observed", "Shannon", "Inverse Simpson", "Bray-Curtis dissimilarity"))) %>% 
  arrange(index)

vOTUs_diversity_index %>% 
  write_tsv("Export/diversity_summary.tsv")


#Fasta_file Sequences--------------------------------------------------------------

id_subset <- 
  vOTU_stats_subset %>% 
  pull(vOTU)

test <- id_subset
vOTUs_sequences.fasta <- fasta_file[test]
writeXStringSet(vOTUs_sequences.fasta, "data/Export/vOTUs_sequences.fasta")

## Only 1 specific sequence

fasta_file3 <-readDNAStringSet("data/mid/dereplication/repr_viral_seqs.fasta")
subset_list <- grep("vOTU09534", names(fasta_file3), value = TRUE)
vOTU09534_sequence <- fasta_file3[subset_list]
writeXStringSet(vOTU09534_sequence, "data/Export/vOTU09534_sequence")

 ## annotations
fasta_file2 <-readDNAStringSet("data/mid/dramv/genes.fna")
subset_list <- grep("vOTU05693", names(fasta_file2), value = TRUE)
vOTU05693_genes.fasta <- fasta_file2[subset_list]
writeXStringSet(vOTU05693_genes.fasta, "data/Export/vOTU05693_genes.fasta")

# Diff_Abundance table -------------------------------------------------------------

diff_abund_vOTU<-
  diff_abund_table %>%
  left_join(meta_variables %>% 
              select(var_name, variable = var_id),
            by = "variable") %>%
  select(!variable) %>%
  dplyr::rename(variable = var_name) %>%
  filter(dataset == "vOTUs") %>%
  filter(any(qval <= 0.05)) %>%
  select(!metadata) %>%
  dplyr::rename(vOTU = feature) %>%
  mutate(value = if_else(is.na(value), "Low", value)) %>%
  mutate(value = if_else(value == "b_high", "High", value))

diff_abund_vOTU%>% 
  write_tsv("data/Export/Tables/diff_abund_vOTU.tsv")


## Number of vOTUs - Chao

family_summary_table <-
  vOTU_stats %>% 
  group_by(Family) %>% 
  summarize(n_vOTUs = n(),
            n_genomes = sum(cluster_contigs),
            chao_n_genomes = vegan::estimateR(cluster_contigs, index = "ACE") %>% enframe() %>% filter(name %in% "S.chao1") %>% pull(value) %>% signif(3),
            genomes_vOTU_ratio = (n_genomes/n_vOTUs) %>% signif(3))%>% 
 mutate(chao1_fraction = (n_vOTUs/chao_n_genomes) %>% signif(3)) %>% 
  filter(n_vOTUs > 20) 

write_tsv(family_summary_table, "Export/family_summary_table.tsv")



