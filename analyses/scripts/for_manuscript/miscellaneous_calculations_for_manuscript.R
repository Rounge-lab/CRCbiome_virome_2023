

## Number of raw paired end reads
sample_data %>% 
  pull(Total_Reads_raw_ATLAS) %>% 
  sum()

## Number of paired end reads
sample_data %>% 
  pull(Total_Reads_QC_ATLAS) %>% 
  sum()

## Number of paired end reads per sample
sample_data %>% 
  pull(Total_Reads_QC_ATLAS) %>% 
  summary()

## IQR
sample_data %>% 
  pull(Total_Reads_QC_ATLAS) %>% 
  IQR()

## Storage time for samples
sample_data %>% 
  pull(sampling_to_extraction) %>% 
  summary()

sample_data %>% 
  summarize(from_sampling_to_extraction = difftime(ekstraksjon_dato, collected_date, units = "weeks") %>% as.integer())

## Storage time for samples
sample_data %>% 
  pull(qubit_total_dna) %>% 
  summary()

## Number of genomes by quality
checkv_res_all %>% 
  count(checkv_quality)

get_some_stats <- function() {
  diversity_summary <- read_tsv("data/Export/diversity_summary.tsv")
  print("mean inverse simpson")
  print(diversity_summary %>% filter(index %in% "invsimpson_diversity") %>% pull(mean))
  print("sd inverse simpson")
  print(diversity_summary %>% filter(index %in% "invsimpson_diversity") %>% pull(sd))
  
  beta_dist <- 
    virus_abundance %>% 
    as.data.frame() %>% 
    column_to_rownames(var = "sample_id") %>% 
    as.matrix() %>% 
    vegan::vegdist()
  print("mean and sd bray curtis dissimilarity")
  beta_dist %>% 
    as.matrix() %>% 
    as.data.frame() %>% 
    rownames_to_column("sample_id_1") %>% 
    pivot_longer(-sample_id_1, names_to = "sample_id_2", values_to = "dissimilarity") %>% 
    ## upper tri
    filter(as.integer(as.factor(sample_id_1)) > as.integer(as.factor(sample_id_2))) %>%
    summarize(mean = mean(dissimilarity),
              sd = sd(dissimilarity),
              min_dissim = min(dissimilarity),
              max_dissim = max(dissimilarity),
              median_dissim = median(dissimilarity)) %>% 
    print()
}
# get_some_stats()

vOTU_stats %>% 
  count(Family) %>% 
  as.data.frame()

# vOTU_stats %>% 
#   filter(Family %in% "Unclassified") %>% 
#   View()

vOTU_taxonomy <- read_csv("data/mid/vcontact/results_vcontact2_graph.csv")

## Number of vOTUs classified by direct and indirect connections
family_summary_table <-
  vOTU_stats %>% 
  group_by(Family) %>% 
  summarize(n_vOTUs = n(),
            n_genomes = sum(cluster_contigs),
            max_genomes = max(cluster_contigs),
            shannon_n_genomes = vegan::diversity(cluster_contigs) %>% signif(3),
            chao_n_genomes = vegan::estimateR(cluster_contigs, index = "ACE") %>% enframe() %>% filter(name %in% "S.chao1") %>% pull(value) %>% signif(3),
            chao_n_genomes_se = vegan::estimateR(cluster_contigs, index = "ACE") %>% enframe() %>% filter(name %in% "se.chao1") %>% pull(value) %>% signif(3),
            ace_n_genomes = vegan::estimateR(cluster_contigs, index = "chao") %>% enframe() %>% filter(name %in% "S.ACE") %>% pull(value) %>% signif(3),
            ace_n_genomes_se = vegan::estimateR(cluster_contigs, index = "chao") %>% enframe() %>% filter(name %in% "se.ACE") %>% pull(value) %>% signif(3),
            genomes_vOTU_ratio = (n_genomes/n_vOTUs) %>% signif(3)) %>% 
  left_join(vOTU_taxonomy %>% 
              mutate(connection_type = str_extract(Level, "^[:alpha:]"),
                     connection_rank = str_extract(Level, "[:digit:]*$") %>% as.integer()) %>% 
              count(connection_type, Family) %>% 
              pivot_wider(names_from = connection_type, values_from = n, values_fill = 0) %>% 
              select(Family, C, N) %>% 
              # filter(C > 0 | N > 0) %>% 
              rename(direct = C, indirect = N) %>% 
              mutate(n_vOTUs = direct + indirect,
                     fraction_direct = (direct/n_vOTUs) %>% signif(3),
                     fraction_indirect = (1-fraction_direct) %>% signif(3)) %>% 
              select(-n_vOTUs)) %>% 
  mutate(Family = fct_reorder(Family, n_vOTUs, .desc = TRUE)) %>% 
  arrange(Family) 

## How many viral genomes were used as input to dereplication?
vOTU_stats %>% 
  pull(cluster_contigs) %>% 
  sum()

18494/49642

vegan::estimateR(vOTU_stats$cluster_contigs)
sum(family_summary_table$chao_n_genomes)
sum(family_summary_table$ace_n_genomes, na.rm = TRUE)

summary(vOTU_stats)

## What fraction of vOTU diversity was discovered by family?
family_summary_table %>% 
  mutate(chao1_fraction = (n_vOTUs/chao_n_genomes) %>% signif(3)) %>% 
  select(Family, n_vOTUs, chao_n_genomes, chao1_fraction, fraction_direct) %>%
  as.data.frame()

family_summary_table %>% 
  mutate(chao1_fraction = (n_vOTUs/chao_n_genomes) %>% signif(3)) %>% 
  ggplot(aes(x = fraction_direct, y = chao1_fraction, color = n_vOTUs > 20, label = Family)) +
  geom_point() +
  ggrepel::geom_label_repel(data = family_summary_table %>% 
                              mutate(chao1_fraction = (n_vOTUs/chao_n_genomes) %>% signif(3)) %>% 
                              filter(n_vOTUs > 20))



family_summary_table %>% 
  ggplot(aes(x = genomes_vOTU_ratio, fraction_direct, label = Family, color = n_vOTUs>20)) +
  geom_point() +
  ggrepel::geom_label_repel(data = family_summary_table %>% filter(n_vOTUs > 20)) +
  theme_bw() +
  labs(y = "fraction reference connections direct",
       x = "ratio of genomes to vOTUs")



family_summary_table %>% 
  pivot_longer(c(fraction_direct, fraction_indirect, n_vOTUs)) %>% 
  mutate(frac = case_when(str_detect(name, "fraction") ~ TRUE,
                          TRUE ~ FALSE)) %>% 
  ggplot(aes(y = Family, x = value, fill = name)) +
  geom_col() +
  scale_x_continuous(trans = "log10") +
  facet_wrap(~frac, ncol = 2, scales = "free_x")
  

## Number of vOTUs classified
vOTU_taxonomy %>% 
  # count(Family) %>% 
  filter(str_detect(Family, "ae$")) %>% 
  # filter(Family %in% c("Intestiviridae", "Steigviridae", "Suoliviridae", "Crevaviridae")) %>% 
  # summarize(n_crass = sum(n[ Family %in% c("Intestiviridae", "Steigviridae", "Suoliviridae", "Crevaviridae")]))
  mutate(connection_type = str_extract(Level, "^[:alpha:]"),
         connection_rank = str_extract(Level, "[:digit:]*$") %>% as.integer()) %>% 
  summarize(n_classified = sum(connection_type %in% c("C", "N")),
            frac_classified = sum(connection_type %in% c("C", "N"))/n(),
            family_classified = sum(str_detect(Family, "dae")),
            frac_family_classified = sum(str_detect(Family, "dae"))/n())

## Number of vOTUs classified
vOTU_stats %>% 
  count(Family) %>% 
  mutate(frac = round(n/sum(n), 3)) %>% 
  filter(str_detect(Family, "dae")) %>% 
  mutate(frac_sub = round(n/sum(n), 2)) %>% 
  arrange(desc(n)) %>% 
  as.data.frame()


## What is the length of the peduoviridae family genomes?
vOTU_stats %>% 
  filter(Family %in% "Peduoviridae") %>% 
  pull(genome_length) %>% 
  summary()


## The "higher order" group probably contains many of the previous siphoviridae, myoviridae and podoviridae
prev_taxonomy <-
  read_csv("data/prev_vcontact/results_vcontact2_graph.csv") %>% 
  rename(vOTU = Scaffold) %>% 
  select(vOTU, Family, Class) %>% 
  rename_with(.cols = -vOTU, .fn = function(x) paste0("old_", x))

vOTU_stats %>% 
  select(vOTU, Family, Class) %>% 
  left_join(prev_taxonomy) %>% 
  mutate(abolished_taxa = case_when(old_Family %in% c("Siphoviridae", "Myoviridae", "Podoviridae") ~ "abolished taxa",
                                    TRUE ~ "others")) %>% 
  filter(Family %in% c("Unclassified", "n.a.")) %>% 
  count(abolished_taxa, old_Family)
  # count(Class, old_Class)

vOTU_stats %>% 
  select(vOTU, Family, Class) %>% 
  left_join(prev_taxonomy) %>% 
  filter(Family %in% c("Unclassified", "n.a.")) %>% 
  count(Family, old_Family)

vOTU_stats %>% 
  select(vOTU, Family, Class) %>% 
  left_join(prev_taxonomy) %>% 
  mutate(abolished_taxa = case_when(old_Family %in% c("Siphoviridae", "Myoviridae", "Podoviridae") ~ "abolished taxa",
                                    TRUE ~ "others")) %>% 
  count(Family, abolished_taxa) %>% 
  filter(abolished_taxa %in% "abolished taxa")

## fraction lytic
vOTU_stats %>% 
  group_by(Family) %>% 
  summarize(fraction_provirus = sum(cluster_proviruses)/sum(cluster_contigs),
            fraction_lytic = sum(cluster_not_proviruses)/sum(cluster_contigs),
            n_lytic = sum(cluster_not_proviruses),
            n_lysogenic = sum(cluster_proviruses),
            n_tot = n()) %>% 
  filter(n_tot > 20) %>% 
  arrange(desc(fraction_lytic)) %>% 
  as.data.frame()


vOTU_stats %>% 
  left_join(amg_summary %>% 
              mutate(vOTU = str_extract(scaffold, "^vOTU[:digit:]*")) %>% 
              select(vOTU) %>%
              unique() %>% 
              mutate(amg_out = "in_amg_table")) %>% 
  count(amg_out)

vOTU_stats %>% 
  left_join(vOTU_amg_summary %>% 
              mutate(vOTU = str_extract(scaffold, "^vOTU[:digit:]*")) %>% 
              select(vOTU) %>%
              unique() %>% 
              mutate(amg_out = "in_amg_table")) %>% 
  count(amg_out)


  

  # mutate(any_annotation = . %>% select(contains("hit")) %>% apply(1, function(x) any(!is.na(x))))

vOTU_stats %>% 
  left_join(dramv_annotations %>% 
              mutate(vOTU = str_extract(scaffold, "^vOTU[:digit:]*")) %>% 
              select(vOTU) %>%
              count(vOTU) %>% 
              mutate(dramv_out = "in_dramv_table")) %>% 
  count(n) #%>% 
  #View()

checkv_res_all %>% 
  filter(checkv_quality %in% c("Complete", "High-quality", "Medium-quality")) %>% 
  mutate(viral_length = case_when(is.na(proviral_length) ~ contig_length,
                                  !is.na(proviral_length) ~ proviral_length)) %>% 
  arrange(viral_length)
  # mutate(viral_over_1500 = viral_length > 1500) %>% 
  # count(viral_over_1500)

virsorter_for_dramv <- read_tsv("data/mid/dereplication/virsorter/final-viral-score.tsv")


## AMGs

## AMGs detected in what fraction of genomes?
amgs_cat_per_vOTU %>% 
  pivot_longer(-vOTU) %>% 
  filter(!str_detect(name, "Wood")) %>% 
  group_by(vOTU) %>% 
  summarize(any_AMG = any(value)) %>% 
  ungroup() %>% 
  skimr::skim(any_AMG)

## How common in viral families
amgs_cat_per_vOTU %>% 
  pivot_longer(-vOTU) %>% 
  filter(!str_detect(name, "Wood")) %>% 
  group_by(vOTU) %>% 
  summarize(any_AMG = any(value)) %>% 
  ungroup() %>% 
  left_join(vOTU_stats %>% select(vOTU, Family)) %>% 
  group_by(Family) %>% 
  skimr::skim(any_AMG) %>% 
  arrange(desc(logical.mean))

## How common in crassviruses?
amgs_cat_per_vOTU %>% 
  pivot_longer(-vOTU) %>% 
  filter(!str_detect(name, "Wood")) %>% 
  group_by(vOTU) %>% 
  summarize(any_AMG = any(value)) %>% 
  ungroup() %>% 
  left_join(vOTU_stats %>% select(vOTU, Family)) %>% 
  mutate(crassvirus = Family %in% c("Crevaviridae", "Intestiviridae", "Steigviridae", "Suoliviridae")) %>% 
  filter(crassvirus) %>% 
  skimr::skim(any_AMG)

amgs_cat_per_vOTU %>% 
  pivot_longer(-vOTU) %>% 
  filter(!str_detect(name, "Wood")) %>% 
  group_by(name) %>% 
  summarize(fraction_vOTUs = sum(value)/n()) 

amgs_cat_per_vOTU %>% 
  pivot_longer(-vOTU) %>% 
  filter(!str_detect(name, "Wood")) %>% 
  left_join(vOTU_stats %>% select(vOTU, Family)) %>% 
  group_by(name, Family) %>% 
  summarize(fraction_vOTUs = sum(value)/n()) %>% 
  View()

## MISC per family
amgs_cat_per_vOTU %>% 
  pivot_longer(-vOTU) %>% 
  filter(str_detect(name, "MISC")) %>% 
  left_join(vOTU_stats %>% select(vOTU, Family)) %>% 
  group_by(Family) %>% 
  summarize(fraction_vOTUs = (sum(value)/n()) %>% signif(3)) %>% 
  as.data.frame()

## MISC per family
amgs_cat_per_vOTU %>% 
  pivot_longer(-vOTU) %>% 
  filter(str_detect(name, "Nitro")) %>% 
  filter(!str_detect(name, "Wood")) %>% 
  left_join(vOTU_stats %>% select(vOTU, Family)) %>% 
  group_by(Family) %>% 
  summarize(fraction_vOTUs = (sum(value)/n()) %>% signif(3)) %>% 
  as.data.frame()
  
amgs_cat_per_vOTU %>% 
  pivot_longer(-vOTU) %>% 
  filter(!str_detect(name, "Wood")) %>% 
  left_join(vOTU_stats %>% select(vOTU, Family)) %>% 
  filter(Family %in% c("Crevaviridae", "Intestiviridae", "Steigviridae", "Suoliviridae")) %>% 
  group_by(name) %>% 
  summarize(fraction_vOTUs = sum(value)/n())



## Abundance stuff
checkV_res_all_by_sample %>% 
  select(sample_id, n_complete, n_high_qual, n_med_qual) %>% 
  pivot_longer(starts_with("n")) %>% 
  group_by(sample_id) %>% 
  summarize(n_identified = sum(value)) %>% 
  ungroup() %>% 
  pull(n_identified) %>% 
  summary()

virus_abundance %>% 
  pivot_longer(-sample_id) %>% 
  group_by(sample_id) %>% 
  summarize(n_identified = sum(value)) %>% 
  ungroup() %>% 
  pull(n_identified) %>% 
  summary()
  
## How many of the vOTUs are detected at >= 1% (this is the threshold we use for differential abundance, for example)
virus_abundance %>% 
  pivot_longer(-sample_id) %>% 
  group_by(name) %>% 
  summarize(frac_pres = sum(value > 0)/n()) %>% 
  ungroup() %>% 
  summarize(n_over_1 = sum(frac_pres>0.01),
            frac_over_1 = sum(frac_pres>0.01)/n())

viral_fraction_by_family <- 
  virus_abundance %>% 
  pivot_longer(-sample_id, names_to = "vOTU") %>% 
  group_by(sample_id) %>% 
  mutate(value = value/sum(value)) %>% 
  left_join(vOTU_stats %>% select(vOTU, Family)) %>% 
  group_by(sample_id, Family) %>% 
  summarize(fraction_by_family = sum(value)) %>% 
  ungroup()


## Fraction of total abundance attributable to viruses with taxonomic annotation
viral_fraction_by_family %>% 
  filter(Family != "n.a.") %>% 
  mutate(crassvirus = Family %in% c("Crevaviridae", "Intestiviridae", "Steigviridae", "Suoliviridae")) %>% 
  group_by(sample_id) %>% 
  summarize(known_frac = sum(fraction_by_family),
            crass_frac = sum(fraction_by_family[ crassvirus]),
            any_crass = sum(fraction_by_family[ crassvirus])>0) %>% 
  ungroup() %>% 
  mutate(crass_frac_of_known = crass_frac/known_frac) %>% 
  skimr::skim(known_frac, crass_frac, any_crass)

chisq.test(100*(rbind(c(1,2,4,2,4),
                      c(2,1,2,1,5))))

## Prevalence of viral families
family_prevalence <- 
  virus_abundance %>% 
  pivot_longer(-sample_id, names_to = "vOTU", values_to = "abundance") %>% 
  left_join(vOTU_stats %>% select(vOTU, Family)) %>% 
  group_by(sample_id, Family) %>% 
  summarize(fam_pres = any(abundance > 0)) %>% 
  group_by(Family) %>% 
  summarize(fam_prev = sum(fam_pres)/n()) %>% 
  arrange(desc(fam_prev))


vOTU_stats %>% 
  count(Order, Family)

## Prevalence of crassviruses
crass_prev <- 
  virus_abundance %>% 
  pivot_longer(-sample_id, names_to = "vOTU", values_to = "abundance") %>% 
  inner_join(vOTU_stats %>% filter(Order %in% "Crassvirales") %>% select(vOTU)) %>% 
  group_by(sample_id) %>% 
  summarize(crass_pres = any(abundance > 0)) %>% 
  ungroup() %>% 
  summarize(crass_prev = sum(crass_pres)/n())

## Crass abundance
virus_abundance %>% 
  pivot_longer(-sample_id, names_to = "vOTU", values_to = "abundance") %>% 
  group_by(sample_id) %>% 
  mutate(rel_abund = abundance/sum(abundance)) %>% 
  inner_join(vOTU_stats %>% filter(Order %in% "Crassvirales") %>% select(vOTU)) %>% 
  summarize(crass_abundance = sum(rel_abund)) %>% 
  ungroup() %>% 
  skimr::skim(crass_abundance)

## Microviridae abundance
virus_abundance %>% 
  pivot_longer(-sample_id, names_to = "vOTU", values_to = "abundance") %>% 
  group_by(sample_id) %>% 
  mutate(rel_abund = abundance/sum(abundance)) %>% 
  inner_join(vOTU_stats %>% filter(Family %in% "Microviridae") %>% select(vOTU)) %>% 
  summarize(micro_abundance = sum(rel_abund)) %>% 
  ungroup() %>% 
  skimr::skim(micro_abundance)

## Are there overall differences in viral family abundance according to center?
viral_fraction_by_family %>% 
  left_join(sample_data %>% left_join(meta_dat) %>% select(sample_id, senter)) %>% 
  ggplot(aes(x = fraction_by_family, y = senter, fill = senter)) +
  geom_boxplot() +
  facet_wrap(~Family, scales = "free_x")

abundance_diff_by_fam_and_host_var <-
  viral_fraction_by_family %>% 
  left_join(sample_data %>% 
              select(sample_id, deltaker_id) %>% 
              left_join(meta_dat) %>% 
              select(-deltaker_id) %>% 
              mutate(across(where(is.double), .fns = function(x) cat_func(x) %>% str_replace("negative", "a_low") %>% str_replace("positive", "b_high"))) %>% 
              pivot_longer(-sample_id)) %>% 
  group_by(Family, name) %>% 
  group_split() %>% 
  lapply(function(x) {
    x %>% 
      filter(!value %in% c("Unknown", "Missing")) %>% 
      kruskal.test(formula = fraction_by_family~value, data = .) %>% 
      tidy() %>% 
      mutate(Family = x$Family[1],
             variable = x$name[1])
  }) %>% 
  bind_rows() 

abundance_diff_by_fam_and_host_var %>% 
  mutate(corrected_p = p.adjust(p.value, method = "BH")) %>% 
  View()

meta_dat %>% 
  skimr::skim(Antibiotics)



# Host associations: alpha diversity --------------------------------------


alpha_div <- 
  read_rds("data/diversity/alpha_div/alpha_diversity.Rds")

diversity_metrics_normality <- 
  read_tsv("data/diversity/alpha_div/diversity_metrics_normality.tsv", show_col_types = FALSE)

alpha_div_assoc <- 
  read_rds("data/diversity/alpha_div/alpha_associations.Rds")

load("data/diversity/beta_div/permanova_tests_lifestyle_demography.RData")

diff_abund_table <- read_tsv("data/diff_abund/maaslin2_vOTUs_amgs_vs_life_demo_clin.tsv", show_col_types = FALSE)


meta_dat %>% 
  left_join(sample_data %>% select(sample_id, deltaker_id, Total_Reads_QC_ATLAS)) %>% 
  # rename(var_1 = PhysAct_Score) %>% 
  # rename(var_1 = Alko) %>% 
  rename(var_1 = Karboh_energi) %>% 
  left_join(alpha_div %>% select(invsimpson, sample_id)) %>% 
  mutate(var_cat = cat_func(var_1),
         reads_cat = cat_func(Total_Reads_QC_ATLAS)) %>% 
  ggplot(aes(x = reads_cat, y = invsimpson, fill = var_cat)) +
  geom_boxplot()

## The vOTU vOTU05693 is worth taking a closer look at


vOTU_taxonomy %>% 
  group_by(VC) %>% 
  mutate(same_cluster = any(Scaffold %in% "vOTU05693")) %>% 
  ungroup() %>% 
  filter(same_cluster)
## Only other unknown in the same cluster


## But what about indirectly connected
gbg_vcontact <- read_csv("data/mid/vcontact/genome_by_genome_overview.csv")
gbg_vcontact %>% 
  filter()
vOTU_taxonomy %>% 
  group_by(VC) %>% 
  mutate(same_cluster = any(Scaffold %in% "vOTU05693")) %>% 
  ungroup() %>% 
  filter(same_cluster) %>% 
  left_join(gbg_vcontact %>% select(VC, preVC)) %>% 
  unique()
## So, the precluster is the same as the cluster  
  
vOTU_stats %>% 
  filter(vOTU %in% "vOTU05693")

