

create_vOTU_based_stats <- function() {
  
  ## Summarize taxonomy, lifecycle, abundance, prevalence
  vOTU_by_taxonomy_lifecycle <- 
    create_virus_taxonomy_lifecycle_table(vir_abund = virus_abundance, 
                                          vir_stats = vOTU_stats, 
                                          vir_lifecycle = vOTU_lifecycle_cat, 
                                          min_n_for_summary = 20) 
  
  vOTU_by_taxonomy_lifecycle %>% 
    write_rds("tables/vOTU_based_stats/vOTU_by_taxonomy_lifecycle.Rds")
  
  
  lifecycle_by_family <-
    vOTU_stats %>% 
    select(vOTU, cluster_not_proviruses, cluster_proviruses) %>% 
    left_join(vOTU_by_taxonomy_lifecycle %>% select(vOTU, Family)) %>% 
    filter(!is.na(Family)) %>% 
    group_by(Family) %>% 
    summarize(n_lyso = sum(cluster_proviruses),
              n_lytic = sum(cluster_not_proviruses),
              n_total = n_lyso + n_lytic,
              fraction_lyso = n_lyso/n_total) %>% 
    ungroup() 
  
  lifecycle_by_family %>% 
    write_rds("tables/vOTU_based_stats/lifecycle_by_family.Rds")
  
  ## Test for differences in proportions in family/lifecycle
  lapply(levels(lifecycle_by_family$Family), function(x) {
    lifecycle_by_family %>%
      select(Family, n_lyso, n_lytic) %>% 
      pivot_longer(starts_with("n_"), names_to = "lifecycle") %>% 
      mutate(tmp_fam = case_when(Family %in% x ~ "interesting_family",
                                 TRUE  ~ "others")) %>% 
      select(-Family) %>% 
      group_by(lifecycle, tmp_fam) %>% 
      summarize(n = sum(value)) %>% 
      pivot_wider(names_from = tmp_fam, values_from = n, values_fill = 0) %>% 
      as.data.frame() %>% 
      column_to_rownames("lifecycle") %>% 
      fisher.test() %>% 
      tidy() %>% 
      mutate(Family = x) %>% 
      left_join(lifecycle_by_family)
  }) %>% 
    bind_rows() %>% 
    mutate(Family = factor(Family, levels = rev(levels(lifecycle_by_family$Family)))) %>% 
    mutate(p_corr = p.adjust(p.value, method = "bonferroni")) %>% 
    write_rds("tables/vOTU_based_stats/family_lifecycle_tests.Rds")
  
  
  amgs_per_fam_frac <-
    amgs_cat_per_vOTU %>% 
    rename(`Carbon Utilization` = `carbon utilization`) %>%
    pivot_longer(-vOTU, names_to = "cat", values_to = "presence") %>% 
    filter(!str_detect(cat, "Wood")) %>% 
    left_join(vOTU_by_taxonomy_lifecycle %>% select(vOTU, cluster_contigs, Family)) %>% 
    filter(!is.na(Family)) %>% 
    group_by(Family, cat) %>% 
    summarize(with_cat = sum(presence)/n(),
              sum_with = sum(presence),
              sum_without = sum(!presence)) %>% 
    ungroup()
  
  ## Test for differences in proportions in family/lifecycle
  lapply(levels(amgs_per_fam_frac$Family), function(fam) {
    lapply(unique(amgs_per_fam_frac$cat), function(amg_cat) {
      amgs_per_fam_frac %>%
        filter(cat %in% amg_cat) %>% 
        select(Family, sum_with, sum_without) %>% 
        pivot_longer(starts_with("sum_"), names_to = "presence") %>% 
        mutate(tmp_fam = case_when(Family %in% fam ~ "interesting_family",
                                   TRUE  ~ "others")) %>% 
        select(-Family) %>% 
        group_by(presence, tmp_fam) %>% 
        summarize(n = sum(value)) %>% 
        pivot_wider(names_from = tmp_fam, values_from = n, values_fill = 0) %>% 
        as.data.frame() %>% 
        column_to_rownames("presence") %>% 
        fisher.test() %>% 
        tidy() %>% 
        mutate(Family = fam,
               cat = amg_cat) %>% 
        left_join(amgs_per_fam_frac)  
    }) %>% 
      bind_rows()
  }) %>% 
    bind_rows() %>% 
    mutate(Family = factor(Family, levels = rev(levels(amgs_per_fam_frac$Family)))) %>% 
    mutate(p_corr = p.adjust(p.value, method = "bonferroni")) %>% 
    write_rds("tables/vOTU_based_stats/family_AMG_tests.Rds")
  
}

