


process_amg_data <- function() {
  
  # load data ---------------------------------------------------------------
  
  vOTU_amg_summary <- read_tsv("data/mid/dramv/distillate/vMAG_stats.tsv") %>% 
    dplyr::rename(scaffold = 1)
  amg_summary <- read_tsv("data/mid/dramv/distillate/amg_summary.tsv")
  
  if (!any(ls() %in% "virus_abundance")) {
    virus_abundance <-
      read_tsv("data/mid/abundance/coverage_table_min75.tsv")  
  }
  
  
  # summarize amg table -----------------------------------------------------
  
  ## Summarize by presence - if a there is at least one amg with some annotation, this is considered "presence"
  amg_cat_by_vOTU <-
    vOTU_amg_summary %>% 
    select(scaffold) %>% 
    left_join(amg_summary %>% 
                count(scaffold, category) %>% 
                mutate(n = n > 0) %>% 
                pivot_wider(names_from = category, values_from = n, values_fill = FALSE) %>% 
                select(-`NA`)) %>% 
    mutate(vOTU = str_extract(scaffold, "^vOTU[:digit:]*")) %>% 
    select(vOTU, everything(), -scaffold) %>% 
    mutate(across(.cols = -vOTU, .fns = function(x) ifelse(x %in% TRUE, TRUE, FALSE)))
  
  amg_subhead_by_vOTU <-
    vOTU_amg_summary %>% 
    select(scaffold) %>% 
    left_join(amg_summary %>% 
                count(scaffold, subheader) %>% 
                mutate(n = n > 0) %>% 
                pivot_wider(names_from = subheader, values_from = n, values_fill = FALSE) %>% 
                select(-`NA`)) %>% 
    mutate(vOTU = str_extract(scaffold, "^vOTU[:digit:]*")) %>% 
    select(vOTU, everything(), -scaffold) %>% 
    mutate(across(.cols = -vOTU, .fns = function(x) ifelse(x %in% TRUE, TRUE, FALSE)))
  
  amg_module_by_vOTU <-
    vOTU_amg_summary %>% 
    select(scaffold) %>% 
    left_join(amg_summary %>% 
                count(scaffold, module) %>% 
                mutate(n = n > 0) %>% 
                pivot_wider(names_from = module, values_from = n, values_fill = FALSE) %>% 
                select(-`NA`)) %>% 
    mutate(vOTU = str_extract(scaffold, "^vOTU[:digit:]*")) %>% 
    select(vOTU, everything(), -scaffold) %>% 
    mutate(across(.cols = -vOTU, .fns = function(x) ifelse(x %in% TRUE, TRUE, FALSE)))
  
  amg_MISC_module_by_vOTU <-
    vOTU_amg_summary %>% 
    select(scaffold) %>% 
    left_join(amg_summary %>% 
                filter(category %in% "MISC") %>% 
                count(scaffold, module) %>% 
                mutate(n = n > 0) %>% 
                pivot_wider(names_from = module, values_from = n, values_fill = FALSE)) %>% 
    mutate(vOTU = str_extract(scaffold, "^vOTU[:digit:]*")) %>% 
    select(vOTU, everything(), -scaffold) %>% 
    mutate(across(.cols = -vOTU, .fns = function(x) ifelse(x %in% TRUE, TRUE, FALSE)))
  
  # combine with presence per sample ----------------------------------------
  
  amg_cat_by_sample <-
    virus_abundance %>% 
    pivot_longer(-ID, names_to = "sample_id", values_to = "abundance") %>% 
    dplyr::rename(vOTU = ID) %>% 
    mutate(vOTU_presence = abundance > 0) %>% 
    filter(vOTU_presence) %>% 
    select(sample_id, vOTU, vOTU_presence) %>% 
    left_join(amg_cat_by_vOTU %>% 
                pivot_longer(-vOTU, names_to = "category", values_to = "amg_presence")) %>% 
    group_by(sample_id, category) %>% 
    summarize(n_otu_amgs = sum(amg_presence, na.rm = TRUE)) %>% 
    filter(!is.na(category)) %>% 
    pivot_wider(names_from = category, values_from = n_otu_amgs, values_fill = 0) 
  
  amg_subhead_by_sample <-
    virus_abundance %>% 
    pivot_longer(-ID, names_to = "sample_id", values_to = "abundance") %>% 
    dplyr::rename(vOTU = ID) %>% 
    mutate(vOTU_presence = abundance > 0) %>% 
    filter(vOTU_presence) %>% 
    select(sample_id, vOTU, vOTU_presence) %>% 
    left_join(amg_subhead_by_vOTU %>% 
                pivot_longer(-vOTU, names_to = "subhead", values_to = "amg_presence")) %>% 
    group_by(sample_id, subhead) %>% 
    summarize(n_otu_amgs = sum(amg_presence, na.rm = TRUE)) %>% 
    filter(!is.na(subhead)) %>% 
    pivot_wider(names_from = subhead, values_from = n_otu_amgs, values_fill = 0) 
  
  amg_module_by_sample <-
    virus_abundance %>% 
    pivot_longer(-ID, names_to = "sample_id", values_to = "abundance") %>% 
    dplyr::rename(vOTU = ID) %>% 
    mutate(vOTU_presence = abundance > 0) %>% 
    filter(vOTU_presence) %>% 
    select(sample_id, vOTU, vOTU_presence) %>% 
    left_join(amg_module_by_vOTU %>% 
                pivot_longer(-vOTU, names_to = "module", values_to = "amg_presence")) %>% 
    group_by(sample_id, module) %>% 
    summarize(n_otu_amgs = sum(amg_presence, na.rm = TRUE)) %>% 
    filter(!is.na(module)) %>% 
    pivot_wider(names_from = module, values_from = n_otu_amgs, values_fill = 0) 
  
  amg_cat_abundance <-
    virus_abundance %>% 
    pivot_longer(-ID, names_to = "sample_id", values_to = "abundance") %>% 
    dplyr::rename(vOTU = ID) %>% 
    filter(abundance > 0) %>% 
    group_by(sample_id) %>% 
    group_split() %>% 
    lapply(function(x) {
      x %>% 
        inner_join(amg_cat_by_vOTU %>% 
                     pivot_longer(-vOTU, names_to = "cat", values_to = "amg_presence") %>% 
                     filter(amg_presence) %>% 
                     select(-amg_presence), by = "vOTU")
    }) %>% 
    bind_rows() %>% 
    group_by(cat, sample_id) %>% 
    summarize(abundance = sum(abundance, na.rm = TRUE)) %>% 
    ungroup() %>% 
    pivot_wider(names_from = cat, values_from = abundance, values_fill = 0)
  
  amg_subhead_abundance <-
    virus_abundance %>% 
    pivot_longer(-ID, names_to = "sample_id", values_to = "abundance") %>% 
    dplyr::rename(vOTU = ID) %>% 
    filter(abundance > 0) %>% 
    group_by(sample_id) %>% 
    group_split() %>% 
    lapply(function(x) {
      x %>% 
        inner_join(amg_subhead_by_vOTU %>% 
                     pivot_longer(-vOTU, names_to = "subhead", values_to = "amg_presence") %>% 
                     filter(amg_presence) %>% 
                     select(-amg_presence), by = "vOTU")
    }) %>% 
    bind_rows() %>% 
    group_by(subhead, sample_id) %>% 
    summarize(abundance = sum(abundance, na.rm = TRUE)) %>% 
    ungroup() %>% 
    pivot_wider(names_from = subhead, values_from = abundance, values_fill = 0)
  
  amg_module_abundance <-
    virus_abundance %>% 
    pivot_longer(-ID, names_to = "sample_id", values_to = "abundance") %>% 
    dplyr::rename(vOTU = ID) %>% 
    filter(abundance > 0) %>% 
    group_by(sample_id) %>% 
    group_split() %>% 
    lapply(function(x) {
      x %>% 
        inner_join(amg_module_by_vOTU %>% 
                     pivot_longer(-vOTU, names_to = "module", values_to = "amg_presence") %>% 
                     filter(amg_presence) %>% 
                     select(-amg_presence), by = "vOTU")
    }) %>% 
    bind_rows() %>% 
    group_by(module, sample_id) %>% 
    summarize(abundance = sum(abundance, na.rm = TRUE)) %>% 
    ungroup() %>% 
    pivot_wider(names_from = module, values_from = abundance, values_fill = 0)
    
  # write to file -----------------------------------------------------------
  
  amg_cat_by_vOTU %>% 
    write_tsv("data/amgs/amg_cat_by_vOTU_mid_min75.tsv")
  amg_module_by_vOTU %>% 
    write_tsv("data/amgs/amg_module_by_vOTU_mid_min75.tsv")
  amg_subhead_by_vOTU %>% 
    write_tsv("data/amgs/amg_subhead_by_vOTU_mid_min75.tsv")
  amg_MISC_module_by_vOTU %>% 
    write_tsv("data/amgs/amg_MISC_module_by_vOTU_mid_min75.tsv")
  
  amg_cat_by_sample %>% 
    write_tsv("data/amgs/amg_cat_by_sample_mid_min75.tsv")
  amg_module_by_sample %>% 
    write_tsv("data/amgs/amg_module_by_sample_mid_min75.tsv")
  amg_subhead_by_sample %>% 
    write_tsv("data/amgs/amg_subhead_by_sample_mid_min75.tsv")
  
  amg_cat_abundance %>% 
    write_tsv("data/amgs/amg_cat_abundance_mid_min75.tsv")
  amg_subhead_abundance %>% 
    write_tsv("data/amgs/amg_subhead_abundance_mid_min75.tsv")
  amg_module_abundance %>% 
    write_tsv("data/amgs/amg_module_abundance_mid_min75.tsv")
}
