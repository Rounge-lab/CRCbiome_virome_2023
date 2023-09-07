

# env ---------------------------------------------------------------------

library(Maaslin2)

# run maaslin -------------------------------------------------------------

diff_abund_results <-
  lapply(c("vOTUs"), function(dataset) {
    if (dataset == "vOTUs") tmp_abund <- virus_abundance
    
    tmp_abund <- tmp_abund %>% mutate(sample_id = str_replace(sample_id, "-", "_"))
    
    if (!dataset %in% list.dirs("data/diff_abund/", full.names = FALSE)) dir.create(paste0("data/diff_abund/", dataset))
    
    meta_dat %>% 
      left_join(sample_meta %>% select(sample_id, deltaker_id)) %>% 
      select(-deltaker_id) %>% 
      mutate(across(where(is.double), .fns = function(x) cat_func(x) %>% str_replace("negative", "a_low") %>% str_replace("positive", "b_high"))) %>% 
      pivot_longer(-sample_id) %>% 
      # filter(name %in% "kjonn") -> x
      group_by(name) %>% 
      group_split() %>% 
      lapply(function(x) {
        tmp_ref <- meta_dat_vars %>% filter(var_id %in% x$name[1]) %>% pull(ref)
        if (is.na(tmp_ref)) {
          x <- x %>% 
            filter(value != "mid")
          tmp_ref <- "a_low"
        }
        x <- 
          x %>% 
          filter(!value %in% c("Missing", "Unknown"))
        
        if (!x$name[1] %in% c("kjonn", "age_cat", "senter")) {
          x <- 
            x %>% 
            left_join(sample_meta %>% left_join(meta_dat, by = "deltaker_id") %>% select(sample_id, kjonn, age_cat, senter), by = "sample_id")
          mod <- paste0("value,kjonn,age_cat,senter")
          refs <- paste0("value,", tmp_ref,
                         ";kjonn,", levels(meta_dat$kjonn)[1], 
                         ";age_cat,", levels(meta_dat$age_cat)[1],
                         ";senter,", levels(meta_dat$senter)[1])
        } else {
          if (x$name[1] %in% "kjonn") {
            x <- 
              x %>% 
              left_join(sample_meta %>% left_join(meta_dat, by = "deltaker_id") %>% select(sample_id, age_cat, senter), by = "sample_id")
            mod <- paste0("value,age_cat,senter")
            refs <- paste0("value,", tmp_ref, 
                          ";age_cat,", levels(meta_dat$age_cat)[1],
                          ";senter,", levels(meta_dat$senter)[1])
            
          }
          if (x$name[1] %in% "age_cat") {
            x <- 
              x %>% 
              left_join(sample_meta %>% left_join(meta_dat, by = "deltaker_id") %>% select(sample_id, kjonn, senter), by = "sample_id")
            mod <- paste0("value,kjonn,senter")
            refs <- paste0("value,", tmp_ref, 
                          ";kjonn,", levels(meta_dat$kjonn)[1],
                          ";senter,", levels(meta_dat$senter)[1])
            
          }
          if (x$name[1] %in% "senter") {
            x <- 
              x %>% 
              left_join(sample_meta %>% left_join(meta_dat, by = "deltaker_id") %>% select(sample_id, kjonn, age_cat), by = "sample_id")
            mod <- paste0("value,kjonn,age_cat")
            refs <- paste0("value,", tmp_ref, 
                          ";kjonn,", levels(meta_dat$kjonn)[1],
                          ";age_cat,", levels(meta_dat$age_cat)[1])
            
          }
          
        }
        
        tmp_maaslin <-
          Maaslin2(
            input_data = x %>% 
              select(sample_id) %>% 
              left_join(tmp_abund) %>% 
              as.data.frame() %>% 
              column_to_rownames("sample_id"),
            input_metadata = x %>% 
              column_to_rownames("sample_id"),
            min_prevalence = ifelse(any(grepl("mock", ls())), 0.001, 0.1),
            normalization = ifelse(dataset %in% "vOTUs", "TSS", "NONE"),
            analysis_method = "LM",
            output = paste0("data/diff_abund/", dataset, "/", x$name[1]),
            ## Add adjustment: Age + sex
            fixed_effects = mod,
            reference = refs,
            plot_heatmap = FALSE,
            plot_scatter = FALSE
          )
        
        tmp_maaslin$results %>% 
          tibble() %>% 
          mutate(variable = x$name[1],
                 dataset = dataset) %>% 
          filter(metadata %in% "value")
        
      }) %>% 
      bind_rows()
  }) %>% 
  bind_rows()

diff_abund_results %>% 
  write_tsv("data/diff_abund/maaslin2_vOTUs_amgs_vs_life_demo_clin.tsv")



# sensitivity analysis ----------------------------------------------------

no_crc_ids <- 
  read_rds("data/participant_data/screening_data.Rds") %>% 
  filter(!str_detect(final_result, "^6. Cancer")) %>% 
  select(deltaker_id)

diff_abund_results <-
  lapply(c("vOTUs"), function(dataset) {
    if (dataset == "vOTUs") tmp_abund <- virus_abundance
    
    tmp_abund <- tmp_abund %>% mutate(sample_id = str_replace(sample_id, "-", "_"))
    
    if (!dataset %in% list.dirs("data/diff_abund/no_crc/", full.names = FALSE)) dir.create(paste0("data/diff_abund/no_crc/", dataset))
    
    meta_dat %>% 
      inner_join(no_crc_ids) %>% 
      left_join(sample_meta %>% select(sample_id, deltaker_id)) %>% 
      select(-deltaker_id) %>% 
      mutate(across(where(is.double), .fns = function(x) cat_func(x) %>% str_replace("negative", "a_low") %>% str_replace("positive", "b_high"))) %>% 
      pivot_longer(-sample_id) %>% 
      # filter(name %in% "kjonn") -> x
      group_by(name) %>% 
      group_split() %>% 
      lapply(function(x) {
        tmp_ref <- meta_dat_vars %>% filter(var_id %in% x$name[1]) %>% pull(ref)
        if (is.na(tmp_ref)) {
          x <- x %>% 
            filter(value != "mid")
          tmp_ref <- "a_low"
        }
        x <- 
          x %>% 
          filter(!value %in% c("Missing", "Unknown"))
        
        if (!x$name[1] %in% c("kjonn", "age_cat", "senter")) {
          x <- 
            x %>% 
            left_join(sample_meta %>% left_join(meta_dat, by = "deltaker_id") %>% select(sample_id, kjonn, age_cat, senter), by = "sample_id")
          mod <- paste0("value,kjonn,age_cat,senter")
          refs <- paste0("value,", tmp_ref,
                         ";kjonn,", levels(meta_dat$kjonn)[1], 
                         ";age_cat,", levels(meta_dat$age_cat)[1],
                         ";senter,", levels(meta_dat$senter)[1])
        } else {
          if (x$name[1] %in% "kjonn") {
            x <- 
              x %>% 
              left_join(sample_meta %>% left_join(meta_dat, by = "deltaker_id") %>% select(sample_id, age_cat, senter), by = "sample_id")
            mod <- paste0("value,age_cat,senter")
            refs <- paste0("value,", tmp_ref, 
                          ";age_cat,", levels(meta_dat$age_cat)[1],
                          ";senter,", levels(meta_dat$senter)[1])
            
          }
          if (x$name[1] %in% "age_cat") {
            x <- 
              x %>% 
              left_join(sample_meta %>% left_join(meta_dat, by = "deltaker_id") %>% select(sample_id, kjonn, senter), by = "sample_id")
            mod <- paste0("value,kjonn,senter")
            refs <- paste0("value,", tmp_ref, 
                          ";kjonn,", levels(meta_dat$kjonn)[1],
                          ";senter,", levels(meta_dat$senter)[1])
            
          }
          if (x$name[1] %in% "senter") {
            x <- 
              x %>% 
              left_join(sample_meta %>% left_join(meta_dat, by = "deltaker_id") %>% select(sample_id, kjonn, age_cat), by = "sample_id")
            mod <- paste0("value,kjonn,age_cat")
            refs <- paste0("value,", tmp_ref, 
                          ";kjonn,", levels(meta_dat$kjonn)[1],
                          ";age_cat,", levels(meta_dat$age_cat)[1])
            
          }
          
        }
        
        tmp_maaslin <-
          Maaslin2(
            input_data = x %>% 
              select(sample_id) %>% 
              left_join(tmp_abund) %>% 
              as.data.frame() %>% 
              column_to_rownames("sample_id"),
            input_metadata = x %>% 
              column_to_rownames("sample_id"),
            min_prevalence = ifelse(any(grepl("mock", ls())), 0.001, 0.1),
            normalization = ifelse(dataset %in% "vOTUs", "TSS", "NONE"),
            analysis_method = "LM",
            output = paste0("data/diff_abund/no_crc/", dataset, "/", x$name[1]),
            ## Add adjustment: Age + sex
            fixed_effects = mod,
            reference = refs,
            plot_heatmap = FALSE,
            plot_scatter = FALSE
          )
        
        tmp_maaslin$results %>% 
          tibble() %>% 
          mutate(variable = x$name[1],
                 dataset = dataset) %>% 
          filter(metadata %in% "value")
        
      }) %>% 
      bind_rows()
  }) %>% 
  bind_rows()

diff_abund_results %>% 
  write_tsv("data/diff_abund/no_crc/maaslin2_vOTUs_amgs_vs_life_demo_clin.tsv")


