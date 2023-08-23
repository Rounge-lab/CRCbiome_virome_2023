

alpha_diversity_calculations <- function() {
  
  ## Calculate alpha diversity
  alpha_div <- 
    virus_abundance %>% 
    as.data.frame() %>% 
    column_to_rownames("sample_id") %>% 
    summarize(shannon = diversity(x = ., index = "shannon"),
              # simpson = diversity(x = ., index = "simpson"),
              invsimpson = diversity(x = ., index = "invsimpson"),
              observed = apply(., 1, function(x) sum(x>0))) %>% 
    tibble() %>% 
    bind_cols("sample_id" = virus_abundance %>% pull(sample_id)) 
  
  alpha_div %>% 
    write_rds("data/diversity/alpha_div/alpha_diversity.Rds")
  
  diversity_metrics_normality <- 
    alpha_div %>% 
    pivot_longer(c(shannon, invsimpson, observed), 
                 names_to = "index", values_to = "diversity") %>% 
    group_by(index) %>% 
    group_split() %>% 
    lapply(function(x) {
      x %>% 
        pull(diversity) %>% 
        shapiro.test() %>% 
        tidy() %>% 
        mutate(index = x$index[1])
    }) %>% 
    bind_rows()
  
  diversity_metrics_normality %>% 
    write_tsv("data/diversity/alpha_div/diversity_metrics_normality.tsv")
  
  ## Measure associations
  alpha_div_assoc <- 
    lapply(meta_variables$var_id, function(i) {
      
      alpha_div %>% 
        left_join(meta_dat %>% 
                    left_join(sample_data %>% select(deltaker_id, sample_id, Total_Reads_raw_ATLAS)) %>% 
                    select(sample_id, all_of(i), Total_Reads_raw_ATLAS), by = "sample_id") %>% 
        rename("var" = all_of(i)) %>% 
        filter(!is.na(var), !var %in% "Missing", !var %in% "Unknown") %>% 
        pivot_longer(cols = c(shannon, invsimpson, observed), 
                     names_to = "index", values_to = "diversity") %>% 
        filter(!is.na(var), !var %in% "Missing", !var %in% "Unknown", !var %in% "4. Other lesions") %>% 
        mutate(var = cat_func(var)) %>% 
        filter(!var %in% c("mid")) %>%
        group_by(index) %>% 
        group_split() %>% 
        lapply(function(x) {
          
          tmp <-
            x %>% 
            lm(diversity~var + Total_Reads_raw_ATLAS, data = .) %>%
            aov() 
          tmp %>% 
            tidy() %>% 
            left_join(tmp %>% omega_squared() %>% tibble() %>% rename(term = Parameter)) %>% 
            mutate(var_id = i,
                   index = x$index[1],
                   n = nrow(x))
        }) %>% 
        bind_rows()
    }) %>% 
    bind_rows() %>% 
    group_by(index) %>% 
    mutate(corrected_p = p.adjust(p = p.value, method = "BH")) %>% 
    ungroup() %>% 
    filter(!grepl("Intercept", term)) %>% 
    filter(!grepl("esidual", term)) %>% 
    left_join(meta_variables, by = "var_id") 
  
  alpha_div_assoc %>% 
    write_rds("data/diversity/alpha_div/alpha_associations.Rds")
  
}

beta_diversity_calculations <- function() {
  ## Percent variation attributed to variable - continuous as is
  effect_size_cont <- 
    lapply(meta_variables$var_id, function(i) {
      
      print(paste("variable ", which(meta_variables$var_id %in% i), " of ", length(meta_variables$var_id)))
      
      tmp_dat <- 
        virus_abundance %>% 
        left_join(meta_dat %>% left_join(sample_data %>% select(deltaker_id, sample_id)) %>% select(sample_id, all_of(i)), by = "sample_id") %>% 
        rename("var" = all_of(i)) %>% 
        filter(!is.na(var), !var %in% "Missing", !var %in% "Unknown")
      if (is.factor(tmp_dat$var)) {
        tmp_dat <-
          tmp_dat %>% 
          filter(!var %in% (table(var) %>% enframe() %>% filter(value < 10) %>% pull(name)))
      }
      
      tmp <- 
        adonis2(tmp_dat %>% select(-c(var, sample_id)) ~ tmp_dat %>% pull(var), method = "bray") %>% 
        adonis_OmegaSq()
      
      tibble("var" = i,
             "R2" = tmp$aov.tab$R2[1],
             "parOmegaSq" = tmp$aov.tab$parOmegaSq[1],
             "p" = tmp$aov.tab$`Pr(>F)`[1], 
             "df" = tmp$aov.tab$Df[1],
             "var_details" = tmp_dat %>% pull(var) %>% typeof(),
             "n" = tmp_dat %>% nrow())
    }) %>% 
    bind_rows()
  
  
  ## Percent variation attributed to variable - tertiles - high v low
  
  effect_size_cat <- 
    
    lapply(meta_variables$var_id, function(i) {
      
      print(paste("variable ", which(meta_variables$var_id %in% i), " of ", length(meta_variables$var_id)))
      
      tmp_dat <- 
        virus_abundance %>% 
        left_join(meta_dat %>% left_join(sample_data %>% select(deltaker_id, sample_id)) %>% select(sample_id, all_of(i)), by = "sample_id") %>% 
        rename("var" = all_of(i)) %>% 
        filter(!is.na(var), !var %in% "Missing", !var %in% "Unknown")
      if (is.factor(tmp_dat$var)) {
        tmp_dat <-
          tmp_dat %>% 
          filter(!var %in% (table(var) %>% enframe() %>% filter(value < 10) %>% pull(name)))
      }
      tmp_dat <- 
        tmp_dat %>% 
        mutate(var = cat_func(var)) %>% 
        filter(!var %in% c("mid"))
      
      tmp <- adonis(tmp_dat %>% select(-c(var, sample_id)) ~ tmp_dat %>% pull(var), method = "bray") %>% 
        adonis_OmegaSq()
      
      tibble("var" = i,
             "R2" = tmp$aov.tab$R2[1],
             "parOmegaSq" = tmp$aov.tab$parOmegaSq[1],
             "p" = tmp$aov.tab$`Pr(>F)`[1], 
             "df" = tmp$aov.tab$Df[1],
             "var_details" = tmp_dat %>% pull(var) %>% typeof(),
             "n" = tmp_dat %>% nrow())
    }) %>% 
    bind_rows()
  
  
  save(effect_size_cont, effect_size_cat, file = "data/diversity/beta_div/permanova_tests_lifestyle_demography.RData")
  
}



alpha_diversity_calculations_no_crc <- function() {
  
  no_crc_ids <- 
    read_rds("data/participant_data/screening_data.Rds") %>% 
    filter(!str_detect(final_result, "6. Cancer")) %>% 
    select(deltaker_id)
  
  alpha_div <- read_rds("data/diversity/alpha_div/alpha_diversity.Rds")
  
  ## Measure associations
  alpha_div_assoc_no_crc <- 
    lapply(meta_variables$var_id, function(i) {
      
      alpha_div %>% 
        left_join(meta_dat %>% 
                    inner_join(no_crc_ids, by = "deltaker_id") %>% 
                    left_join(sample_data %>% select(deltaker_id, sample_id, Total_Reads_raw_ATLAS), by = "deltaker_id") %>% 
                    select(sample_id, all_of(i), Total_Reads_raw_ATLAS), by = "sample_id") %>% 
        rename("var" = all_of(i)) %>% 
        filter(!is.na(var), !var %in% "Missing", !var %in% "Unknown") %>% 
        pivot_longer(cols = c(shannon, invsimpson, observed), 
                     names_to = "index", values_to = "diversity") %>% 
        filter(!is.na(var), !var %in% "Missing", !var %in% "Unknown", !var %in% "4. Other lesions") %>% 
        mutate(var = cat_func(var)) %>% 
        filter(!var %in% c("mid")) %>%
        group_by(index) %>% 
        group_split() %>% 
        lapply(function(x) {
          
          tmp <-
            x %>% 
            lm(diversity~var + Total_Reads_raw_ATLAS, data = .) %>%
            aov() 
          tmp %>% 
            tidy() %>% 
            left_join(tmp %>% omega_squared() %>% tibble() %>% rename(term = Parameter)) %>% 
            mutate(var_id = i,
                   index = x$index[1],
                   n = nrow(x))
        }) %>% 
        bind_rows()
    }) %>% 
    bind_rows() %>% 
    group_by(index) %>% 
    mutate(corrected_p = p.adjust(p = p.value, method = "BH")) %>% 
    ungroup() %>% 
    filter(!grepl("Intercept", term)) %>% 
    filter(!grepl("esidual", term)) %>% 
    left_join(meta_variables, by = "var_id") 
  
  alpha_div_assoc_no_crc %>% 
    write_rds("data/diversity/alpha_div/alpha_associations_no_crc.Rds")
  
}

beta_diversity_calculations_no_crc <- function() {
  
  no_crc_ids <- 
    read_rds("data/participant_data/screening_data.Rds") %>% 
    filter(!str_detect(final_result, "6. Cancer")) %>% 
    select(deltaker_id)
  
  ## Percent variation attributed to variable - continuous as is
  effect_size_cont_no_crc <- 
    lapply(meta_variables$var_id, function(i) {
      
      print(paste("variable ", which(meta_variables$var_id %in% i), " of ", length(meta_variables$var_id)))
      
      tmp_dat <- 
        virus_abundance %>% 
        left_join(meta_dat %>% 
                    inner_join(no_crc_ids, by = "deltaker_id") %>% 
                    left_join(sample_data %>% select(deltaker_id, sample_id)) %>% 
                    select(sample_id, all_of(i)), by = "sample_id") %>% 
        rename("var" = all_of(i)) %>% 
        filter(!is.na(var), !var %in% "Missing", !var %in% "Unknown")
      if (is.factor(tmp_dat$var)) {
        tmp_dat <-
          tmp_dat %>% 
          filter(!var %in% (table(var) %>% enframe() %>% filter(value < 10) %>% pull(name)))
      }
      
      tmp <- 
        adonis(tmp_dat %>% select(-c(var, sample_id)) ~ tmp_dat %>% pull(var), method = "bray") %>% 
        adonis_OmegaSq()
      tibble("var" = i,
             "R2" = tmp$aov.tab$R2[1],
             "parOmegaSq" = tmp$aov.tab$parOmegaSq[1],
             "p" = tmp$aov.tab$`Pr(>F)`[1], 
             "df" = tmp$aov.tab$Df[1],
             "var_details" = tmp_dat %>% pull(var) %>% typeof(),
             "n" = tmp_dat %>% nrow())
    }) %>% 
    bind_rows()
  
  
  ## Percent variation attributed to variable - tertiles - high v low
  
  effect_size_cat_no_crc <- 
    
    lapply(meta_variables$var_id, function(i) {
      
      print(paste("variable ", which(meta_variables$var_id %in% i), " of ", length(meta_variables$var_id)))
      
      tmp_dat <- 
        virus_abundance %>% 
        left_join(meta_dat %>% 
                    inner_join(no_crc_ids, by = "deltaker_id") %>%
                    left_join(sample_data %>% select(deltaker_id, sample_id)) %>% 
                    select(sample_id, all_of(i)), by = "sample_id") %>% 
        rename("var" = all_of(i)) %>% 
        filter(!is.na(var), !var %in% "Missing", !var %in% "Unknown")
      if (is.factor(tmp_dat$var)) {
        tmp_dat <-
          tmp_dat %>% 
          filter(!var %in% (table(var) %>% enframe() %>% filter(value < 10) %>% pull(name)))
      }
      tmp_dat <- 
        tmp_dat %>% 
        mutate(var = cat_func(var)) %>% 
        filter(!var %in% c("mid"))
      
      tmp <- adonis(tmp_dat %>% select(-c(var, sample_id)) ~ tmp_dat %>% pull(var), method = "bray") %>% 
        adonis_OmegaSq()
      
      tibble("var" = i,
             "R2" = tmp$aov.tab$R2[1],
             "parOmegaSq" = tmp$aov.tab$parOmegaSq[1],
             "p" = tmp$aov.tab$`Pr(>F)`[1], 
             "df" = tmp$aov.tab$Df[1],
             "var_details" = tmp_dat %>% pull(var) %>% typeof(),
             "n" = tmp_dat %>% nrow())
    }) %>% 
    bind_rows()
  
  
  save(effect_size_cont_no_crc, effect_size_cat_no_crc, file = "data/diversity/beta_div/permanova_tests_lifestyle_demography_no_crc.RData")
  
}