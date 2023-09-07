

create_host_associations_figures <- function() {
  
  # load some data ----------------------------------------------------------
  alpha_div <- 
    read_rds("data/diversity/alpha_div/alpha_diversity.Rds")
  
  alpha_div_assoc <- 
    read_rds("data/diversity/alpha_div/alpha_associations.Rds")
  
  load("data/diversity/beta_div/permanova_tests_lifestyle_demography.RData")
  
  diff_abund_table <- read_tsv("data/diff_abund/maaslin2_vOTUs_amgs_vs_life_demo_clin.tsv", show_col_types = FALSE)
  
  
  ## Alpha plot
  invsimpson_assoc_cat <-
    alpha_div_assoc %>% 
    select(-c(var_name, dataset, ref)) %>% 
    left_join(meta_variables) %>% 
    filter(!grepl("Reads", term),
           index %in% "invsimpson") %>%
    rename(p = corrected_p) %>% 
    mutate(dataset = case_when(dataset %in% "diet and lifestyle" ~ "lifestyle",
                               TRUE ~ dataset)) %>% 
    mutate(dataset = "something") %>% 
    plot_pve(effect_size_var = "Omega2_partial", sumSq = TRUE, colors = four_cat_col[2])
  
  
  ## Beta plot
  effect_size_beta_cat <-
    effect_size_cat %>% 
    left_join(meta_variables %>% rename(var = var_id)) %>% 
    mutate(dataset = case_when(dataset %in% "diet and lifestyle" ~ "lifestyle",
                               TRUE ~ dataset)) %>% 
    # filter(var_name != "WCRF") %>% 
    mutate(dataset = "something") %>% 
    plot_pve(colors = four_cat_col[4])
  
  ## Diff abund plot
  create_diff_abund_plot_table <- function() {
    n_diff_abund <- 
      diff_abund_table %>% 
      filter(dataset %in% "vOTUs") %>% 
      mutate(direction = case_when(qval < 0.05 & coef < 0 ~ "negative",
                                   qval < 0.05 & coef > 0 ~ "positive",
                                   TRUE ~ "no_assoc") %>% factor(levels = c("positive", "no_assoc", "negative"))) %>% 
      group_by(variable) %>% 
      count(direction) %>% 
      ungroup() %>% 
      mutate(number = ifelse(!direction %in% "no_assoc", n, 0)) %>% 
      mutate(n_for_plot = case_when(direction %in% "neg" ~ -number,
                                    TRUE ~ number)) %>% 
      left_join(meta_variables %>% rename(variable = var_id)) %>% 
      left_join(diff_abund_table %>% filter(dataset %in% "vOTUs") %>% select(variable, value, N) %>% unique() %>% group_by(variable, N) %>% summarize(n_groups = n()+1)) %>% 
      mutate(var_label = paste(var_name, " (g=", n_groups, ", n=", N, ")", sep = ""))
    n_diff_abund %>% 
      mutate(var_label = factor(var_label, levels = (n_diff_abund %>% group_by(var_label) %>% summarize(sign_n = sum(number)) %>% ungroup() %>% arrange(sign_n) %>% pull(var_label))))
  }
  
  n_diff_abund <- create_diff_abund_plot_table()
  
  diff_abund_plot <-
    n_diff_abund %>% 
    ggplot(aes(x = number, y = var_label, fill = direction)) +
    theme_bw()+
    geom_col() +
    guides(fill = guide_legend()) +
    scale_fill_manual(values = c("negative" = three_cat_col[1], "positive" = three_cat_col[3]), drop = TRUE)+
    labs(x = "differentially abundant vOTUs", y = "") +
    theme(legend.position = c(0.7, 0.2))
  
  
  tmp <- 
    diff_abund_table %>% 
    rename(vOTU = feature) %>% 
    left_join(vOTU_stats) %>% 
    left_join(meta_variables %>% rename(variable = var_id, var_dat = dataset) %>% select(variable, var_name, var_dat)) %>% 
    filter(dataset %in% "vOTUs",
           variable %in% c("Smoking", "PhysAct_Score", "Fiber")) %>% 
    mutate(var_name = factor(var_name, levels = c("Smoking", "Physical activity", "Fiber (g/day)"))) %>% 
    # filter(!variable %in% c("Mettet_energi", "C_enum_energi", "Trans_u_energi")) %>% 
    mutate(y_ax = -log10(qval))
  
  volcano_plot <-
    tmp %>% 
    ggplot(aes(x = coef, y = y_ax)) +
    # geom_point(color = tibble(v = c("demography", "diet", "lifestyle"),
    #                           col_ = c("#885DB2", "#B2885D", "#5DB288")) %>% filter(v %in% x) %>% pull(col_)) +
    geom_point(size = 0.5) +
    facet_wrap(~var_name, ncol = 1) +
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = 2) +
    theme_bw() +
    labs(y = "-log10(qval)",
         x = "log2FC") +
    scale_y_continuous(limits = c(0, 4.5))
  
  
  
  tmp %>% 
    ggplot(aes(x = coef, y = y_ax)) +
    # geom_point(color = tibble(v = c("demography", "diet", "lifestyle"),
    #                           col_ = c("#885DB2", "#B2885D", "#5DB288")) %>% filter(v %in% x) %>% pull(col_)) +
    geom_point(size = 0.5, data = tmp %>% filter(qval > 0.05)) +
    geom_point(data = tmp %>% filter(qval < 0.05), aes(color = Family, size = 2-fraction_annotated)) +
    ggrepel::geom_label_repel(data = tmp %>% filter(qval < 0.05) %>% filter(vOTU %in% c("vOTU05693", "vOTU12609")), aes(label = vOTU), nudge_y = 1.5) +
    facet_wrap(~var_name, ncol = 1) +
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = 2) +
    theme_bw() +
    labs(y = "-log10(qval)",
         x = "log2FC") +
    scale_y_continuous(limits = c(0, 5.5))
  
  
  diversity_composition_diff_abund_plot <- 
      grid.arrange(invsimpson_assoc_cat + labs(subtitle = "a)",
                                               x = expression("effect size ("*omega*"²)")) + 
                     # theme(legend.position = c(0.7, 0.2),
                     theme(legend.position = "none",
                           legend.background = element_rect(fill = "white", color = "black"),
                           legend.title = element_blank(),
                           text = element_text(size = 5))+
                     scale_x_continuous(limits = c(NA,0.02)),
                   effect_size_beta_cat + 
                     labs(subtitle = "b)",
                          x = expression("effect size ("*omega*"²)")) + 
                     theme(legend.position = "none",
                     # theme(legend.position = c(0.7, 0.2),
                           legend.background = element_rect(fill = "white", color = "black"),
                           legend.title = element_blank(),
                     text = element_text(size = 5)) +
                     scale_x_continuous(limits = c(NA,0.004)),
                   diff_abund_plot + 
                     labs(subtitle = "c)") +
                     theme(legend.position = c(0.5, 0.2),
                           legend.background = element_rect(fill = "white", color = "black"),
                           text = element_text(size = 5)),
                   volcano_plot + 
                     labs(subtitle = "d)") +
                     theme(text = element_text(size = 5)),
                   layout_matrix = rbind(c(1,1,2,2,3,3,4), c(1,1,2,2,3,3,4)))
  
  ggsave("figures/diversity_composition/diversity_composition_plot.png", plot = diversity_composition_diff_abund_plot, height = 180/15*7, width = 180, units = "mm")
  ggsave("figures/diversity_composition/diversity_composition_plot.pdf", plot = diversity_composition_diff_abund_plot, height = 180/15*7, width = 180, units = "mm")
  
  
  lab_func <- function(x) {
    gsub("value", "", gsub("b_h", "H", x))
  }
  
  # volcano_plot <-
  full_volcano <- 
    diff_abund_table %>% 
    left_join(meta_variables %>% rename(variable = var_id, var_dat = dataset) %>% select(variable, var_name, var_dat)) %>% 
    filter(dataset %in% "vOTUs") %>% 
    group_by(var_name) %>% 
    mutate(sum_sign = sum(qval < 0.05)) %>% 
    ungroup() %>% 
    mutate(var_name = fct_reorder(.f = var_name, .x = sum_sign, .fun = max, .desc = TRUE)) %>% 
    # filter(!variable %in% c("Mettet_energi", "C_enum_energi", "Trans_u_energi")) %>% 
    mutate(y_ax = -log10(qval)) %>% 
    ggplot(aes(x = coef, y = y_ax)) +
    # geom_point(color = tibble(v = c("demography", "diet", "lifestyle"),
    #                           col_ = c("#885DB2", "#B2885D", "#5DB288")) %>% filter(v %in% x) %>% pull(col_)) +
    geom_point() +
    facet_wrap(~var_name+name, labeller = as_labeller(lab_func), ncol = 4) +
    # facet_wrap(~var_name+name) +
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = 2) +
    theme_bw() +
    labs(y = "-log10(qval)",
         x = "log2FC") +
    scale_y_continuous(limits = c(0, 4.5))
  
  full_volcano %>% 
    ggsave(filename = "figures/diversity_composition/differential_abundance_full_vOTU_all_vars.png", plot = ., height = 12, width = 8)
  full_volcano %>% 
    ggsave(filename = "figures/diversity_composition/differential_abundance_full_vOTU_all_vars.pdf", plot = ., height = 12, width = 8)
  
  
  ## Assess sensitivity analyses: No CRC
  alpha_div_assoc_sens_comp <- 
    alpha_div_assoc %>% 
    filter(term %in% "var") %>% 
    mutate(sens = "full") %>% 
    left_join(read_rds("data/diversity/alpha_div/alpha_associations_no_crc.Rds") %>% 
                filter(term %in% "var") %>% 
                mutate(sens = "no_crc") %>% 
                rename(no_crc_es = Omega2_partial) %>% 
                select(no_crc_es, no_crc_corr_p = corrected_p, var_id, index))
    
  alpha_sensitivity_comp <- 
    alpha_div_assoc_sens_comp %>% 
    mutate(significance_codes = case_when(corrected_p < 0.05 & no_crc_corr_p < 0.05 ~ "both significant",
                                          !corrected_p < 0.05 & no_crc_corr_p < 0.05 ~ "no crc significant",
                                          corrected_p < 0.05 & !no_crc_corr_p < 0.05 ~ "original significant",
                                          TRUE ~ "none significant") %>% 
             factor(levels = c("none significant", "original significant", "no crc significant", "both significant"))) %>%
    ggplot(aes(x = Omega2_partial, y = no_crc_es, color = significance_codes)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_point(size = 3) +
    theme_bw() +
    labs(x = expression("Original effect size ("*omega*"²)"), 
         y = expression("Effect size ("*omega*"²) without CRC"),
         # shape = "Original FDR < 0.05",
         color = "FDR significance") +
    scale_color_manual(values = c("gray", four_cat_col[c(2,1)]))
  
  beta_sensitivity_comp <-
    effect_size_cat %>% 
    left_join(effect_size_cat_no_crc %>% select(var, no_crc_es = parOmegaSq, no_crc_p = p)) %>% 
    mutate(significance_codes = case_when(p < 0.05 & no_crc_p < 0.05 ~ "both significant",
                                          !p < 0.05 & no_crc_p < 0.05 ~ "no crc significant",
                                          p < 0.05 & !no_crc_p < 0.05 ~ "original significant",
                                          TRUE ~ "none significant") %>% 
             factor(levels = c("none significant", "original significant", "no crc significant", "both significant"))) %>% 
    ggplot(aes(x = parOmegaSq, y = no_crc_es, color = significance_codes)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_point(size = 3) +
    theme_bw() +
    labs(x = expression("Original effect size ("*omega*"²)"), 
         y = expression("Effect size ("*omega*"²) without CRC"),
         # shape = "Original FDR < 0.05",
         color = "FDR significance") +
    scale_color_manual(values = c("gray", four_cat_col[c(3,1)]))
  
  diff_abund_sensitivity_comp_data <-
    diff_abund_table %>% 
    filter(dataset %in% "vOTUs") %>% 
    left_join(read_tsv("data/diff_abund/no_crc/maaslin2_vOTUs_amgs_vs_life_demo_clin.tsv") %>% 
                select(no_crc_coef = coef, no_crc_qval = qval, variable, value, feature)) %>% 
    left_join(meta_variables %>% rename(variable = var_id, var_dat = dataset) %>% select(variable, var_name, var_dat)) %>% 
    group_by(var_name) %>% 
    mutate(sum_sign = sum(qval < 0.05)) %>% 
    ungroup() %>% 
    mutate(var_name = fct_reorder(.f = var_name, .x = sum_sign, .fun = max, .desc = TRUE)) %>% 
    mutate(significance_codes = case_when(qval < 0.05 & no_crc_qval < 0.05 ~ "both significant",
                                          !qval < 0.05 & no_crc_qval < 0.05 ~ "no crc significant",
                                          qval < 0.05 & !no_crc_qval < 0.05 ~ "original significant",
                                          TRUE ~ "none significant") %>% 
             factor(levels = c("none significant", "original significant", "no crc significant", "both significant")))
  
  diff_abund_comp_plot <-
    diff_abund_sensitivity_comp_data %>% 
    filter(significance_codes != "none significant") %>% 
    ggplot(aes(x = coef, y = no_crc_coef, color = significance_codes)) +
    geom_point(data = diff_abund_sensitivity_comp_data %>% filter(significance_codes == "none significant"), alpha = 0.1) +
    geom_point() +
    facet_wrap(~var_name + name, labeller = as_labeller(lab_func), ncol = 4) +
    scale_color_manual(values = c(four_cat_col[1:2], "gray", four_cat_col[3])) +
    theme_bw() +
    labs(x = "Original coefficient", y = "Coefficient without CRC", color = "FDR significance")
  
  host_association_sensitivity_plot <-
    grid.arrange(alpha_sensitivity_comp + theme(legend.position = "none") + labs(subtitle = "a)"),
               beta_sensitivity_comp + theme(legend.position = "none") + labs(subtitle = "b)"),
               diff_abund_comp_plot + labs(subtitle = "c)"),
               layout_matrix = rbind(c(1,3,3), c(2,3,3)))
  ggsave("figures/diversity_composition/host_associations_sensitivity_plot.pdf", plot = host_association_sensitivity_plot, height = 210, width = 270, units = "mm")
  ggsave("figures/diversity_composition/host_associations_sensitivity_plot.png", plot = host_association_sensitivity_plot, height = 210, width = 270, units = "mm")
  
}



