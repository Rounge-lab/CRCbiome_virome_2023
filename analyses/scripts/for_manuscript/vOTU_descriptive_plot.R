

make_vOTU_descriptive_plots <- function() {
   
   ###Taxonomy
   taxa_summary <-
      vOTU_stats %>% 
      create_taxa_summary(min_fam_freq = 20) %>% 
      filter(Classification != "Classified") %>% 
      mutate(Family = factor(Classification, 
                             levels = c(c("Unknown", "Higher order", "Other"),
                                        levels(Classification)[ -which(levels(Classification) %in% c("Unknown", "Higher order", "Other"))])))
   vOTU_by_taxonomy_lifecycle <- 
      read_rds("tables/vOTU_based_stats/vOTU_by_taxonomy_lifecycle.Rds")
   
   lifecycle_by_family <- 
      read_rds("tables/vOTU_based_stats/lifecycle_by_family.Rds")
   
   family_lifecycle_tests <- 
      read_rds("tables/vOTU_based_stats/family_lifecycle_tests.Rds")
   
   family_lifecycle_tests <-
      family_lifecycle_tests %>% 
      mutate(Family_lyt_reorder = factor(Family, levels = (family_lifecycle_tests %>% 
                                                              mutate(fraction_lytic = n_lytic/n_total) %>% 
                                                              filter(!Family %in% c("Higher order", "Unknown")) %>% 
                                                              arrange(desc(fraction_lytic)) %>% 
                                                              pull(Family) %>% 
                                                              as.character() %>% 
                                                              c("Higher order", "Unknown") %>% 
                                                              rev())))
   
   
   ## Summarize taxonomy, lifecycle, abundance, prevalence
   vOTU_by_taxonomy_lifecycle <-
      read_rds("tables/vOTU_based_stats/vOTU_by_taxonomy_lifecycle.Rds")
   
   amgs_per_fam_tests <- 
      read_rds("tables/vOTU_based_stats/family_AMG_tests.Rds")
   
   ## Plot numbers per taxonomy
   taxa_bar_horiz <-
      taxa_summary %>% 
      filter(!Classification %in% "Classified") %>% 
      ggplot(aes(x = n, y = Family, fill = Family)) +
      geom_col(aes(x = n_genomes), fill = "lightgray", width = 0.6) +
      geom_col() +
      scale_fill_manual(values = c(colorRampPalette(colors = c("#27E2AA", "#4D27E2"))(3)[3:2],
                                   "darkgray",
                                   colorRampPalette(colors = c("#BCE227", "#E2275F"))(nlevels(vOTU_by_taxonomy_lifecycle$Family)-3)), guide = "none") +
      scale_x_continuous(trans = "log10") +
      theme_bw() +
      labs(x = "vOTUs (genomes)")
         
   
   ## Genome size per family
   g_sizes <- 
      vOTU_stats %>% 
      select(vOTU, genome_length, genes) %>% 
      left_join(vOTU_by_taxonomy_lifecycle) %>%
      mutate(Family = factor(Family, levels = levels(taxa_summary$Family))) %>% 
      filter(!is.na(Family)) %>%
      mutate(genome_length = genome_length/1000) %>% 
      ggplot(aes(x = genome_length, y = Family, fill = Family)) +
      geom_boxplot(outlier.size = 0.25, lwd = .25) +
      theme_bw() +
      labs(x = "vOTU length (kbp)", fill = "", y = "") +
      scale_x_continuous(trans = "log10") +
      scale_fill_manual(values = c(colorRampPalette(colors = c("#27E2AA", "#4D27E2"))(3)[3:2],
                                   "darkgray",
                                   colorRampPalette(colors = c("#BCE227", "#E2275F"))(nlevels(vOTU_by_taxonomy_lifecycle$Family)-3)), guide = "none") +
      theme(axis.text.y = element_blank())
   
   lifecycle_lytic_plot <-
      family_lifecycle_tests %>% 
      mutate(genomes_per_cluster = cut(log10(n_total), breaks = c(0,1,2,3,5), labels = c("1-10 genomes", "10-100 genomes", "100-1000 genomes", ">1000 genomes"))) %>% 
      mutate(fraction_lytic = n_lytic/n_total) %>% 
      mutate(Family = factor(Family, levels = levels(taxa_summary$Family))) %>% 
      mutate(fraction_lytic = fraction_lytic*100) %>% 
      ggplot(aes(x = fraction_lytic, y = Family, color = Family)) +
      geom_vline(xintercept = family_lifecycle_tests %>% summarize(tot_frac = sum(n_lytic)/sum(n_total)*100) %>% pull(tot_frac), color = "black", linetype = 2) +
      # geom_point(aes(size = genomes_per_cluster)) +
      geom_point(size = 2) +
      scale_color_manual(values = c(colorRampPalette(colors = c("#27E2AA", "#4D27E2"))(3)[3:2],
                                    "darkgray",
                                    colorRampPalette(colors = c("#BCE227", "#E2275F"))(nlevels(vOTU_by_taxonomy_lifecycle$Family)-3)), guide = "none") +
      theme_bw() +
      theme(legend.position = "none",
            legend.background = element_rect(color = "black"),
            legend.title = element_blank()) +
      labs(x = "lytic genomes (%)", y = "", size = "") +
      scale_x_continuous(limits = c(0,100)) +
      scale_size_manual(values = c(0.75,2,4))
   
   ## Gene annotation per family
   fraction_annotated <-
      vOTU_stats %>% 
      mutate(perc_annotated = fraction_annotated*100) %>% 
      select(vOTU, perc_annotated, genes) %>% 
      left_join(vOTU_by_taxonomy_lifecycle) %>%
      mutate(Family = factor(Family, levels = levels(taxa_summary$Family))) %>% 
      filter(!is.na(Family),
             !is.na(perc_annotated)) %>% 
      ggplot(aes(x = perc_annotated, y = Family, fill = Family)) +
      theme_bw() +
      geom_boxplot(outlier.size = 0.25, lwd =.25) +
      labs(x = "vOTU gene annotation (%)", fill = "", y = "") +
      scale_fill_manual(values = c(colorRampPalette(colors = c("#27E2AA", "#4D27E2"))(3)[3:2],
                                   "darkgray",
                                   colorRampPalette(colors = c("#BCE227", "#E2275F"))(nlevels(vOTU_by_taxonomy_lifecycle$Family)-3)), guide = "none") +
      theme(axis.text.y = element_blank())
   
   ## AMG fraction per family
   amgs_per_family_plot <-
      amgs_per_fam_tests %>%
      mutate(significance = case_when(p_corr < 0.001 ~ "***",
                                      p_corr < 0.01 ~ "**",
                                      p_corr < 0.05 ~ "*",
                                      TRUE ~ ""),
             direction = case_when(estimate >= 1 ~ "positive",
                                   estimate < 1 ~ "negative")) %>% 
      mutate(Family_lyt_reorder = factor(Family, levels = levels(family_lifecycle_tests$Family_lyt_reorder))) %>% 
      mutate(Family = factor(Family, levels = levels(taxa_summary$Family))) %>% 
      mutate(cat = factor(cat, levels = amgs_per_fam_tests %>% 
                             select(Family, cat, with_cat) %>% 
                             pivot_wider(names_from = cat, values_from = with_cat, values_fill = 0) %>% 
                             column_to_rownames("Family") %>% 
                             order_by_cluster_columns())) %>% 
      # ggplot(aes(y = Family, x = cat, fill = with_cat)) +
      ggplot(aes(y = Family, x = cat, fill = with_cat, label = significance)) +
      geom_tile() +
      geom_text(size = 5*0.36) +
      theme_bw() +
      labs(x = "", fill = "", y = "") +
      theme(panel.grid = element_blank(),
            line = element_line(linewidth = 0.5)) +
      # theme(axis.text.x = element_text(angle = 40, hjust = 1),
      #       panel.grid = element_blank(),
      #       panel.border = element_blank()) +
      scale_fill_gradient(low = "white", high = "#3A84C5", limits = c(0,1))
   
   prevalence_abundance_plot <-
      vOTU_by_taxonomy_lifecycle %>% 
      filter(cluster_contigs > 2,
             !Family %in% "Unknown") %>%
      # filter(cluster_contigs > 4) %>%
      ggplot(aes(x = prevalence, y = mean_present)) +
      geom_point(aes(color = Family), size = 1) +
      geom_density2d(data = vOTU_by_taxonomy_lifecycle %>% filter(cluster_contigs > 2) %>% select(-Family), aes(x = prevalence, y = mean_present), color = "gray", alpha = 0.5) +
      scale_x_continuous(trans = "log10", labels = function(x) sprintf('%.3g', x)) +
      scale_y_continuous(trans = "log10") +
      scale_color_manual(values = c(rev(colorRampPalette(colors = c("#BCE227", "#E2275F"))(nlevels(vOTU_by_taxonomy_lifecycle$Family)-3)),
                                    "darkgray",
                                    colorRampPalette(colors = c("#27E2AA", "#4D27E2"))(3)[2])) +
      theme_bw() +
      theme(legend.position = "none") +
      labs(y = "mean prevalent abundance", x = "prevalence") +
      theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1)) +
      facet_wrap(~Family, nrow = 2)
   
   
   vOTU_fig <- 
      gridExtra::grid.arrange(taxa_bar_horiz + 
                                 theme(legend.key.size = unit(.3, "cm"), 
                                       text = element_text(size = 8)) + 
                                 labs(subtitle = "a)"), 
                              g_sizes + 
                                 theme_bw() +
                                 theme(text = element_text(size = 8),
                                       axis.text.y = element_blank()) +
                                 labs(subtitle = "b)"),
                              lifecycle_lytic_plot + 
                                 theme(text = element_text(size = 8),
                                       axis.text.y = element_blank()) +
                                 labs(y = "", subtitle = "c)"),
                              fraction_annotated +
                                 theme(text = element_text(size = 8),
                                       axis.text.y = element_blank()) +
                                 labs(y = "", subtitle = "d)"),
                              amgs_per_family_plot +
                                 theme(axis.text.x = element_text(),
                                       axis.text.y = element_blank(),
                                       text = element_text(size = 8),
                                       legend.key.size = unit(.3, "cm")) +
                                 labs(subtitle = "e)"),
                              prevalence_abundance_plot + 
                                 theme(legend.position = "none", 
                                       text = element_text(size = 8)) + 
                                 labs(subtitle = "f)"), 
                              layout_matrix = rbind(c(1,1,1,2,2,3,3,4,4,5,5,5),
                                                    c(6,6,6,6,6,6,6,6,6,6,6,6))) 
   vOTU_fig %>% 
      ggsave(filename = "figures/vOTU_summary/vOTU_summary_fig.png", plot = ., height = 170/9*8, width = 170, units = "mm")
   vOTU_fig %>% 
      ggsave(filename = "figures/vOTU_summary/vOTU_summary_fig.pdf", plot = ., height = 170/9*8, width = 170, units = "mm")
}


create_amg_plots <- function() {
   
   vOTU_by_taxonomy_lifecycle <- 
      read_rds("tables/vOTU_based_stats/vOTU_by_taxonomy_lifecycle.Rds")
   
   ## AMGs by family
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
   
   
   amgs_per_family_plot <-
      amgs_per_fam_frac %>%
      mutate(Family = factor(Family, levels = amgs_per_fam_frac %>% 
                                select(-starts_with("sum")) %>% 
                                pivot_wider(names_from = Family, values_from = with_cat, values_fill = 0) %>% 
                                column_to_rownames("cat") %>% 
                                order_by_cluster_columns())) %>% 
      mutate(cat = factor(cat, levels = amgs_per_fam_frac %>% 
                             select(-starts_with("sum")) %>% 
                             pivot_wider(names_from = cat, values_from = with_cat, values_fill = 0) %>% 
                             column_to_rownames("Family") %>% 
                             order_by_cluster_columns())) %>% 
      ggplot(aes(y = Family, x = cat, fill = with_cat)) +
      geom_tile() +
      theme_bw() +
      labs(x = "", fill = "", y = "") +
      theme(axis.text.x = element_text(angle = 40, hjust = 1),
            panel.grid = element_blank(),
            panel.border = element_blank()) +
      scale_fill_gradient(low = "#EFF1F1", high = "#4D27E2", limits = c(0,1))
   
   amgs_per_family_plot %>% 
      ggsave("figures/vOTU_summary/amgs_per_family_heat.png", plot = ., height = 4, width = 3.5)
   
   
   amgs_MISC_module_by_vOTU <- 
      read_tsv("data/amgs/amg_MISC_module_by_vOTU_mid_min75.tsv")
   
   amgs_MISC_per_fam_frac <-
      amgs_MISC_module_by_vOTU %>% 
      # rename(`Carbon Utilization` = `carbon utilization`) %>%
      pivot_longer(-vOTU, names_to = "cat", values_to = "presence") %>% 
      filter(!str_detect(cat, "Wood")) %>% 
      left_join(vOTU_by_taxonomy_lifecycle %>% select(vOTU, cluster_contigs, Family)) %>% 
      filter(!is.na(Family)) %>% 
      group_by(Family, cat) %>% 
      summarize(with_cat = sum(presence)/n(),
                sum_with = sum(presence),
                sum_without = sum(!presence)) %>% 
      ungroup()
   
   amgs_MISC_per_fam_frac
   
   amgs_MISC_per_family_plot <-
      amgs_MISC_per_fam_frac %>%
      mutate(Family = factor(Family, levels = amgs_MISC_per_fam_frac %>% 
                                select(-starts_with("sum")) %>% 
                                pivot_wider(names_from = Family, values_from = with_cat, values_fill = 0) %>% 
                                column_to_rownames("cat") %>% 
                                order_by_cluster_columns())) %>% 
      mutate(cat = factor(cat, levels = amgs_MISC_per_fam_frac %>% 
                             select(-starts_with("sum")) %>% 
                             pivot_wider(names_from = cat, values_from = with_cat, values_fill = 0) %>% 
                             column_to_rownames("Family") %>% 
                             order_by_cluster_columns())) %>% 
      ggplot(aes(x = Family, y = cat, fill = with_cat, label = signif(with_cat, 2))) +
      geom_tile() +
      geom_text() +
      theme_bw() +
      labs(x = "", fill = "", y = "") +
      theme(axis.text.x = element_text(angle = 40, hjust = 1),
            panel.grid = element_blank(),
            panel.border = element_blank()) +
      scale_fill_gradient(low = "#EFF1F1", high = "#4D27E2", trans = "log10")
   
   upset_plot <- 
      amgs_cat_per_vOTU %>% 
      select(!`carbon utilization (Woodcroft)`) %>%
      rename(`Carbon Utilization` = `carbon utilization`) %>%
      mutate(across(-starts_with("vOTU"), .fns = function(x) x*1)) %>% 
      as.data.frame() %>% 
      column_to_rownames("vOTU") %>% 
      upset(
         nsets = 6, 
         order.by = "freq",
         matrix.color =  col_pal_distinct[2],
         main.bar.color = col_pal_distinct[3],
         sets.bar.color = col_pal_distinct[6])
   
   pdf("figures/vOTU_summary/upset_amg_cat.pdf")
   upset_plot
   dev.off()
   
   
   
}




# (family_lifecycle_distribution_plot <-
#       lifecycle_by_family %>% 
#       ggplot(aes(x = n_total, y = fraction_lyso, label = Family, color = Family)) +
#       geom_text_repel(box.padding = .5,
#                       max.overlaps = Inf,
#                       xlim = c(NA,log10(20)),
#                       size = 3) +
#       geom_point() +
#       scale_x_continuous(trans = "log10", limits = c(1,NA)) +
#       scale_y_continuous(limits = c(0,1)) +
#       scale_color_manual(values = c(colorRampPalette(colors = c("#BCE227", "#E2275F"))(nlevels(vOTU_by_taxonomy_lifecycle$Family)-3),
#                                     "gray",
#                                     colorRampPalette(colors = c("#27E2AA", "#4D27E2"))(3)[2:3])) +
#       theme_bw() +
#       labs(x = "viral contigs",
#            y = "fraction lysogenic"))
# 
# 
# lifecycle_lysogenic_OR_plot <- 
#    family_lifecycle_tests %>% 
#    ggplot(aes(x = estimate, xmin = conf.low, xmax = conf.high, y = Family, color = Family)) +
#    geom_vline(xintercept = 1, color = "black", linetype = 2) +
#    geom_pointrange() +
#    scale_x_continuous(trans = "log10") +
#    scale_color_manual(values = rev(c(colorRampPalette(colors = c("#BCE227", "#E2275F"))(nlevels(vOTU_by_taxonomy_lifecycle$Family)-3),
#                                      "gray",
#                                      colorRampPalette(colors = c("#27E2AA", "#4D27E2"))(3)[2:3]))) +
#    theme_bw() +
#    theme(legend.position = "none") +
#    labs(x = "OR (lysogenic lifecycle for family versus rest)")
# 
# 
# 
# ## Barplot
# (amgs_per_fam_bar <-
#       amgs_per_fam_frac %>% 
#       select(-c(sum_with, sum_without)) %>% 
#       mutate(without_cat = 1-with_cat) %>% 
#       pivot_longer(c(with_cat, without_cat), names_to = "pres_abs", values_to = "fraction") %>% 
#       ggplot(aes(x = fraction, y = Family, fill = pres_abs)) +
#       geom_col() +
#       facet_wrap(~cat, nrow = 1) +
#       theme_bw() +
#       scale_fill_manual(values = unname(col_pal_distinct[c(3,2)])))
# 
# ## Heatmap
# (amgs_per_fam_heat <-
#       amgs_per_fam_frac %>% 
#       ggplot(aes(x = cat, y = Family, fill = with_cat)) +
#       geom_tile() +
#       # facet_wrap(~cat, nrow = 1) +
#       theme_bw() +
#       # scale_fill_continuous(values = unname(col_pal_distinct[c(3,2)])) +
#       theme(axis.text.x = element_text(angle = 30, hjust = 1)))
# 
# ## AMGs by lifecycle
# amgs_per_lifecyc <-
#    amgs_cat_per_vOTU %>% 
#    pivot_longer(-vOTU, names_to = "cat", values_to = "presence") %>% 
#    filter(!str_detect(cat, "Wood")) %>% 
#    left_join(vOTU_by_taxonomy_lifecycle %>% select(vOTU, cluster_contigs, lifecycle_propensity_adj)) %>% 
#    filter(!is.na(lifecycle_propensity_adj)) %>% 
#    group_by(lifecycle_propensity_adj, cat) %>% 
#    summarize(with_cat = sum(presence)/n(),
#              sum_with = sum(presence),
#              sum_without = sum(!presence)) %>% 
#    ungroup()
# 
# ## Barplot
# (amgs_per_lifecyc_bar <- 
#       amgs_per_lifecyc %>% 
#       select(-c(sum_with, sum_without)) %>% 
#       mutate(without_cat = 1-with_cat) %>% 
#       pivot_longer(c(with_cat, without_cat), names_to = "pres_abs", values_to = "fraction") %>% 
#       ggplot(aes(x = fraction, y = lifecycle_propensity_adj, fill = pres_abs)) +
#       geom_col() +
#       facet_wrap(~cat, nrow = 1) +
#       theme_bw() +
#       labs(y = "lifecycle propensity") +
#       scale_fill_manual(values = unname(col_pal_distinct[c(3,2)])))
# 
# ## Heatmap
# (amgs_per_lifecyc_heat <- 
#       amgs_per_lifecyc %>% 
#       ggplot(aes(x = cat, y = lifecycle_propensity_adj, fill = with_cat)) +
#       geom_tile() +
#       # facet_wrap(~cat, nrow = 1) +
#       theme_bw() +
#       # scale_fill_continuous(values = unname(col_pal_distinct[c(3,2)])) +
#       labs(y = "lifecycle propensity") +
#       theme(axis.text.x = element_text(angle = 30, hjust = 1)))
# 
# 
# 

 
   



# ## AMGs, family, lifecycle
# (amg_per_lifecycle_taxonomy_plot <- 
#       amgs_per_lifecyc %>% 
#       mutate(var = "lifecycle") %>% 
#       rename(group = lifecycle_propensity_adj) %>% 
#       bind_rows(amgs_per_fam_frac %>%  
#                    mutate(var = "family") %>% 
#                    rename(group = Family)) %>% 
#       ggplot(aes(x = group, y = cat, fill = with_cat)) +
#       geom_tile() +
#       facet_wrap(~var, nrow = 1, scales = "free_x") +
#       theme_bw() +
#       labs(x = "", fill = "") +
#       theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
#       scale_fill_gradient(low = "#EFF1F1", high = "#4D27E2", limits = c(0,1)))
# 
# amg_per_lifecycle_taxonomy_plot %>% 
#    ggsave(filename = "figures/vOTU_summary/amg_cats_by_fam_and_lifecyc.png", plot = ., height = 3, width = 8.5)
# amg_per_lifecycle_taxonomy_plot %>% 
#    ggsave(filename = "figures/vOTU_summary/amg_cats_by_fam_and_lifecyc.png", plot = ., height = 3, width = 8.5)





# vOTU_by_taxonomy_lifecycle_test %>%
#    pivot_wider(names_from = lifecycle_propensity_adj, values_from = n, values_fill = 0) %>% 
#    left_join(family_lifecycle_tests %>% select(Family, p.value))
# 
# vOTU_by_taxonomy_lifecycle %>% 
#    filter(Family %in% "Steigviridae") %>% 
#    left_join(vOTU_stats, by = "vOTU") %>% 
#    View()
# 
# vOTU_by_taxonomy_lifecycle %>% 
#    filter(Family %in% "Intestiviridae") %>% 
#    left_join(vOTU_stats, by = "vOTU") %>% 
#    View()
# 
# vOTU_by_taxonomy_lifecycle %>% 
#    filter(Family %in% "Suoliviridae") %>% 
#    left_join(vOTU_stats, by = "vOTU") %>% 
#    View()

#chisq.test (vOTU_by_taxonomy_lifecycle$Family, vOTU_by_taxonomy_lifecycle$lifecycle_propensity_adj)


# vOTU_stats %>% 
#    group_by(Family) %>% 
#    summarize(n_provirus = sum(cluster_proviruses),
#              n_not_provirus = sum(cluster_not_proviruses)) %>% 
#    mutate(fraction = n_provirus/(n_provirus+n_not_provirus)) %>% 
#    as.data.frame()
# 
# vOTU_lifecycle_cat %>% 
#    # filter(cluster_contigs >= 30) %>% 
#    ggplot(aes(x = cluster_not_proviruses, cluster_proviruses, color = lifecycle_propensity_adj)) +
#    geom_point()
# 
# vOTU_lifecycle_cat %>% 
#    #count(lifecycle_propensity_unadj, lifecycle_propensity_adj) 
#    ggplot(aes(x = lifecycle_propensity_adj, cluster_contigs, fill = lifecycle_propensity_adj)) +
#    theme_bw() +
#    scale_y_continuous(trans="log10")+
#    geom_boxplot()



# 
# 
# make_vOTU_descriptive_plots <- function() {
#    
#    ###Taxonomy
#    taxa_summary <-
#       vOTU_stats %>% 
#       create_taxa_summary(min_fam_freq = 20)
#    
#    vOTU_by_taxonomy_lifecycle <- 
#       read_rds("tables/vOTU_based_stats/vOTU_by_taxonomy_lifecycle.Rds")
#    
#    lifecycle_by_family <- 
#       read_rds("tables/vOTU_based_stats/lifecycle_by_family.Rds")
#    
#    family_lifecycle_tests <- 
#       read_rds("tables/vOTU_based_stats/family_lifecycle_tests.Rds")
#    
#    family_lifecycle_tests <-
#       family_lifecycle_tests %>% 
#       mutate(Family_lyt_reorder = factor(Family, levels = (family_lifecycle_tests %>% 
#                                                               mutate(fraction_lytic = n_lytic/n_total) %>% 
#                                                               filter(!Family %in% c("Higher order", "Unknown")) %>% 
#                                                               arrange(desc(fraction_lytic)) %>% 
#                                                               pull(Family) %>% 
#                                                               as.character() %>% 
#                                                               c("Higher order", "Unknown") %>% 
#                                                               rev())))
#    
#    taxa_family_classified <-
#       taxa_summary %>% 
#       ggplot(aes(x = level, y = rel_n, fill = Classification)) +
#       geom_col(width = 0.4, color = "black")+
#       labs(x = "", 
#            y = "vOTU fraction") +
#       scale_fill_manual(values = c("gray", 
#                                    colorRampPalette(colors = c("#E2275F", "#BCE227"))(nlevels(taxa_summary$Classification)-4),
#                                    colorRampPalette(colors = c("#27E2AA", "#4D27E2"))(3))) +
#       theme_bw()
#    
#    
#    ## Summarize taxonomy, lifecycle, abundance, prevalence
#    vOTU_by_taxonomy_lifecycle <-
#       read_rds("tables/vOTU_based_stats/vOTU_by_taxonomy_lifecycle.Rds")
#    
#    prevalence_abundance_plot <-
#       vOTU_by_taxonomy_lifecycle %>% 
#       filter(cluster_contigs > 2,
#              !Family %in% "Unknown") %>%
#       # filter(cluster_contigs > 4) %>%
#       ggplot(aes(x = prevalence, y = mean_present)) +
#       geom_point(aes(color = Family), size = 1) +
#       geom_density2d(data = vOTU_by_taxonomy_lifecycle %>% filter(cluster_contigs > 2) %>% select(-Family), aes(x = prevalence, y = mean_present), color = "gray", alpha = 0.5) +
#       scale_x_continuous(trans = "log10", labels = function(x) sprintf('%.3g', x)) +
#       scale_y_continuous(trans = "log10") +
#       scale_color_manual(values = c(colorRampPalette(colors = c("#BCE227", "#E2275F"))(nlevels(vOTU_by_taxonomy_lifecycle$Family)-3),
#                                     "gray",
#                                     colorRampPalette(colors = c("#27E2AA", "#4D27E2"))(3)[2])) +
#       theme_bw() +
#       theme(legend.position = "none") +
#       labs(y = "mean prevalent abundance", x = "prevalence") +
#       theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1)) +
#       facet_wrap(~Family, ncol = 3)
#    
#    # ggsave("figures/vOTU_summary/mean_prevalent_abundance_by_family.png", plot = prevalence_abundance_plot, height = 4.5, width = 5)
#    
#    lifecycle_lytic_plot <-
#       family_lifecycle_tests %>% 
#       mutate(genomes_per_cluster = cut(log10(n_total), breaks = c(0,1,2,3,5), labels = c("1-10 genomes", "10-100 genomes", "100-1000 genomes", ">1000 genomes"))) %>% 
#       mutate(fraction_lytic = n_lytic/n_total) %>% 
#       ggplot(aes(x = fraction_lytic, y = Family_lyt_reorder, color = Family)) +
#       geom_vline(xintercept = family_lifecycle_tests %>% summarize(tot_frac = sum(n_lytic)/sum(n_total)) %>% pull(tot_frac), color = "black", linetype = 2) +
#       geom_point(aes(size = genomes_per_cluster)) +
#       scale_color_manual(values = rev(c(colorRampPalette(colors = c("#BCE227", "#E2275F"))(nlevels(vOTU_by_taxonomy_lifecycle$Family)-3),
#                                         "gray",
#                                         colorRampPalette(colors = c("#27E2AA", "#4D27E2"))(3)[2:3])), guide = "none") +
#       theme_bw() +
#       # theme(legend.position = c(0.3,0.75),
#       theme(legend.position = "none",
#             legend.background = element_rect(color = "black"),
#             legend.title = element_blank()) +
#       labs(x = "lytic genome fraction", y = "", size = "") +
#       scale_x_continuous(limits = c(0,1)) +
#       scale_size_manual(values = c(0.75,2,4))
#    
#    # ggsave("figures/vOTU_summary/lifecycle_by_family.png", plot = family_lifecycle_distribution_plot, height = 4, width = 5.5)
#    
#    amgs_per_fam_tests <- read_rds("tables/vOTU_based_stats/family_AMG_tests.Rds")
#    
#    amgs_per_family_plot <-
#       amgs_per_fam_tests %>%
#       mutate(significance = case_when(p_corr < 0.001 ~ "***",
#                                       p_corr < 0.01 ~ "**",
#                                       p_corr < 0.05 ~ "*",
#                                       TRUE ~ ""),
#              direction = case_when(estimate >= 1 ~ "positive",
#                                    estimate < 1 ~ "negative")) %>% 
#       mutate(Family_lyt_reorder = factor(Family, levels = levels(family_lifecycle_tests$Family_lyt_reorder))) %>% 
#       mutate(Family = factor(Family, levels = amgs_per_fam_tests %>% 
#                                 select(Family, cat, with_cat) %>% 
#                                 pivot_wider(names_from = Family, values_from = with_cat, values_fill = 0) %>% 
#                                 column_to_rownames("cat") %>% 
#                                 order_by_cluster_columns())) %>% 
#       mutate(cat = factor(cat, levels = amgs_per_fam_tests %>% 
#                                 select(Family, cat, with_cat) %>% 
#                              pivot_wider(names_from = cat, values_from = with_cat, values_fill = 0) %>% 
#                              column_to_rownames("Family") %>% 
#                              order_by_cluster_columns())) %>% 
#       # ggplot(aes(y = Family, x = cat, fill = with_cat)) +
#       ggplot(aes(y = Family_lyt_reorder, x = cat, fill = with_cat, label = significance)) +
#       geom_tile() +
#       geom_text(size = 5*0.36) +
#       theme_bw() +
#       labs(x = "", fill = "", y = "") +
#       theme(panel.grid = element_blank()) +
#       # theme(axis.text.x = element_text(angle = 40, hjust = 1),
#       #       panel.grid = element_blank(),
#       #       panel.border = element_blank()) +
#       scale_fill_gradient(low = "white", high = "#4D27E2", limits = c(0,1))
#    
#    g_sizes <- 
#       vOTU_stats %>% 
#       # filter(completeness > 90) %>% 
#       select(vOTU, genome_length, genes) %>% 
#       left_join(vOTU_by_taxonomy_lifecycle) %>%
#       left_join(family_lifecycle_tests %>% select(Family, Family_lyt_reorder)) %>% 
#       filter(!is.na(Family_lyt_reorder)) %>% 
#       ggplot(aes(x = genome_length, y = Family, fill = Family)) +
#       geom_boxplot(outlier.size = 0.25) +
#       theme_bw() +
#       labs(x = "genome length", fill = "", y = "") +
#       scale_x_continuous(trans = "log10") +
#       scale_fill_manual(values = c("gray", 
#                                    colorRampPalette(colors = c("#E2275F", "#BCE227"))(nlevels(taxa_summary$Classification)-4),
#                                    colorRampPalette(colors = c("#27E2AA", "#4D27E2"))(3)), guide = "none")# +
#       theme(axis.text.y = element_blank())
#    
#    fraction_annotated <- 
#       vOTU_stats %>% 
#       # filter(completeness > 90) %>% 
#       select(vOTU, fraction_annotated, genes) %>% 
#       left_join(vOTU_by_taxonomy_lifecycle) %>%
#       left_join(family_lifecycle_tests %>% select(Family, Family_lyt_reorder)) %>% 
#       filter(!is.na(Family_lyt_reorder),
#              !is.na(fraction_annotated)) %>% 
#       ggplot(aes(x = fraction_annotated, y = Family_lyt_reorder, fill = Family)) +
#       geom_boxplot(outlier.size = 0.25) +
#       theme_bw() +
#       labs(x = "fraction of genes annotated", fill = "", y = "") +
#       # scale_x_continuous(trans = "log10") +
#       scale_fill_manual(values = c(colorRampPalette(colors = c("#BCE227", "#E2275F"))(nlevels(vOTU_by_taxonomy_lifecycle$Family)-3),
#                                    "gray",
#                                    colorRampPalette(colors = c("#27E2AA", "#4D27E2"))(3)[2:3]), guide = "none") +
#       theme(axis.text.y = element_blank())
#    
#    taxa_bar_horiz <-
#       taxa_summary %>% 
#       filter(!Classification %in% "Classified") %>% 
#       ggplot(aes(x = n, y = Classification, fill = Classification)) +
#       geom_col() +
#       scale_fill_manual(values = c(colorRampPalette(colors = c("#BCE227", "#E2275F"))(nlevels(vOTU_by_taxonomy_lifecycle$Family)-3),
#                                    "gray",
#                                    colorRampPalette(colors = c("#27E2AA", "#4D27E2"))(3)[2:3]), guide = "none") +
#       scale_x_continuous(trans = "log10")
#    
#    ## vOTU_Figure
#    # vOTU_fig <- 
#    #    gridExtra::grid.arrange(taxa_family_classified + theme(legend.key.size = unit(.4, "cm"), legend.text = element_text(size = 8)) + labs(subtitle = "a)"), 
#    #                            prevalence_abundance_plot + theme(legend.position = "none") + labs(subtitle = "c)"), 
#    #                            lifecycle_lytic_plot + scale_y_discrete(position = "right") + theme(plot.margin = margin(l = 50)) + labs(y = "", subtitle = "b)"),
#    #                            layout_matrix = cbind(c(1,1,3,3),
#    #                                                  c(2,2,2,2),
#    #                                                  c(2,2,2,2))) 
#    # vOTU_fig %>% 
#    #    ggsave(filename = "figures/vOTU_summary/vOTU_summary_fig.png", plot = ., height = 6, width = 9)
#    # vOTU_fig %>% 
#    #    ggsave(filename = "figures/vOTU_summary/vOTU_summary_fig.pdf", plot = ., height = 6, width = 9)
#    # # vOTU_fig %>% 
#    # #    ggsave(filename = "figures/vOTU_summary/vOTU_summary_fig.eps", plot = ., height = 6, width = 9)
#    # 
#    vOTU_fig <- 
#       gridExtra::grid.arrange(taxa_family_classified + 
#                                  theme(legend.key.size = unit(.3, "cm"), 
#                                        text = element_text(size = 8)) + 
#                                  labs(subtitle = "a)"), 
#                               prevalence_abundance_plot + 
#                                  theme(legend.position = "none", 
#                                        text = element_text(size = 8)) + 
#                                  labs(subtitle = "d)"), 
#                               lifecycle_lytic_plot + 
#                                  scale_y_discrete(position = "right") + 
#                                  theme(plot.margin = margin(l = 20), 
#                                        text = element_text(size = 8)) + 
#                                  labs(y = "", subtitle = "b)"),
#                               amgs_per_family_plot +
#                                  theme(axis.text.x = element_text(),
#                                        axis.text.y = element_blank(),
#                                        text = element_text(size = 8),
#                                        panel.border = element_rect(color = "black", fill = NA)) +
#                                  labs(subtitle = "c)"),
#                               layout_matrix = rbind(c(1,1,2,2),
#                                                     c(3,4,2,2))) 
#    vOTU_fig %>% 
#       ggsave(filename = "figures/vOTU_summary/vOTU_summary_fig_1.png", plot = ., height = 170/9*6, width = 170, units = "mm")
#    vOTU_fig %>% 
#       ggsave(filename = "figures/vOTU_summary/vOTU_summary_fig.pdf", plot = ., height = 6, width = 9)
#    
#    vOTU_fig <- 
#       gridExtra::grid.arrange(taxa_bar_horiz + 
#                                  theme(legend.key.size = unit(.3, "cm"), 
#                                        text = element_text(size = 8)) + 
#                                  labs(subtitle = "a)"), 
#                               g_sizes + 
#                                  theme_bw() +
#                                  theme(text = element_text(size = 8),
#                                        axis.text.y = element_blank()) +
#                                  labs(subtitle = "b)"),
#                               lifecycle_lytic_plot + 
#                                  theme(text = element_text(size = 8),
#                                        axis.text.y = element_blank()) +
#                                  labs(y = "", subtitle = "c)"),
#                               fraction_annotated +
#                                  theme(text = element_text(size = 8),
#                                        axis.text.y = element_blank()) +
#                                  labs(y = "", subtitle = "d)"),
#                               amgs_per_family_plot +
#                                  theme(axis.text.x = element_text(),
#                                        axis.text.y = element_blank(),
#                                        text = element_text(size = 8),
#                                        legend.key.size = unit(.3, "cm"),
#                                        panel.border = element_rect(color = "black", fill = NA)) +
#                                  labs(subtitle = "e)"),
#                               prevalence_abundance_plot + 
#                                  facet_wrap(~Family, nrow = 2) +
#                                  theme(legend.position = "none", 
#                                        text = element_text(size = 8)) + 
#                                  labs(subtitle = "f)"), 
#                               layout_matrix = rbind(c(1,2,3,4,5),
#                                                     c(6,6,6,6,6))) 
#    vOTU_fig %>% 
#       ggsave(filename = "figures/vOTU_summary/vOTU_summary_fig_1.png", plot = ., height = 170/9*6, width = 170, units = "mm")
#    vOTU_fig %>% 
#       ggsave(filename = "figures/vOTU_summary/vOTU_summary_fig.pdf", plot = ., height = 6, width = 9)
# }