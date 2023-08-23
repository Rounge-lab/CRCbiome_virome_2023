


create_technical_fig <- function() {
  
  
  # fraction_viral_contigs <- read_tsv("data/abundance_tables/fraction_viral_contigs.tsv")
  
  summary_data_by_sample <-
    sample_data %>% 
    # select(sample_id, sampling_to_extraction, qubit_total_dna, Total_Reads_QC_ATLAS, n_contigs, N50) %>% 
    select(sample_id, sampling_to_extraction, nanodrop_total_dna, Total_Reads_QC_ATLAS, n_contigs, N50) %>% 
    left_join(checkV_res_all_by_sample %>% 
                mutate(viral_contigs = n_complete+n_high_qual+n_med_qual) %>% 
                select(viral_contigs, sample_id) %>% 
                # select(putative_viral_scaffolds, viral_contigs, sample_id) %>% 
                mutate(sample_id = gsub("-", "_", sample_id))) %>% 
    # left_join(fraction_viral_contigs %>% select(sample_id, fraction_viral_contigs)) %>% 
    left_join(virus_abundance %>% 
                pivot_longer(-sample_id, names_to = "ID") %>% 
                group_by(sample_id) %>% 
                summarize(n_vOTU = sum(value>0),
                          invsimpson = vegan::diversity(value, index = "invsimpson")) %>% 
                mutate(sample_id = gsub("-", "_", sample_id)))
  
  
  technical_measures_corrs <-
    lapply(seq(ncol(summary_data_by_sample[,-1])), function(i) {
      lapply(seq(ncol(summary_data_by_sample[,-1])), function(ii) {
        tmp <- 
          cor.test(summary_data_by_sample %>% pull(i+1),
                   summary_data_by_sample %>% pull(ii+1),
                   use = "complete", method = "spearman") %>% 
          tidy() %>% 
          mutate(var_1 = names(summary_data_by_sample)[i+1],
                 var_2 = names(summary_data_by_sample)[ii+1])
      }) %>% 
        bind_rows()
    }) %>% 
    bind_rows() %>% 
    mutate(across(.cols = c(var_1, var_2), 
                  .fns = function(x) factor(x, 
                                            levels = names(summary_data_by_sample)[-1],
                                            labels = c("Storage time", "DNA concentration", "Sequencing depth", 
                                                       "Metagenome contigs", "Assembly N50", "Viral genomes", 
                                                       "Observed vOTUs", "Alpha diversity")))) %>% 
    mutate(significance = case_when(p.value < 0.001 ~ "***",
                                    p.value < 0.01 ~ "**",
                                    p.value < 0.05 ~ "*",
                                    TRUE ~ ""))
  
  
  # sample_data %>% 
  #   pivot_longer(c(nanodrop_total_dna, qubit_total_dna)) %>% 
  #   ggplot(aes(x = factor(ekstraksjonsbatch), y = value, fill = name)) +
  #   geom_boxplot()
  # 
  # extr_data %>% 
  #   mutate(sequenced = paste("S_", rekv_nr, sep = "") %in% sample_data$sample_id) %>% 
  #   pivot_longer(c(nanodrop_total_dna, qubit_total_dna)) %>% 
  #   ggplot(aes(x = factor(ekstraksjonsbatch), y = value, fill = sequenced)) +
  #   geom_boxplot() +
  #   facet_wrap(~name, ncol = 1)
  
  technical_cor_plot <-
    technical_measures_corrs %>% 
    mutate(var_n1 = as.integer(var_1),
           var_n2 = as.integer(var_2)) %>%
    filter(var_n1 < var_n2) %>% 
    ggplot(aes(x = var_1, y = var_2, fill = estimate, label = round(estimate, digits = 2))) +
    geom_tile() +
    geom_text(size = 3) +
    theme_bw() +
    # scale_fill_viridis_c(name = "rho") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          panel.grid = element_blank()) +
    labs(x = "",
         y = "") +
    scale_fill_gradient2(name = "", low = "#6b519d", mid = "white", high =  "#5fc0bf", midpoint = 0, limits = c(-1,1)) +
    coord_fixed()
  
  ggsave("figures/misc/technical_cor_plot.png", plot = technical_cor_plot, width = 6, height = 5.5)
  ggsave("figures/misc/technical_cor_plot.pdf", plot = technical_cor_plot, width = 6, height = 5.5)
  
  ###Boxplot of contigs by stage
  
  contigs_by_stage_plot <-
    checkV_res_all_by_sample %>%
    mutate(viral_filtered = n_complete + n_high_qual+n_med_qual) %>%
    select(sample_id, n_contigs, putative_viral_scaffolds, viral_filtered) %>%
    left_join(virus_abundance %>% pivot_longer(-sample_id) %>% group_by(sample_id) %>% summarize(observed75 = sum(value > 0))) %>% 
    pivot_longer(-sample_id) %>%
    mutate(name = factor(name, 
                         levels = rev(c("n_contigs","putative_viral_scaffolds", "viral_filtered", "observed75")), 
                         labels = rev(c("Metagenome contigs", "Putative viral genomes", "Viral genomes", "Observed vOTUs")))) %>%
    ggplot(aes(x = value, fill = name, y = name)) +
    scale_fill_manual(name="Contigs", labels = c("Metagenome contigs", "Putative viral genomes", "Viral genomes", "Observed vOTUs"), 
                      values = rev(c("#52a1b5","#5fc0bf","#9791c6","#814591"))) +
    geom_boxplot() +
    theme_bw() +
    scale_x_continuous(trans="log10") +
    theme(legend.position = "none") +
    # theme(legend.position = "none",
    #       axis.text.y = element_text(size = 8, angle = 40, hjust = 1)) +
    labs(y = "",
         x = "contigs/genomes per sample")
  
  technical_data_summary_plot <- 
    summary_data_by_sample %>% 
    pivot_longer(-sample_id) %>% 
    mutate(name = factor(name, 
                         levels = names(summary_data_by_sample)[-1],
                         labels = c("Storage time", "DNA concentration", "Sequencing depth", 
                                    "Metagenome contigs", "Assembly N50", "Viral genomes", 
                                    "Observed vOTUs", "Alpha diversity"))) %>% 
    ggplot(aes(x = value, fill = name)) +
    geom_histogram() +
    theme_bw() +
    facet_wrap(~name, scales = "free", ncol = 1) +
    labs(x = "", y = "number of samples") +
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.grid.major.y = element_blank()) +
    scale_fill_manual(values = get_theme_gradient(low_col = "turquoise", n_cols = ncol(summary_data_by_sample)-1))
  
  
  combined_technical_plot <-
    grid.arrange(technical_data_summary_plot + labs(subtitle = "a)"),
                 contigs_by_stage_plot + labs(subtitle = "b)") +theme(plot.margin = margin(r = 100)),
                 technical_cor_plot + labs(subtitle = "c)"),
                 layout_matrix = cbind(c(1,1,1), c(2,3,3), c(2,3,3)))
  
  ggsave("figures/contig_stats/combined_technical_plot.png", plot = combined_technical_plot, height = 8, width = 10)
  ggsave("figures/contig_stats/combined_technical_plot.pdf", plot = combined_technical_plot, height = 8, width = 10)
}