

plot_contig_qual <- function() {
  
  # supplementary figure viral contigs quality ------------------------------
  
  checkV_res_all_by_sample <-
    checkv_res_all %>%
    group_by(sample_id) %>%
    summarize(avg_contig_len = mean(contig_length),
              putative_viral_scaffolds = n(),
              n_complete = sum(checkv_quality %in% "Complete"),
              n_high_qual = sum(checkv_quality %in% "High-quality"),
              n_med_qual = sum(checkv_quality %in% "Medium-quality"),
              n_low_qual = sum(checkv_quality %in% "Low-quality"),
              n_non_deter = sum(checkv_quality %in% "Not-determined")) %>%
    ungroup() %>%
    left_join(sample_data %>% select(sample_id, n_contigs) %>% mutate(sample_id = gsub("_", "-", sample_id)))
  
  
  ## S1:CheckV
  n_contigs_per_sample <-
    checkV_res_all_by_sample %>%
    pivot_longer(-c(sample_id, avg_contig_len, putative_viral_scaffolds, n_contigs), names_to = "checkv_quality", values_to = "n") %>% 
    mutate(checkv_quality = factor(checkv_quality, 
                                   levels = c("n_complete", "n_high_qual", "n_med_qual", "n_low_qual", "n_non_deter"),
                                   labels = c("Complete", "High-quality","Medium-quality","Low-quality", "Not-determined"))) %>%
    ggplot(aes(x = checkv_quality, y = n, fill = checkv_quality)) +
    scale_fill_manual(name = "Quality", values = c("#1c6462", "#1c6462","#1c6462", "red", "red"))+
    geom_boxplot() +
    theme_bw()+
    theme(legend.position = "none") +
    scale_y_continuous(trans = "log10")+
    labs(x= "", y= "viral genomes per sample")
  
  contig_length_per_quality <-
    checkv_res_all %>%
    mutate(checkv_quality = factor(checkv_quality, 
                                   levels = c("Complete", "High-quality","Medium-quality","Low-quality", "Not-determined"))) %>%
    ggplot(aes(x = checkv_quality, y = contig_length, fill = checkv_quality)) +
    scale_fill_manual(name = "Quality", values = c("#1c6462", "#1c6462","#1c6462", "red", "red"))+
    geom_boxplot() +
    theme_bw()+
    theme(legend.position = "none") +
    scale_y_continuous(trans = "log10")+
    labs(x= "", y= "viral genome length")
  
  contig_qual_plot <- 
    grid.arrange(n_contigs_per_sample + labs(subtitle = "a)"),
                 contig_length_per_quality + labs(subtitle = "b)"), ncol = 2)
  
  ggsave("figures/contig_stats/contig_quality.png", plot = contig_qual_plot, height = 3, width = 10)
  ggsave("figures/contig_stats/contig_quality.pdf", plot = contig_qual_plot, height = 3, width = 10)
}
