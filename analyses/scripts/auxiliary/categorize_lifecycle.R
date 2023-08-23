


# analyze -----------------------------------------------------------------

categorize_lifecycle <- function(ref_prop_provirus) {
  
  vOTU_lifecycle_cat <-
    vOTU_stats %>% 
    group_by(vOTU) %>% 
    mutate(prop_provirus = cluster_proviruses/cluster_contigs,
           pval_provirus_prop = case_when(cluster_contigs > 4 ~ binom.test(x = cluster_proviruses, n = cluster_contigs, p = ref_prop_provirus, alternative = "two.sided")$p.value)) %>% 
    ungroup() %>% 
    mutate(padj_provirus_prop = p.adjust(pval_provirus_prop, method = "fdr")) %>% 
    select(vOTU, starts_with("cluster"), prop_provirus, pval_provirus_prop, padj_provirus_prop) %>% 
    mutate(lifecycle_propensity_unadj = case_when(is.na(pval_provirus_prop) ~ "unknown",
                                                  pval_provirus_prop < 0.05 & prop_provirus == 0 ~ "exclusively lytic",
                                                  pval_provirus_prop < 0.05 & prop_provirus == 1 ~ "exclusively lysogenic",
                                                  pval_provirus_prop < 0.05 & prop_provirus < ref_prop_provirus ~ "predominantly lytic",
                                                  pval_provirus_prop < 0.05 & prop_provirus > ref_prop_provirus ~ "predominantly lysogenic",
                                                  pval_provirus_prop > 0.05 ~ "no difference"),
           lifecycle_propensity_adj = case_when(is.na(padj_provirus_prop) ~ "unknown",
                                                padj_provirus_prop < 0.05 & prop_provirus == 0 ~ "exclusively lytic",
                                                padj_provirus_prop < 0.05 & prop_provirus == 1 ~ "exclusively lysogenic",
                                                padj_provirus_prop < 0.05 & prop_provirus < ref_prop_provirus ~ "predominantly lytic",
                                                padj_provirus_prop < 0.05 & prop_provirus > ref_prop_provirus ~ "predominantly lysogenic",
                                                padj_provirus_prop > 0.05 ~ "no difference")) %>% 
    select(vOTU, starts_with("cluster"), prop_provirus, pval_provirus_prop, padj_provirus_prop, lifecycle_propensity_unadj, lifecycle_propensity_adj)
  
  vOTU_lifecycle_cat %>% 
    write_tsv("data/vOTU_stats/vOTU_lyfecycle_cat.tsv")
}
