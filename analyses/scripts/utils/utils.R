

## colors
four_cat_col <- c("#885DB2", "#B25E5D", "#87B25D", "#5DB2B2")
three_cat_col <- c("#885DB2", "#B2885D", "#5DB288")

## project theme colors
col_pal_distinct <- c(dark_purple = "#6b519d", 
                      light_purple = "#9791c6", 
                      turquoise = "#5fc0bf",
                      greenblue = "#52a1b5",
                      green = "#1c6462",
                      gray = "#CEC9D6")

## Function to create gradient from gray to specified color
get_theme_gradient <- function(low_col = "gray", high_col = "dark_purple", n_cols = 4) {
  scales::seq_gradient_pal(low = col_pal_distinct[low_col], high = col_pal_distinct[high_col], space = "Lab")(seq(0,1,length.out = n_cols))  
}

## Function to make a tibble that can be used to create a PCoA plot
dist_to_PCoA <- function(distance_matrix, group_var = "group", PCoA_dims = c(1,2)) {
  pcoa_obj <- betadisper(d = distance_matrix, group = group_var)
  
  pcoa_obj$vectors %>% 
    data.frame() %>% 
    select(PCoA_dims) %>% 
    rownames_to_column("sample_id") %>%
    tibble() %>% 
    bind_cols(group = pcoa_obj$group) %>% 
    pivot_longer(-c(sample_id, group), names_to = "PCoA", values_to = "PCo_value") %>%
    left_join(pcoa_obj$centroids %>% 
                data.frame() %>% 
                select(PCoA_dims) %>% 
                rownames_to_column("group") %>% 
                tibble() %>% 
                pivot_longer(-c(group), names_to = "PCoA", values_to = "centroid")) %>% 
    left_join(pcoa_obj$eig %>% 
                enframe(name = "PCoA") %>% 
                mutate(rel_eig = value/sum(value)) %>% 
                select(PCoA, rel_eig)) %>% 
    pivot_longer(c(PCo_value, centroid), names_to = "var_type") #%>% 
}

## Function to plot two PCoA dimensions
plot_pcoa <- function(x, dim_1 = "PCoA1", dim_2 = "PCoA2") {
  tmp_labs <- x %>% 
    filter(PCoA %in% c(dim_1, dim_2)) %>% 
    mutate(lab = paste(PCoA, " (", round(rel_eig*100, digits = 1), "%)", sep = "")) %>% 
    group_by(lab) %>% 
    slice(1) %>% 
    ungroup() %>% 
    pull(lab)
  
  x <- x %>% 
    select(-rel_eig) %>% 
    pivot_wider(names_from = PCoA, values_from = value) %>%
    rename(dim1 = any_of(dim_1),
           dim2 = any_of(dim_2))
  x %>% 
    ggplot(aes(x = dim1, y = dim2, color = group, group = sample_id)) +
    stat_ellipse(data = x %>% filter(var_type %in% "PCo_value"), aes(x = dim1, y = dim2, color = group, group = group), level = .75, show.legend = FALSE) +
    geom_line(alpha = .1, show.legend = FALSE) +
    geom_point(data = x %>% filter(var_type %in% "centroid"), aes(x = dim1, y = dim2, fill = group), color = "black", shape = 21, size = 3) +
    coord_fixed() +
    theme_bw() +
    labs(color = "", fill = "", x = dim_1, y = dim_2)  +
    labs(x = tmp_labs[1],
         y = tmp_labs[2])
}

## micEco function to obtain omega^2 - effect size - for permanova
source("analyses/scripts/utils/adonis_Omega.R")


## Function to categorize continuous variables into tertiles - except in cases where > 33% of cases are 0, when tertiles are based on the positive cases.
cat_func <- function(x) {
  if(is.numeric(x)) {
    print("cat")
    ## median for those with more than 0 - if these comprise more than 33%
    if (sum(x %in% 0) >= 0.33*sum(!is.na(x))) {
      cut(x, 
          breaks = c(-Inf, 0, quantile(x[ x > 0], probs = seq(0, 1, length.out = 3), na.rm = TRUE, names = FALSE)[2], Inf), 
          labels = c("negative", "mid", "positive"))  
    } else {
      ## Otherwise, tertiles for whole range
      cut(x, 
          breaks = c(-Inf, quantile(x, probs = seq(0, 1, length.out = 4), na.rm = TRUE, names = FALSE)[2:3], Inf), 
          labels = c("negative", "mid", "positive"))  
    }
  } else {
    x
  }
}

## Plot effect sizes by variable
plot_pve <- function(pve_tab, effect_size_var = "parOmegaSq", colors = four_cat_col, sumSq = T) {
  if (!sumSq) {
    pve_tab <- pve_tab %>% mutate(df = "")
  }
  pve_tab <-
    pve_tab %>% 
    mutate(significance = case_when(p  == 0.001 ~ "***",
                                    p < 0.01  ~ "**",
                                    p < 0.05 ~ "*",
                                    TRUE ~ "")) %>% 
    mutate(var_lab = case_when(sumSq ~ paste0(var_name, " (df=", df, ", n=", n, ")"),
                               !sumSq ~ paste0(var_name, " (n=", n, ")"))) %>% 
    rename(eff_size = which(names(.) %in% effect_size_var))
  pve_tab %>% 
    ggplot(aes(x = eff_size, y = fct_reorder(var_lab, eff_size), fill = dataset, label = significance)) +
    geom_col() +
    geom_text(hjust = -0.5, vjust = .8) +
    scale_fill_manual(values = colors, name = "") +
    theme_bw() +
    scale_x_continuous(limits = c(NA, max(pve_tab$eff_size)+0.1*max(pve_tab$eff_size))) +
    labs(x = paste0("effect size (", expression("\u03a9"),"²)"),
         y = "")
}


## Create taxa summary
create_taxa_summary <- function(tax_dat, min_fam_freq = 20) {
  tmp_taxa_summary <-
    tax_dat %>% 
    mutate(Classification = case_when(Family %in% "n.a." ~ "Unknown",
                                      Family %in% c("Unclassified") ~ "Higher order",
                                      TRUE ~ "Classified")) %>% 
    group_by(Classification) %>% 
    summarize(n = n(),
              n_genomes = sum(cluster_contigs)) %>% 
    ungroup() %>% 
    # count(Classification) %>% 
    mutate(rel_n = n/sum(n),
           level="Overall") %>% 
    bind_rows(tax_dat %>% 
                filter(!Family %in% c("n.a.", "Unclassified")) %>% 
                mutate(Classification = Family) %>% 
                group_by(Classification) %>% 
                summarize(n = n(),
                          n_genomes = sum(cluster_contigs)) %>% 
                ungroup() %>% 
                # count(Classification) %>% 
                mutate(rel_n = n/sum(n),
                       level="Classified")) %>% 
    mutate(Classification = case_when(n < min_fam_freq ~ "Other",
                                      TRUE ~ Classification)) %>% 
    group_by(Classification, level) %>% 
    summarize(n = sum(n),
              n_genomes = sum(n_genomes),
              rel_n = sum(rel_n)) %>% 
    ungroup()
  
  tmp_taxa_summary %>% 
    mutate(level = factor(level, levels = c("Overall", "Classified"))) %>%
    mutate(Classification = factor(Classification, levels = c("Other",
                                                              (tmp_taxa_summary %>% 
                                                                 filter(level %in% "Classified", !Classification %in% c("Other")) %>% 
                                                                 arrange(n) %>% 
                                                                 pull(Classification)),
                                                              "Classified", "Higher order", "Unknown")))   
}


create_virus_taxonomy_lifecycle_table <- function(vir_abund, vir_stats, vir_lifecycle, min_n_for_summary = 20) {
  tmp_table <- 
    vir_abund %>% 
    pivot_longer(-sample_id, names_to = "vOTU", values_to = "abundance") %>% 
    group_by(vOTU) %>% 
    summarize(mean_present = mean(abundance[abundance > 0]),
              prevalence = sum(abundance>0)/n()) %>% 
    left_join(vir_stats %>% select(vOTU, Family, cluster_contigs)) %>% 
    left_join(vir_lifecycle %>% select(vOTU, lifecycle_propensity_adj)) %>% 
    group_by(Family) %>% 
    mutate(fam_count = n()) %>% 
    ungroup() %>% 
    mutate(Family = case_when(fam_count < min_n_for_summary ~ "Other",
                              Family %in% "n.a." ~ "Unknown",
                              Family %in% "Unclassified" ~ "Higher order",
                              TRUE ~ Family)) 
  
  fam_counts <- tmp_table %>% 
    count(Family) %>% 
    arrange(desc(n))
  
  fam_lvls <- pull(fam_counts, Family)
  for (i in c("Other", "Higher order", "Unknown")) {
    if (i %in% fam_lvls) {
      fam_lvls <- c(fam_lvls[-which(fam_lvls %in% i)], i)
    }
  }
  tmp_table %>% 
    mutate(Family = factor(Family, levels = fam_lvls))
}

## Order column names of matrix/dataframe by hierarchical clustering
order_by_cluster_columns <- function(x) {
  
  tmp_0 <- names(x)[which(apply(x == 0, 2, all))]
  x <- x[,!apply(x == 0, 2, all)]
  
  tmp <- x %>% scale() %>% t() %>% dist() %>% hclust()
  c(names(x)[tmp$order], tmp_0) 
  
}

## Order column names of matrix/dataframe by hierarchical clustering
order_by_cluster_variable <- function(x, var_tobe_clustered, var_to_cluster_on, value_to_cluster_on) {
  
  x %>% 
    select(all_of(c(var_tobe_clustered, var_to_cluster_on, value_to_cluster_on))) %>% 
    rename(x_var = all_of(var_tobe_clustered),
           y_var = all_of(var_to_cluster_on),
           n_var = all_of(value_to_cluster_on)) %>% 
    pivot_wider(names_from = x_var, values_from = n_var, values_fill = 0) %>% 
    as.data.frame() %>% 
    column_to_rownames("y_var") %>% 
    order_by_cluster_columns()
}
