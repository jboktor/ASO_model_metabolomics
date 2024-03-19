source("src/load_packages.R")
source("src/functions.R") 
# load("files/analyte_class_colors.RData")
analyte_class_colors <- readRDS(glue("files/analyte_class_colors_grouped.rds"))

corr_imputed_r <-
  read_excel("data/Correlations/caltech PD mice model Cor.060921.xlsx", 
             sheet = "imputed data - Spearman r")

df.analyte.map <-
  read_excel("files/caltech.PD mice.040820.xlsx", 
             sheet = "interaction.anova") %>% 
  janitor::clean_names() %>% 
  select(metabolite, metabolite_name, analyte_class) %>% 
  distinct()
df.analyte.map.slim <- 
  select(df.analyte.map, metabolite, analyte_class) %>% 
  dplyr::rename(featA = metabolite)

corr_plasma_df <- 
  corr_imputed_r %>% 
  dplyr::rename(featA = row.names) %>% 
  pivot_longer(-featA,  names_to = "featB") %>% 
  separate(featA, c("featA", "tissue_A"), sep = "___") %>% 
  separate(featB, c("featB", "tissue_B"), sep = "___") %>% 
  filter(tissue_A == "Plasma",
         featA == featB)
saveRDS(corr_plasma_df, "files/plasma_correlation_data.rds")

class_filter <- 
  corr_plasma_df %>% 
  left_join(df.analyte.map.slim, by = "featA") %>% 
  filter(tissue_B != "Plasma") %>% 
  drop_na(value) %>% 
  group_by(analyte_class, tissue_B) %>% 
  dplyr::summarise(count = n()) %>% 
  filter(count > 3) %>% 
  pull(analyte_class) %>% 
  unique()


# Tissue Specific X-axis
for (class in class_filter){
  
  plot <- corr_plasma_df %>% 
    left_join(df.analyte.map.slim, by = "featA") %>% 
    filter(tissue_B != "Plasma",
           analyte_class %in% class) %>% 
    dplyr::mutate(tissue_B = factor(tissue_B, levels = tissue_order)) %>% 
    dplyr::mutate(tissue_B_class = case_when(tissue_B %in% c("Brainstem", "Cortex", "Nigra", "Striatum") ~ "Brain",
                                             tissue_B %in% c("Duodenum", "Duodenum Content", "Cecum", "Colon", 
                                                             "Colon Content") ~ "Gut")) %>% 
    ggplot(aes(tissue_B, y = value)) + 
    geom_rect(aes(xmin=Inf, xmax=-Inf, ymin=0.25, ymax=Inf), fill = "#e9e9e9") +
    geom_rect(aes(xmin=Inf, xmax=-Inf, ymin=-0.25, ymax=-Inf), fill = "#e9e9e9") +
    geom_violin(size= 0.3, alpha = 0.9) +
    geom_point(aes(color = tissue_B_class), position = position_jitterdodge(jitter.width = 0.2), alpha = 0.4) +
    geom_boxplot(width = 0.15, outlier.alpha = 0, size = 0.3) +
    labs(y = "Spearman's Rho", color = NULL) +
    facet_wrap(~analyte_class) +
    scale_color_nejm() +
    my_clean_theme() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
  print(plot)
  
  plt_name <- paste0("data/Correlations/plasma_correlations_", class, ".svg")
  #gsave(plot, filename = plt_name, width = 6, height  = 4)
  
  }


# Brain Gut X-axis
for (class in class_filter){
  
  plot <- corr_plasma_df %>% 
    left_join(df.analyte.map.slim, by = "featA") %>% 
    filter(tissue_B != "Plasma",
           analyte_class %in% class) %>% 
    dplyr::mutate(tissue_B = factor(tissue_B, levels = tissue_order)) %>% 
    dplyr::mutate(tissue_B_class = case_when(tissue_B %in% c("Brainstem", "Cortex", "Nigra", "Striatum") ~ "Brain",
                                             tissue_B %in% c("Duodenum", "Duodenum Content", "Cecum", "Colon", 
                                                             "Colon Content") ~ "Gut")) %>% 
    ggplot(aes(tissue_B_class, y = value)) + 
    geom_rect(aes(xmin=Inf, xmax=-Inf, ymin=0.25, ymax=Inf), fill = "#e9e9e9") +
    geom_rect(aes(xmin=Inf, xmax=-Inf, ymin=-0.25, ymax=-Inf), fill = "#e9e9e9") +
    geom_violin(size= 0.3, alpha = 0.9) +
    geom_point(aes(color = tissue_B_class), position = position_jitterdodge(jitter.width = 0.3), alpha = 0.4) +
    geom_boxplot(width = 0.15, outlier.alpha = 0, size = 0.3, alpha = 0.9) +
    labs(y = "Spearman's Rho", color = NULL) +
    facet_wrap(~analyte_class) +
    scale_color_nejm() +
    my_clean_theme() +
     theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
  print(plot)
  
  plt_name <- paste0("data/Correlations/plasma_correlations_tissue_source_", class, ".svg")
  #gsave(plot, filename = plt_name, width = 3.2, height  = 4)
  
}


# Plot features w/ strongest correlations across all tissues 

corr_plasma_df <- readRDS("files/plasma_correlation_data.rds")

strong_cors <- 
  corr_plasma_df %>% 
  left_join(df.analyte.map.slim, by = "featA") %>% 
  filter(tissue_B != "Plasma") %>% 
  filter(value > 0.75 | value < -0.75) %>% 
  mutate(label_col = paste(featB,  tissue_B, sep = " "))

plasma_corr_top <- 
  strong_cors %>% 
  mutate(tissue_B=factor(tissue_B, levels=tissue_order)) %>%
  ggplot(aes(x=value, y= reorder(featA, value))) +
  geom_point(aes(fill=analyte_class), shape=21, size=3) +
  geom_segment(aes(x=0, xend=value, y=featA, yend=featA), color="gray") +
  theme_bw() +
  scale_y_discrete(position = "right") +
  facet_grid(rows = vars(tissue_B), scales = "free_y", space = "free", switch = "y") +
  labs(x="Spearman's Rho", title = "Metabolite Correlations to Plasma") +
  scale_fill_manual(values = analyte_class_colors, name ="Analyte Class") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        strip.text.y.left = element_text(angle = 0),
        panel.grid.major.y = element_blank())
plasma_corr_top
ggsave(plasma_corr_top, filename = glue("data/Correlations/{Sys.Date()}_plasma_correlations_top_metabolite_corrs.png"),
       width = 6, height  = 8)

  
  