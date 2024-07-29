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



# Creating an upset plot visual for analytes that are correlated w/ plasma across multiple tissues

#' format a datafame with metabolites as rows and external tissue sources
#' as columns, values within the dataframe are 1 if there is a strong correlation
#' with the plasma and 0 if not

cor_mat <- strong_cors %>%
  select(featB, tissue_B) %>% 
  mutate(detected = 1) %>%
  pivot_wider(names_from = "tissue_B", 
              values_from = 'detected', values_fill = 0) %>% 
  column_to_rownames(var = "featB") %>% 
  # mutate(Plasma = 1) %>% 
  as.data.frame()


svg(
  glue("data/Venn_Diagrams/{Sys.Date()}_Upset_Plasma-Metabolite_corrs.svg"),
  width = 8,
  height = 4.5
)

upset(
  cor_mat,
  nintersects = 30,
  nsets = 10,
  order.by = "freq",
  decreasing = T,
  mb.ratio = c(0.6, 0.4),
  number.angles = 0,
  text.scale = 1.1,
  point.size = 2.8,
  line.size = 1
)

dev.off()


# Load necessary libraries
library(tidyverse)
library(tidygraph)
library(ggraph)

# Create nodes data
nodes <- strong_cors %>%
  select(featA, tissue_A) %>%
  dplyr::rename(name = featA, group = tissue_A) %>%
  distinct() %>%
  bind_rows(
    strong_cors %>%
      select(tissue_B) %>%
      dplyr::rename(name = tissue_B) %>%
      mutate(group = name) %>%
      distinct()
  ) %>%
  distinct() %>% 
  # Tissue nodes size 6, Metabolite nodes size 2
  mutate(type = ifelse(name %in% unique(edges$to), "tissue", "analyte")) %>% 
  glimpse()

# Create edges data
edges <- strong_cors %>%
  select(featA, tissue_B, value) %>%
  dplyr::rename(from = featA, to = tissue_B, correlation = value)

tissue_cols <- rep("lightgrey", length(unique(strong_cors$tissue_B)))
names(tissue_cols) <- unique(strong_cors$tissue_B)

graph_node_cols <- c(
  analyte_class_colors,
  tissue_cols
)

# Create graph object
graph <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)

network_plot <- ggraph(graph, layout = 'stress') +
  geom_edge_link(aes(edge_color = correlation), width = 1.5, alpha = 0.5, show.legend = TRUE) +
  geom_node_point(aes(size = type, shape = type, fill = type), color = "black") +
  geom_node_label(aes(label = name), repel = TRUE, alpha = 0.8) +
  theme_void() +
  scale_edge_color_gradient(low = "blue", high = "red") +
  scale_size_manual(values = c("tissue" = 8, "analyte" = 4)) +
  scale_fill_manual(values = c("tissue" = "white", "analyte" = "black")) +
  scale_shape_manual(values = c("tissue" = 21, "analyte" = 23)) +
  labs(edge_color = "Spearman's Rho",
       color = NULL) +
  theme(plot.title = element_text(hjust = 0.5))


ggsave(
  glue("data/Correlations/{Sys.Date()}_plasma_correlations_top_metabolite_network.svg"),
  network_plot, 
  width = 8, height  = 6
  )


# Legend: 
#' Network visual of metabolites selected in A). A metabolite is connected to a 
#' tissue node if it has a strong correlation with it's abundance in Plasma samples.
#' TMAO is highly connected, with strong correlations with Plasma levels in six different
#' specimen types.  

