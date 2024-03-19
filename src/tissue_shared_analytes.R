# Metabolomics Analysis - Venn Diagrams/Heatmap/Upset plots

source("src/load_packages.R")
library(glue)

df.interact <-
  read_excel("files/caltech.PD mice.022223.xlsx", sheet = "interaction.regcoef") %>% 
  janitor::clean_names() 

tissue_levels <- c(
  "Plasma",
  "Brainstem",
  "Cortex",
  "Nigra",
  "Striatum",
  "Duodenum",
  "Duodenum Content",
  "Cecum",
  "Colon",
  "Colon Content"
)

# -----------------------------------------------------------------------------
#             Microbiota/Genotype/Interaction Significant Analytes 
# -----------------------------------------------------------------------------

df.interact.A <-
  read_excel("files/caltech.PD mice.022223.xlsx", 
             sheet = "interaction.regcoef") %>% 
  janitor::clean_names() %>%
  dplyr::filter(imputation == "imputed") 

# -----------------------------------------------------------------------------
#                              Significance Heatmap 
# -----------------------------------------------------------------------------

df.sig.M <- 
  df.interact.A %>% filter(p_value <= 0.05, term == "MicrobiotaSPF") %>% 
  select(metabolite, tissue) %>% 
  pivot_wider(names_from = "tissue", 
              values_from = "tissue", 
              values_fill = "0",
              values_fn = function(x) "1")
df.heatmap.M <- df.sig.M %>% 
  pivot_longer(!metabolite, names_to = "tissue", values_to = "detected") %>% 
  mutate(term = "Microbiota Effect \n(SPF / GF)")

df.sig.G <- 
  df.interact.A %>% filter(p_value <= 0.05, term == "GenotypeASO") %>% 
  select(metabolite, tissue) %>% 
  pivot_wider(names_from = "tissue", 
              values_from = "tissue", 
              values_fill = "0",
              values_fn = function(x) "1")
df.heatmap.G <- df.sig.G %>% 
  pivot_longer(!metabolite, names_to = "tissue", values_to = "detected") %>% 
  mutate(term = "Genotype Effect \n(ASO / WT)")

df.sig.interact <- 
  df.interact.A %>% filter(p_value <= 0.05, term == "MicrobiotaSPF:GenotypeASO") %>% 
  select(metabolite, tissue) %>% 
  pivot_wider(names_from = "tissue", 
              values_from = "tissue", 
              values_fill = "0",
              values_fn = function(x) "1")
df.heatmap.I <- df.sig.interact %>% 
  pivot_longer(!metabolite, names_to = "tissue", values_to = "detected") %>% 
  mutate(term = "Microbiota x Genotype \nInteraction")

df.heatmap <- 
  bind_rows(df.heatmap.M, df.heatmap.G, df.heatmap.I) %>% 
  mutate(tissue = factor(tissue, levels = tissue_levels))

shared.sig.heatmap <- df.heatmap %>% 
  ggplot(aes(x=tissue, y=metabolite, fill=detected)) +
  geom_tile() +
  theme_classic() +
  labs(x = NULL, y = "Analytes") +
  scale_fill_manual(values = c("0" = "#f1f1f1", "1" = "#434343")) +
  facet_wrap(~term) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "none"
  )

ggsave(shared.sig.heatmap, filename = glue("data/Heatmaps/{Sys.Date()}_Shared_Significance_Heatmap_facet.svg"), 
       width = 6, height = 9)



# -----------------------------------------------------------------------------
#                              Upset Plot 
# -----------------------------------------------------------------------------

df.upset.micro <-
  df.heatmap %>%
  filter(term == "Microbiota Effect \n(SPF / GF)") %>%
  pivot_wider(names_from = 'tissue', values_from = 'detected') %>%
  column_to_rownames(var = "metabolite") %>% 
  select(-c(term)) %>% 
  mutate_if(is.character, as.numeric) %>% 
  as.data.frame()

svg(
  glue("data/Venn_Diagrams/{Sys.Date()}_Upset_Microbiota.svg"),
  width = 8,
  height = 4.5
)

upset(
  df.upset.micro,
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

# dev.copy(png,glue("data/Venn_Diagrams/{Sys.Date()}_Upset_Microbiota.png"), 
#          width = 600, height = 300)
dev.off()

df.upset.geno <-
  df.heatmap %>%
  filter(term == "Genotype Effect \n(ASO / WT)") %>%
  pivot_wider(names_from = 'tissue', values_from = 'detected') %>%
  column_to_rownames(var = "metabolite") %>% 
  select(-c(term)) %>% 
  mutate_if(is.character, as.numeric) %>% 
  as.data.frame()

svg(
  glue("data/Venn_Diagrams/{Sys.Date()}_Upset_Genotype.svg"),
  width = 8,
  height = 4.5
)

upset(
  df.upset.geno,
  nintersects =30,
  nsets = 10,
  order.by = "freq",
  decreasing = T,
  mb.ratio = c(0.6, 0.4),
  number.angles = 0,
  text.scale = 1.1,
  point.size = 2.8,
  line.size = 1
)

# dev.copy(png,glue("data/Venn_Diagrams/{Sys.Date()}_Upset_Genotype.png"),
#          width = 600, height = 300)
dev.off()



svg(
  glue("data/Venn_Diagrams/{Sys.Date()}_Upset_Microbiota_x_Genotype.svg"),
  width = 8,
  height = 4.5
)

df.upset.interaction <-
  df.heatmap %>%
  filter(term == "Microbiota x Genotype \nInteraction") %>%
  pivot_wider(names_from = 'tissue', values_from = 'detected') %>%
  column_to_rownames(var = "metabolite") %>% 
  select(-c(term)) %>% 
  mutate_if(is.character, as.numeric) %>% 
  as.data.frame()

upset(
  df.upset.interaction,
  nintersects =30,
  nsets = 10,
  order.by = "freq",
  decreasing = T,
  mb.ratio = c(0.6, 0.4),
  number.angles = 0,
  text.scale = 1.1,
  point.size = 2.8,
  line.size = 1
)

# dev.copy(png,glue("data/Venn_Diagrams/{Sys.Date()}_Upset_Microbiota_x_Genotype.png"), 
#          width = 600, height = 300)
dev.off()

write.csv(df.upset.micro, file ="files/Upset_Microbiota.csv", quote = FALSE)
write.csv(df.upset.geno, file ="files/Upset_Genotype.csv", quote = FALSE)
write.csv(df.upset.interaction, file ="files/Upset_Microbiota_x_Genotype.csv", quote = FALSE)


# -----------------------------------------------------------------------------
#                       VENN DIAGRAMS  -  Brain Tissue 
# -----------------------------------------------------------------------------


# Microbiota
venn.diagram(
  x = list(
    df.interact.A %>% filter(tissue == "Brainstem", p_value < 0.05, term == "MicrobiotaSPF") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Cortex", p_value < 0.05, term == "MicrobiotaSPF") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Nigra", p_value < 0.05, term == "MicrobiotaSPF") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Striatum", p_value < 0.05, term == "MicrobiotaSPF") %>% select(metabolite_name) %>% unlist()
  ),
  category.names = c("Brainstem", "Cortex" , "Nigra", "Striatum"),
  filename = glue("data/Venn_Diagrams/{Sys.Date()}_venn_brain_Microbiota.png"),
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
  margin = 0.1,
  ext.text = TRUE
)


# Genotype
venn.diagram(
  x = list(
    df.interact.A %>% filter(tissue == "Brainstem", p_value < 0.05, term == "GenotypeASO") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Cortex", p_value < 0.05, term == "GenotypeASO") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Nigra", p_value < 0.05, term == "GenotypeASO") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Striatum", p_value < 0.05, term == "GenotypeASO") %>% select(metabolite_name) %>% unlist()
  ),
  category.names = c("Brainstem", "Cortex" , "Nigra", "Striatum"),
  filename = glue("data/Venn_Diagrams/{Sys.Date()}_venn_brain_Genotype.png"),
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
  margin = 0.1,
  ext.text = TRUE
)



# Genotype
venn.diagram(
  x = list(
    df.interact.A %>% filter(tissue == "Brainstem", p_value < 0.05, term == "MicrobiotaSPF:GenotypeASO") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Cortex", p_value < 0.05, term == "MicrobiotaSPF:GenotypeASO") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Nigra", p_value < 0.05, term == "MicrobiotaSPF:GenotypeASO") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Striatum", p_value < 0.05, term == "MicrobiotaSPF:GenotypeASO") %>% select(metabolite_name) %>% unlist()
  ),
  category.names = c("Brainstem", "Cortex" , "Nigra", "Striatum"),
  filename = glue("data/Venn_Diagrams/{Sys.Date()}_venn_brain_Microbiota_X_Genotype.png"),
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
  margin = 0.1,
  ext.text = TRUE
)




# -----------------------------------------------------------------------------
#                              Intestinal Tissue 
# -----------------------------------------------------------------------------

# Microbiota
venn.diagram(
  x = list(
    df.interact.A %>% filter(tissue == "Cecum", p_value < 0.05, term == "MicrobiotaSPF") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Colon", p_value < 0.05, term == "MicrobiotaSPF") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Colon Content", p_value < 0.05, term == "MicrobiotaSPF") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Duodenum", p_value < 0.05, term == "MicrobiotaSPF") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Duodenum Content", p_value < 0.05, term == "MicrobiotaSPF") %>% select(metabolite_name) %>% unlist()
  ),
  category.names = c("Cecum", "Colon", "Colon Content", "Duodenum", "Duodenum Content"),
  filename = glue("data/Venn_Diagrams/{Sys.Date()}_venn_gut_Microbiota.png"),
  output = TRUE ,
  imagetype="png" ,
  height = 960 , 
  width = 960 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.4,
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
  margin = 0.3,
  ext.text = TRUE
)


# Genotype
venn.diagram(
  x = list(
    df.interact.A %>% filter(tissue == "Cecum", p_value < 0.05, term == "GenotypeASO") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Colon", p_value < 0.05, term == "GenotypeASO") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Colon Content", p_value < 0.05, term == "GenotypeASO") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Duodenum", p_value < 0.05, term == "GenotypeASO") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Duodenum Content", p_value < 0.05, term == "GenotypeASO") %>% select(metabolite_name) %>% unlist()
  ),
  category.names = c("Cecum", "Colon", "Colon Content", "Duodenum", "Duodenum Content"),
  filename = glue("data/Venn_Diagrams/{Sys.Date()}_venn_gut_Genotype.png"),
  output = TRUE ,
  imagetype="png" ,
  height = 960 , 
  width = 960 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.4,
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
  margin = 0.3,
  ext.text = TRUE
)



# Microbiota X Genotype
venn.diagram(
  x = list(
    df.interact.A %>% filter(tissue == "Cecum", p_value < 0.05, term == "MicrobiotaSPF:GenotypeASO") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Colon", p_value < 0.05, term == "MicrobiotaSPF:GenotypeASO") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Colon Content", p_value < 0.05, term == "MicrobiotaSPF:GenotypeASO") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Duodenum", p_value < 0.05, term == "MicrobiotaSPF:GenotypeASO") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Duodenum Content", p_value < 0.05, term == "MicrobiotaSPF:GenotypeASO") %>% select(metabolite_name) %>% unlist()
  ),
  category.names = c("Cecum", "Colon", "Colon Content", "Duodenum", "Duodenum Content"),
  filename = glue("data/Venn_Diagrams/{Sys.Date()}_venn_gut_Microbiota_X_Genotype.png"),
  output = TRUE ,
  imagetype="png" ,
  height = 960 , 
  width = 960 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.4,
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
  margin = 0.3,
  ext.text = TRUE
)



# 
# brainstem <-
#   read_excel("files/caltech.PD mice_G_by_M_results.xlsx", sheet = "Brainstem") %>% 
#   janitor::clean_names() %>% mutate(tissue = "Brainstem")
# cecum <-
#   read_excel("files/caltech.PD mice_G_by_M_results.xlsx", sheet = "Cecum") %>% 
#   janitor::clean_names() %>% mutate(tissue = "Cecum")
# colon <-
#   read_excel("files/caltech.PD mice_G_by_M_results.xlsx", sheet = "Colon") %>% 
#   janitor::clean_names() %>% mutate(tissue = "Colon")
# colonContent <-
#   read_excel("files/caltech.PD mice_G_by_M_results.xlsx", sheet = "Colon Content") %>% 
#   janitor::clean_names() %>% mutate(tissue = "Colon Content")
# cortex <-
#   read_excel("files/caltech.PD mice_G_by_M_results.xlsx", sheet = "Cortex") %>% 
#   janitor::clean_names() %>% mutate(tissue = "Cortex")
# duodenum <-
#   read_excel("files/caltech.PD mice_G_by_M_results.xlsx", sheet = "Duodenum") %>% 
#   janitor::clean_names() %>% mutate(tissue = "Duodenum")
# duodenumContent <-
#   read_excel("files/caltech.PD mice_G_by_M_results.xlsx", sheet = "Duodenum Content") %>% 
#   janitor::clean_names() %>% mutate(tissue = "Duodenum Content")
# nigra <-
#   read_excel("files/caltech.PD mice_G_by_M_results.xlsx", sheet = "Nigra") %>% 
#   janitor::clean_names() %>% mutate(tissue = "Nigra")
# plasma <-
#   read_excel("files/caltech.PD mice_G_by_M_results.xlsx", sheet = "Plasma") %>% 
#   janitor::clean_names() %>% mutate(tissue = "Plasma")
# striatum <-
#   read_excel("files//caltech.PD mice_G_by_M_results.xlsx", sheet = "Striatum") %>% 
#   janitor::clean_names() %>% mutate(tissue = "Striatum")
# 
# 
# df.GbyM.all <-
#   full_join(
#     dplyr::select(brainstem, c(metabolite, tissue)) %>% dplyr::rename(Brainstem = tissue),
#     dplyr::select(cecum, c(metabolite, tissue)) %>% dplyr::rename(Cecum = tissue), 
#     by="metabolite") %>% 
#   full_join(dplyr::select(colon, c(metabolite, tissue)) %>% dplyr::rename(Colon = tissue), by="metabolite") %>% 
#   full_join(dplyr::select(colonContent, c(metabolite, tissue)) %>% dplyr::rename(Colon_content = tissue), by="metabolite") %>% 
#   full_join(dplyr::select(cortex, c(metabolite, tissue)) %>% dplyr::rename(Cortex = tissue), by="metabolite") %>% 
#   full_join(dplyr::select(duodenum, c(metabolite, tissue)) %>% dplyr::rename(Duodenum = tissue), by="metabolite") %>% 
#   full_join(dplyr::select(duodenumContent, c(metabolite, tissue)) %>% dplyr::rename(Duodenum_content = tissue), by="metabolite") %>% 
#   full_join(dplyr::select(nigra, c(metabolite, tissue)) %>% dplyr::rename(Nigra = tissue), by="metabolite") %>% 
#   full_join(dplyr::select(plasma, c(metabolite, tissue)) %>% dplyr::rename(Plasma = tissue), by="metabolite") %>% 
#   full_join(dplyr::select(striatum, c(metabolite, tissue)) %>% dplyr::rename(Striatum = tissue), by="metabolite")


