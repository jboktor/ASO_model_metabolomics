# Metabolomics Analysis - Venn Diagrams/Heatmap/Upset plots

source("src/load_packages.R")

df.interact <-
  read_excel("files/caltech.PD mice.040820.xlsx", sheet = "interaction.regcoef") %>% 
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
  read_excel("files/caltech.PD mice.040820.xlsx", 
             sheet = "interaction.regcoef") %>% 
  janitor::clean_names() %>%
  dplyr::filter(imputation == "imputed") 

# -----------------------------------------------------------------------------
#                              Significance Heatmap 
# -----------------------------------------------------------------------------

df.sig.M <- 
  df.interact.A %>% filter(pr_t <= 0.05, variable == "MicrobiotaSPF") %>% 
  select(metabolite, tissue) %>% 
  pivot_wider(names_from = "tissue", 
              values_from = "tissue", 
              values_fill = "0",
              values_fn = function(x) "1")
df.heatmap.M <- df.sig.M %>% 
  pivot_longer(!metabolite, names_to = "tissue", values_to = "detected") %>% 
  mutate(variable = "Microbiota Effect \n(SPF / GF)")

df.sig.G <- 
  df.interact.A %>% filter(pr_t <= 0.05, variable == "GenotypeASO") %>% 
  select(metabolite, tissue) %>% 
  pivot_wider(names_from = "tissue", 
              values_from = "tissue", 
              values_fill = "0",
              values_fn = function(x) "1")
df.heatmap.G <- df.sig.G %>% 
  pivot_longer(!metabolite, names_to = "tissue", values_to = "detected") %>% 
  mutate(variable = "Genotype Effect \n(ASO / WT)")

df.sig.interact <- 
  df.interact.A %>% filter(pr_t <= 0.05, variable == "MicrobiotaSPF:GenotypeASO") %>% 
  select(metabolite, tissue) %>% 
  pivot_wider(names_from = "tissue", 
              values_from = "tissue", 
              values_fill = "0",
              values_fn = function(x) "1")
df.heatmap.I <- df.sig.interact %>% 
  pivot_longer(!metabolite, names_to = "tissue", values_to = "detected") %>% 
  mutate(variable = "Microbiota x Genotype \nInteraction")

df.heatmap <- 
  bind_rows(df.heatmap.M, df.heatmap.G, df.heatmap.I) %>% 
  mutate(tissue = factor(tissue, levels = tissue_levels))

shared.sig.heatmap <- df.heatmap %>% 
  ggplot(aes(x=tissue, y=metabolite, fill=detected)) +
  geom_tile() +
  theme_classic() +
  labs(x = NULL, y = "Analytes") +
  scale_fill_manual(values = c("0" = "#f1f1f1", "1" = "#434343")) +
  facet_wrap(~variable) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "none"
  )

ggsave(shared.sig.heatmap, filename = "data/Heatmaps/Shared_Significance_Heatmap_facet.png", 
       width = 6, height = 9, dpi = 2400)


# -----------------------------------------------------------------------------
#                              Upset Plot 
# -----------------------------------------------------------------------------


df.upset.micro <-
  df.heatmap %>%
  filter(variable == "Microbiota Effect \n(SPF / GF)") %>%
  pivot_wider(names_from = 'tissue', values_from = 'detected') %>%
  select(-c(metabolite, variable)) %>% 
  mutate_if(is.character, as.numeric) %>% 
  as.data.frame()

upset.micro <-
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

dev.copy(png,'data/Venn_Diagrams/Upset_Microbiota.png', 
         width = 600, height = 300)
dev.off()

df.upset.geno <-
  df.heatmap %>%
  filter(variable == "Genotype Effect \n(ASO / WT)") %>%
  pivot_wider(names_from = 'tissue', values_from = 'detected') %>%
  select(-c(metabolite, variable)) %>% 
  mutate_if(is.character, as.numeric) %>% 
  as.data.frame()

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

dev.copy(png,'data/Venn_Diagrams/Upset_Genotype.png', 
         width = 600, height = 300)
dev.off()

df.upset.interaction <-
  df.heatmap %>%
  filter(variable == "Microbiota x Genotype \nInteraction") %>%
  pivot_wider(names_from = 'tissue', values_from = 'detected') %>%
  select(-c(metabolite, variable)) %>% 
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

dev.copy(png,'data/Venn_Diagrams/Upset_Microbiota_x_Genotype.png', 
         width = 600, height = 300)
dev.off()

# df.upset <-
#   df.heatmap %>%
#   filter(variable == "Microbiota x Genotype \nInteraction") %>%
#   pivot_wider(names_from = 'tissue', values_from = 'detected') %>%
#   # select(-c(metabolite, variable)) %>% 
#   mutate_at(c("Brainstem", "Cecum", "Colon", "Colon Content", "Cortex", "Duodenum",
#               "Duodenum Content", "Nigra", "Plasma", "Striatum"), as.numeric) %>% 
#   as.data.frame()





# -----------------------------------------------------------------------------
#                       VENN DIAGRAMS  -  Brain Tissue 
# -----------------------------------------------------------------------------


# Microbiota
venn.diagram(
  x = list(
    df.interact.A %>% filter(tissue == "Brainstem", pr_t < 0.05, variable == "MicrobiotaSPF") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Cortex", pr_t < 0.05, variable == "MicrobiotaSPF") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Nigra", pr_t < 0.05, variable == "MicrobiotaSPF") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Striatum", pr_t < 0.05, variable == "MicrobiotaSPF") %>% select(metabolite_name) %>% unlist()
  ),
  category.names = c("Brainstem", "Cortex" , "Nigra", "Striatum"),
  filename = 'data/Venn_Diagrams/venn_brain_Microbiota.png',
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
    df.interact.A %>% filter(tissue == "Brainstem", pr_t < 0.05, variable == "GenotypeASO") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Cortex", pr_t < 0.05, variable == "GenotypeASO") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Nigra", pr_t < 0.05, variable == "GenotypeASO") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Striatum", pr_t < 0.05, variable == "GenotypeASO") %>% select(metabolite_name) %>% unlist()
  ),
  category.names = c("Brainstem", "Cortex" , "Nigra", "Striatum"),
  filename = 'data/Venn_Diagrams/venn_brain_Genotype.png',
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
    df.interact.A %>% filter(tissue == "Brainstem", pr_t < 0.05, variable == "MicrobiotaSPF:GenotypeASO") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Cortex", pr_t < 0.05, variable == "MicrobiotaSPF:GenotypeASO") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Nigra", pr_t < 0.05, variable == "MicrobiotaSPF:GenotypeASO") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Striatum", pr_t < 0.05, variable == "MicrobiotaSPF:GenotypeASO") %>% select(metabolite_name) %>% unlist()
  ),
  category.names = c("Brainstem", "Cortex" , "Nigra", "Striatum"),
  filename = 'data/Venn_Diagrams/venn_brain_Microbiota_X_Genotype.png',
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
    df.interact.A %>% filter(tissue == "Cecum", pr_t < 0.05, variable == "MicrobiotaSPF") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Colon", pr_t < 0.05, variable == "MicrobiotaSPF") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Colon Content", pr_t < 0.05, variable == "MicrobiotaSPF") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Duodenum", pr_t < 0.05, variable == "MicrobiotaSPF") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Duodenum Content", pr_t < 0.05, variable == "MicrobiotaSPF") %>% select(metabolite_name) %>% unlist()
  ),
  category.names = c("Cecum", "Colon", "Colon Content", "Duodenum", "Duodenum Content"),
  filename = 'data/Venn_Diagrams/venn_gut_Microbiota.png',
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
    df.interact.A %>% filter(tissue == "Cecum", pr_t < 0.05, variable == "GenotypeASO") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Colon", pr_t < 0.05, variable == "GenotypeASO") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Colon Content", pr_t < 0.05, variable == "GenotypeASO") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Duodenum", pr_t < 0.05, variable == "GenotypeASO") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Duodenum Content", pr_t < 0.05, variable == "GenotypeASO") %>% select(metabolite_name) %>% unlist()
  ),
  category.names = c("Cecum", "Colon", "Colon Content", "Duodenum", "Duodenum Content"),
  filename = 'data/Venn_Diagrams/venn_gut_Genotype.png',
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
    df.interact.A %>% filter(tissue == "Cecum", pr_t < 0.05, variable == "MicrobiotaSPF:GenotypeASO") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Colon", pr_t < 0.05, variable == "MicrobiotaSPF:GenotypeASO") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Colon Content", pr_t < 0.05, variable == "MicrobiotaSPF:GenotypeASO") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Duodenum", pr_t < 0.05, variable == "MicrobiotaSPF:GenotypeASO") %>% select(metabolite_name) %>% unlist(),
    df.interact.A %>% filter(tissue == "Duodenum Content", pr_t < 0.05, variable == "MicrobiotaSPF:GenotypeASO") %>% select(metabolite_name) %>% unlist()
  ),
  category.names = c("Cecum", "Colon", "Colon Content", "Duodenum", "Duodenum Content"),
  filename = 'data/Venn_Diagrams/venn_gut_Microbiota_X_Genotype.png',
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


