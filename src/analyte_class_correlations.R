# Correlation Analysis

# Create mean analyte class value for each mouse 

source("src/load_packages.R")
source("src/functions.R") 
load("files/analyte_class_colors.RData")
library(glue)

# load data ----
df.quant.raw <-
  read_excel("files/5370_Cal Tech PD Mouse Model_Gut brain axis study_09112020.xlsx", 
             sheet = "q500 data")
df.quant.LODhits <-
  read_excel("files/5370_Cal Tech PD Mouse Model_Gut brain axis study_09112020.xlsx", 
             sheet = "%<LOD")
df.quant.LODvals <-
  read_excel("files/5370_Cal Tech PD Mouse Model_Gut brain axis study_09112020.xlsx", 
             sheet = "LOD Thresholds")

df.analyte.map <-
  read_excel("files/caltech.PD mice.022223.xlsx", 
             sheet = "interaction.anova") %>% 
  janitor::clean_names() %>% 
  select(metabolite, metabolite_name, analyte_class) %>% 
  distinct()
df.analyte.map.slim <- select(df.analyte.map, metabolite, analyte_class)

sampleInfo <- df.quant.raw %>% 
  select(`Sample Bar Code`:Animal)
sampleInfo.list <- colnames(sampleInfo)[-1]
sampleInfo.list.all <- colnames(sampleInfo)

# Impute data ----


# calculate correlations
df.quant <- df.quant.raw # imputed data

df.quant.long <- df.quant.raw %>% # imputed data
  pivot_longer(!sampleInfo.list.all, names_to  ="metabolite") %>% 
  left_join(df.analyte.map, by = "metabolite")

# # summarize analyte class means per tissue
# df.analyte.summary <- df.quant.long %>% 
#   group_by(analyte_class, `Additional Information`) %>%
#   dplyr::summarize(analyte_average = mean(value, na.rm = TRUE)) %>% 
#   drop_na(analyte_class)
# 
# df.quant.plasma <- df.analyte.summary %>% 
#   dplyr::filter(`Additional Information` == "plasma")
# df.quant.brainstem <- df.analyte.summary %>% 
#   dplyr::filter(`Additional Information` == "Brainstem")
# df.quant.cortex <- df.analyte.summary %>% 
#   dplyr::filter(`Additional Information` == "Cortex")
# df.quant.Striatum <- df.analyte.summary %>% 
#   dplyr::filter(`Additional Information` == "Striatum")
# df.quant.SN <- df.analyte.summary %>% 
#   dplyr::filter(`Additional Information` == "Substantia Nigra")
# df.quant.colon <- df.analyte.summary %>% 
#   dplyr::filter(`Additional Information` == "colon")
# df.quant.duodenum <- df.analyte.summary %>% 
#   dplyr::filter(`Additional Information` == "duodenum")
# df.quant.cecum <- df.analyte.summary %>% 
#   dplyr::filter(`Additional Information` == "cecum")
# df.quant.colonCon <- df.analyte.summary %>% 
#   dplyr::filter(`Additional Information` == "colon content")
# df.quant.duodenumCon <- df.analyte.summary %>% 
#   dplyr::filter(`Additional Information` == "duodenum content")


# summarize analyte class means per tissue
df.analyte.summary <- df.quant.long %>% 
  group_by(analyte_class, `Additional Information`) %>%
  dplyr::summarize(analyte_average = mean(value, na.rm = TRUE)) %>% 
  drop_na(analyte_class)

zero_var_filt <- function(df){
  df[, sapply(df, function(x) length(unique(x)) > 1)]
}


df.quant.plasma <- df.quant %>% 
  dplyr::filter(`Additional Information` == "plasma") %>% 
  column_to_rownames(var = "Sample Bar Code") %>% 
  dplyr::select(-sampleInfo.list) %>% 
  zero_var_filt()
df.quant.brainstem <- df.quant %>% 
  dplyr::filter(`Additional Information` == "Brainstem") %>% 
  column_to_rownames(var = "Sample Bar Code") %>% 
  dplyr::select(-sampleInfo.list) %>% 
  zero_var_filt()
df.quant.cortex <- df.quant %>% 
  dplyr::filter(`Additional Information` == "Cortex") %>% 
  column_to_rownames(var = "Sample Bar Code") %>% 
  dplyr::select(-sampleInfo.list) %>% 
  zero_var_filt()
df.quant.Striatum <- df.quant %>% 
  dplyr::filter(`Additional Information` == "Striatum") %>% 
  column_to_rownames(var = "Sample Bar Code") %>% 
  dplyr::select(-sampleInfo.list) %>% 
  zero_var_filt()
df.quant.SN <- df.quant %>% 
  dplyr::filter(`Additional Information` == "Substantia Nigra") %>% 
  column_to_rownames(var = "Sample Bar Code") %>% 
  dplyr::select(-sampleInfo.list) %>% 
  zero_var_filt()
df.quant.colon <- df.quant %>% 
  dplyr::filter(`Additional Information` == "colon") %>% 
  column_to_rownames(var = "Sample Bar Code") %>% 
  dplyr::select(-sampleInfo.list) %>% 
  zero_var_filt()
df.quant.duodenum <- df.quant %>% 
  dplyr::filter(`Additional Information` == "duodenum") %>% 
  column_to_rownames(var = "Sample Bar Code") %>% 
  dplyr::select(-sampleInfo.list) %>% 
  zero_var_filt()
df.quant.cecum <- df.quant %>% 
  dplyr::filter(`Additional Information` == "cecum") %>% 
  column_to_rownames(var = "Sample Bar Code") %>% 
  dplyr::select(-sampleInfo.list) %>% 
  zero_var_filt()
df.quant.colonCon <- df.quant %>% 
  dplyr::filter(`Additional Information` == "colon content") %>% 
  column_to_rownames(var = "Sample Bar Code") %>% 
  dplyr::select(-sampleInfo.list) %>% 
  zero_var_filt()
df.quant.duodenumCon <- df.quant %>% 
  dplyr::filter(`Additional Information` == "duodenum content") %>% 
  column_to_rownames(var = "Sample Bar Code") %>% 
  dplyr::select(-sampleInfo.list) %>% 
  zero_var_filt()

df.cor.cecum <-
  corr.test(
    x = df.quant.plasma,
    y = df.quant.cecum,
    method = "spearman",
    adjust = "BH"
  )

df.cor.duodenum <-
  corr.test(
    x = df.quant.plasma,
    y = df.quant.duodenum,
    method = "spearman",
    adjust = "BH"
  )
df.cor.duodenumCon <-
  corr.test(
    x = df.quant.plasma,
    y = df.quant.duodenumCon,
    method = "spearman",
    adjust = "BH"
  )

df.cor.colon <-
  corr.test(
    x = df.quant.plasma,
    y = df.quant.colon,
    method = "spearman",
    adjust = "BH"
  )
df.cor.colonCon <-
  corr.test(
    x = df.quant.plasma,
    y = df.quant.colonCon,
    method = "spearman",
    adjust = "BH"
  )




df.cor.brainstem <-
  corr.test(
    x = df.quant.plasma,
    y = df.quant.brainstem,
    method = "spearman",
    adjust = "BH"
  )
df.cor.cortex <-
  corr.test(
    x = df.quant.plasma,
    y = df.quant.cortex,
    method = "spearman",
    adjust = "BH"
  )
df.cor.sn <-
  corr.test(
    x = df.quant.plasma,
    y = df.quant.SN,
    method = "spearman",
    adjust = "BH"
  )

# doParallel::registerDoParallel()
df.cor.str <-
  corr.test(
    x = df.quant.plasma,
    y = df.quant.Striatum,
    method = "spearman",
    adjust = "BH"
  )


# df.cor.str <-
#   cor.test(
#     x = df.quant.plasma,
#     y = df.quant.Striatum,
#     method = "spearman",
#     na.action = na.exclude,
#     alternative = "two.sided"
#   )

#_______________________________________________________________________________
#                              Plotting ----
# ______________________________________________________________________________


df.analyte.map.slimX <- dplyr::rename(df.analyte.map.slim,
                                      "tissue_X" = metabolite,
                                      "analyte_class_X" = analyte_class)
df.analyte.map.slimY <- dplyr::rename(df.analyte.map.slim,
                                      "tissue_Y" = metabolite,
                                      "analyte_class_Y" = analyte_class)

df.cor.cecum.plt <- 
  df.cor.cecum$r %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "tissue_X") %>% 
  pivot_longer(!tissue_X, names_to = "tissue_Y", values_to = "rho") %>% 
  dplyr::filter(tissue_X == tissue_Y) %>% 
  left_join(df.analyte.map.slimX, by = "tissue_X") %>% 
  left_join(df.analyte.map.slimY, by = "tissue_Y")

# df.cor.cecum.plt %>% 
#   ggplot(aes(x=analyte_class_X, y = rho)) +
#   geom_violin() +
#   geom_point(aes(fill = analyte_class_X), 
#              position = position_jitterdodge(jitter.width = 0.1), 
#              shape = 21, size = 0.5) +
#   scale_fill_manual(values = analyte_class_colors) +
#   my_clean_theme() +
#   theme(legend.position = "none")
  
tst <-
  df.cor.cecum.plt %>% 
  group_by(analyte_class_X) %>%
  dplyr::summarize(class_cor_average = mean(rho, na.rm = TRUE))



cross_tissue_class_corrs <- function(rho_df){
  
  df.cor <- 
    rho_df %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "tissue_X") %>% 
    pivot_longer(!tissue_X, names_to = "tissue_Y", values_to = "rho") %>% 
    dplyr::filter(tissue_X == tissue_Y) %>% 
    left_join(df.analyte.map.slimX, by = "tissue_X") %>% 
    left_join(df.analyte.map.slimY, by = "tissue_Y") %>% 
    group_by(analyte_class_X) %>%
    dplyr::summarize(class_cor_average = mean(rho, na.rm = TRUE))
  return(df.cor)
}

plasma.cecum.cors <-
  cross_tissue_class_corrs(df.cor.cecum$r) %>% mutate(tissue = "Cecum")
plasma.duodenum.cors <-
  cross_tissue_class_corrs(df.cor.duodenum$r) %>% mutate(tissue = "Duodenum")
plasma.duodenumCon.cors <-
  cross_tissue_class_corrs(df.cor.duodenumCon$r) %>% mutate(tissue = "Duodenum Content")


df.cor.final <- bind_rows(plasma.cecum.cors,
                          plasma.duodenum.cors,
                          plasma.duodenumCon.cors)


df.cor.final %>%
  ggplot(aes(x = tissue, y = class_cor_average, fill = analyte_class_X, group = analyte_class_X)) +
  geom_point(shape = 21) +
  scale_fill_manual(values = analyte_class_colors) +
  my_clean_theme() +
  geom_line(color = "gray") +
  theme(legend.position = "none")

