source("src/load_packages.R")
source("src/functions.R")
analyte_class_colors <- readRDS(glue("files/analyte_class_colors_grouped.rds"))


# DF prep
df.interact.A <-
  read_excel("files/caltech.PD mice.022223.xlsx", 
             sheet = "interaction.regcoef") %>% 
  janitor::clean_names() %>%
  dplyr::filter(imputation == "imputed") 

# #_________________
# # Delete me 
# library(plotly)
# library(glue)
# 
# # df.interact.A$metabolite %>% unique
# # df.interact.A %>% 
# #   filter(metabolite == "TMAO")
# 
# p <- df.interact.A %>%
#   # filter(tissue == "Striatum", analyte_class == "Ceramides") %>%
#   filter(metabolite == "TMAO") %>%
#   ggplot(aes(x=fct_reorder(metabolite, estimate), y=estimate)) +
#   geom_point(aes(color = term)) +
#   theme_bw() +
#   facet_wrap(~tissue) +
#   theme(axis.text.x = element_text(angle=45))
# ggplotly(p)
# 
# df_met <- read.csv(
#   glue("files/caltechdata.csv"), header = F)
# 
# df_met %>% glimpse
# 
# #_______________

df.barplot <- 
  df.interact.A %>%
  filter(term != "BW") %>% 
  group_by(tissue, term, analyte_class) %>%
  dplyr::summarise(mean = mean(estimate), n = n()) %>%
  ungroup() %>%
  dplyr::mutate(analyte_class = factor(analyte_class, levels = names(analyte_class_colors))) %>% 
  dplyr::mutate(tissue = factor(
    tissue,
    levels = tisse_order
  ))
df.barplot.sig <- 
  df.interact.A %>%
  # PVALUE FILTER
  dplyr::filter(p_value < 0.05) %>%
  # REMOVE BW metrics
  dplyr::filter(term != "BW") %>%
  group_by(tissue, term, analyte_class) %>%
  dplyr::summarise(mean = mean(estimate), n = n()) %>%
  ungroup() %>%
  dplyr::mutate(analyte_class = as.factor(analyte_class)) %>% 
  dplyr::mutate(tissue = factor(
    tissue,
    levels = c(
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
  ))

# # Analyte Class color palette
# analyte_class_colors <- 
#   colorRampPalette(colormash)(length(base::unique(
#     df.barplot.sig$analyte_class)))
# names(analyte_class_colors) <- unique(
#   df.barplot.sig$analyte_class)
# save(analyte_class_colors, file = "files/analyte_class_colors.RData")




#-------------------------------------------------------------------------------
#                            Plotting all Estimates
#-------------------------------------------------------------------------------
estimate_all <- 
  df.barplot %>%
  filter(term != "(Intercept)") %>%
  ggplot(aes(x = mean,
             y = fct_rev(tissue),
             fill = analyte_class)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  geom_bar(size = 10^-1, stat = "identity", 
           width = 0.6, color = "black") +
  theme_bw() +
  labs(y = NULL, x = "Analyte Class Mean Estimate",
       fill = "Class") +
  facet_wrap( ~ term, labeller = labeller(term = facelabs)) +
  scale_fill_manual(values = analyte_class_colors) +
  simple_theme()

ggsave(
  estimate_all,
  filename = glue("data/Analyte_class_estimates/{Sys.Date()}_Analyte_Class_Estimate_average_all.svg"),
  width = 15,
  height = 4.5,
  dpi = 600
)

# -------------------------------------------------------------------------------
#                            Plotting all Estimates - free x-axis
#-------------------------------------------------------------------------------
# analyte_class_colors

estimate_all_free <-
  df.barplot %>%
  filter(term != "(Intercept)") %>%
  ggplot(aes(x = mean,
             y = fct_rev(tissue),
             fill = analyte_class)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  geom_bar(size = 10^-1, stat = "identity", width = 0.6, color = "black") +
  theme_bw() +
  labs(y = NULL, x = "Analyte Class Mean Estimate",
       fill = "Class") +
  facet_wrap( ~ term, labeller = labeller(term = facelabs), scales = "free_x") +
  scale_fill_manual(values = analyte_class_colors) +
  # scale_fill_manual(values = safecols) +
  # scale_fill_viridis_d(option = "H") +
  simple_theme()
colorblindr::cvd_grid(estimate_all_free)

ggsave(
  estimate_all_free,
  filename = glue("data/Analyte_class_estimates/{Sys.Date()}_Analyte_Class_Estimate_average_all_free.svg"),
  width = 15,
  height = 4.5,
  dpi = 600
)
# 
# library(RColorBrewer)
# n <- 25
# colrs <- brewer.pal.info[brewer.pal.info$colorblind == TRUE, ]
# col_vec = unlist(mapply(brewer.pal, colrs$maxcolors, rownames(colrs)))
# safecols <- sample(col_vec, n)
# # area <- rep(1,n)
# # pie(area, col = co
# 

#-------------------------------------------------------------------------------
#                     Plotting Significant Estimates
#-------------------------------------------------------------------------------

estimate_sig <- 
  df.barplot.sig %>%
  filter(term != "(Intercept)") %>%
  ggplot(aes(x = mean,
             y = fct_rev(tissue),
             fill = analyte_class)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  geom_bar(size = 10^-1, stat = "identity", width = 0.6, color = "black") +
  theme_bw() +
  labs(y = NULL, x = "Analyte Class Mean Estimate",
       fill = "Class") +
  facet_wrap( ~ term, labeller = labeller(term = facelabs)) +
  scale_fill_manual(values = analyte_class_colors) +
  simple_theme()

ggsave(
  estimate_sig,
  filename = glue("data/Analyte_class_estimates/{Sys.Date()}_Analyte_Class_Estimate_average_significant.svg"),
  width = 15,
  height = 4.5,
  dpi = 600
)


#-------------------------------------------------------------------------------
#               Plotting Significant Estimates - free x-axis
#-------------------------------------------------------------------------------
estimate_sig_free <- 
  df.barplot.sig %>%
  dplyr::mutate(analyte_class = factor(analyte_class, levels = names(analyte_class_colors))) %>% 
  filter(term != "(Intercept)") %>%
  ggplot(aes(x = mean,
             y = fct_rev(tissue),
             fill = analyte_class)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  geom_bar(size = 10^-1, stat = "identity", width = 0.6, color = "black") +
  theme_bw() +
  labs(y = NULL, x = "Analyte Class Mean Estimate",
       fill = "Class") +
  facet_wrap( ~ term, labeller = labeller(term = facelabs), scales = "free_x") +
  scale_fill_manual(values = analyte_class_colors) +
  simple_theme()
# estimate_sig_free

ggsave(
  estimate_sig_free,
  filename = glue("data/Analyte_class_estimates/{Sys.Date()}_Analyte_Class_Estimate_average_significant_free.svg"),
  width = 15,
  height = 4.5,
  dpi = 600
)
ggsave(
  estimate_sig_free,
  filename = glue("data/Analyte_class_estimates/{Sys.Date()}_Analyte_Class_Estimate_average_significant_free.png"),
  width = 15,
  height = 4.5,
  dpi = 600
)

#-------------------------------------------------------------------------------
#               Plotting Significant Estimates - free x-axis
#                 Microbiota X Genotype Interaction only
#-------------------------------------------------------------------------------

estimate_sig_free_interact <- 
  df.barplot.sig %>%
  filter(term == "MicrobiotaSPF:GenotypeASO") %>%
  ggplot(aes(x = mean,
             y = fct_rev(tissue),
             fill = analyte_class)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  geom_bar(size = 10^-1, stat = "identity", width = 0.6, color = "black") +
  theme_bw() +
  labs(y = NULL, x = "Analyte Class Mean Estimate",
       fill = "Class") +
  facet_wrap( ~ term, labeller = labeller(term = facelabs), scales = "free_x") +
  scale_fill_manual(values = analyte_class_colors) +
  simple_theme()

ggsave(
  estimate_sig_free_interact,
  filename = glue("data/Analyte_class_estimates/{Sys.Date()}_Analyte_Class_Estimate_average_significant_free_interactionOnly.svg"),
  width = 8,
  height = 5
  )



#-------------------------------------------------------------------------------
#                     Tissue Specific Class Enrichment
#-------------------------------------------------------------------------------


#----------------------------------
#        Plotting Loop
#----------------------------------

# for(tissue_target in unique(df.barplot.sig$tissue)){
for(tissue_target in "Plasma"){
    
  cat("Tissue: ", tissue_target, "\n")
  # tissue_rlang <- rlang::sym(tissue)
  
  plot <- df.barplot.sig %>%
    filter(tissue == tissue_target,
           term != "(Intercept)") %>%
    ggplot(aes(x=mean, y = fct_reorder(analyte_class, desc(mean)), fill = analyte_class)) +
    theme_bw() +
    geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
    facet_wrap(~term, nrow = 1, scales = "free_x",
               labeller = labeller(term = facelabs)) +
    geom_segment(aes(x=0, xend=mean,
                     y=fct_reorder(analyte_class, desc(mean)),
                     yend=fct_reorder(analyte_class, desc(mean))), color = "grey") +
    geom_point(aes(size = n), shape=21, stroke = 0.2) +
    labs(y = NULL, x = "Analyte Class Mean Estimate",
         fill = "Class") +
    scale_size_continuous(
      limits = c(0, 50),
      breaks = c(5, 10, 20, 30, 40)) +
    scale_fill_manual(values = analyte_class_colors) +
    simple_theme()
  
  plot.legend <- cowplot::plot_grid(cowplot::get_legend(plot))
  plot <- plot + theme(legend.position = "none")
  print(plot)

  ggsave(plot, filename =
           glue("data/Analyte_class_estimates/{Sys.Date()}_Subplot_{tissue_target}_Analyte_class_enrichment.svg"),
         width = 9, height = 4)
  ggsave(plot.legend, filename =
           glue("data/Analyte_class_estimates/{Sys.Date()}_Subplot_{tissue_target}_Analyte_class_enrichment_LEGEND.svg"),
         width = 5, height = 10)
  ggsave(plot, filename =
           glue("data/Analyte_class_estimates/{Sys.Date()}_Subplot_{tissue_target}_Analyte_class_enrichment.png"),
         width = 9, height = 4)
  ggsave(plot.legend, filename =
           glue("data/Analyte_class_estimates/{Sys.Date()}_Subplot_{tissue_target}_Analyte_class_enrichment_LEGEND.png"),
         width = 5, height = 10)
  }



#----------------------------------
#        Brain Tissue
#----------------------------------

brain.summary <- 
  df.barplot.sig %>%
  filter(tissue %in% c("Brainstem", "Cortex", 
                       "Striatum", "Nigra"),
         term != "(Intercept)") %>%
  ggplot(aes(x=mean, y = fct_reorder(analyte_class, desc(mean)), fill = tissue)) +
  theme_bw() +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  facet_wrap(~term, nrow = 1, scales = "free_x",
             labeller = labeller(term = facelabs)) +
  geom_segment(aes(x=0, xend=mean,
                   y=fct_reorder(analyte_class, desc(mean)),
                   yend=fct_reorder(analyte_class, desc(mean))), color = "grey") +
  geom_point(aes(size = n), shape=21, stroke = 0.2) +
  labs(y = NULL, x = "Analyte Class Mean Estimate",
       fill = "Tissue") +
  scale_fill_manual(values = tissue_cols) +
  simple_theme()

brain.summary <- brain.summary + theme(legend.position = "bottom")
print(brain.summary)

ggsave(brain.summary, filename =
         glue("data/Analyte_class_estimates/{Sys.Date()}_Subplot_BrainSummary_Analyte_class_enrichment.svg"),
       width = 10, height = 6)

#----------------------------------
#        Gut Tissue
#----------------------------------

gut.summary <- 
  df.barplot.sig %>%
  filter(tissue %in% c("Duodenum", "Duodenum Content", 
                       "Colon", "Colon Content", "Cecum"),
         term != "(Intercept)") %>%
  ggplot(aes(x=mean, y = fct_reorder(analyte_class, desc(mean)), fill = tissue)) +
  theme_bw() +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  facet_wrap(~term, nrow = 1, scales = "free_x",
             labeller = labeller(term = facelabs)) +
  geom_segment(aes(x=0, xend=mean,
                   y=fct_reorder(analyte_class, desc(mean)),
                   yend=fct_reorder(analyte_class, desc(mean))), color = "grey") +
  geom_point(aes(size = n), shape=21, stroke = 0.2) +
  labs(y = NULL, x = "Analyte Class Mean Estimate",
       fill = "Tissue") +
  scale_fill_manual(values = tissue_cols) +
  simple_theme()

gut.summary <- gut.summary + theme(legend.position = "bottom")
print(gut.summary)

ggsave(gut.summary, filename =
         glue("data/Analyte_class_estimates/{Sys.Date()}_Subplot_GutSummary_Analyte_class_enrichment.svg"),
       width = 10, height = 6)


