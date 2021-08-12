
source("src/load_packages.R")
source("src/functions.R")

# DF prep
df.interact.A <-
  read_excel("files/caltech.PD mice.040820.xlsx", 
             sheet = "interaction.regcoef") %>% 
  janitor::clean_names() %>%
  dplyr::filter(imputation == "imputed") 

df.barplot <- 
  df.interact.A %>%
  group_by(tissue, variable, analyte_class) %>%
  dplyr::summarise(mean = mean(estimate), n = n()) %>%
  ungroup() %>%
  dplyr::mutate(analyte_class = as.factor(analyte_class)) %>% 
  dplyr::mutate(tissue = factor(
    tissue,
    levels = tisse_order
  ))
df.barplot.sig <- 
  df.interact.A %>%
  # PVALUE FILTER
  dplyr::filter(pr_t < 0.05) %>%
  group_by(tissue, variable, analyte_class) %>%
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

# Analyte Class color palette
analyte_class_colors <- 
  colorRampPalette(colormash)(length(base::unique(
    df.barplot.sig$analyte_class)))
names(analyte_class_colors) <- unique(
  df.barplot.sig$analyte_class)

save(analyte_class_colors, file = "files/analyte_class_colors.RData")


#-------------------------------------------------------------------------------
#                            Plotting all Estimates
#-------------------------------------------------------------------------------
estimate_all <- 
  df.barplot %>%
  filter(variable != "(Intercept)") %>%
  ggplot(aes(x = mean,
             y = fct_rev(tissue),
             fill = analyte_class)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  geom_bar(size = 10^-1, stat = "identity", 
           width = 0.6, color = "black") +
  theme_bw() +
  labs(y = NULL, x = "Analyte Class Mean Estimate",
       fill = "Class") +
  facet_wrap( ~ variable, labeller = labeller(variable = facelabs)) +
  scale_fill_manual(values = analyte_class_colors) +
  simple_theme()

ggsave(
  estimate_all,
  filename = "data/Analyte_class_estimates/Analyte_Class_Estimate_average_all.svg",
  width = 15,
  height = 4.5,
  dpi = 600
)

#-------------------------------------------------------------------------------
#                            Plotting all Estimates - free x-axis
#-------------------------------------------------------------------------------
estimate_all_free <- 
  df.barplot %>%
  filter(variable != "(Intercept)") %>%
  ggplot(aes(x = mean,
             y = fct_rev(tissue),
             fill = analyte_class)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  geom_bar(size = 10^-1, stat = "identity", width = 0.6, color = "black") +
  theme_bw() +
  labs(y = NULL, x = "Analyte Class Mean Estimate",
       fill = "Class") +
  facet_wrap( ~ variable, labeller = labeller(variable = facelabs), scales = "free_x") +
  scale_fill_manual(values = analyte_class_colors) +
  simple_theme()

ggsave(
  estimate_all_free,
  filename = "data/Analyte_class_estimates/Analyte_Class_Estimate_average_all_free.svg",
  width = 15,
  height = 4.5,
  dpi = 600
)


#-------------------------------------------------------------------------------
#                     Plotting Significant Estimates
#-------------------------------------------------------------------------------

estimate_sig <- 
  df.barplot.sig %>%
  filter(variable != "(Intercept)") %>%
  ggplot(aes(x = mean,
             y = fct_rev(tissue),
             fill = analyte_class)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  geom_bar(size = 10^-1, stat = "identity", width = 0.6, color = "black") +
  theme_bw() +
  labs(y = NULL, x = "Analyte Class Mean Estimate",
       fill = "Class") +
  facet_wrap( ~ variable, labeller = labeller(variable = facelabs)) +
  scale_fill_manual(values = analyte_class_colors) +
  simple_theme()

ggsave(
  estimate_sig,
  filename = "data/Analyte_class_estimates/Analyte_Class_Estimate_average_significant.svg",
  width = 15,
  height = 4.5,
  dpi = 600
)


#-------------------------------------------------------------------------------
#               Plotting Significant Estimates - free x-axis
#-------------------------------------------------------------------------------
estimate_sig_free <- 
  df.barplot.sig %>%
  filter(variable != "(Intercept)") %>%
  ggplot(aes(x = mean,
             y = fct_rev(tissue),
             fill = analyte_class)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  geom_bar(size = 10^-1, stat = "identity", width = 0.6, color = "black") +
  theme_bw() +
  labs(y = NULL, x = "Analyte Class Mean Estimate",
       fill = "Class") +
  facet_wrap( ~ variable, labeller = labeller(variable = facelabs), scales = "free_x") +
  scale_fill_manual(values = analyte_class_colors) +
  simple_theme()


ggsave(
  estimate_sig_free,
  filename = "data/Analyte_class_estimates/Analyte_Class_Estimate_average_significant_free.svg",
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
  filter(variable == "MicrobiotaSPF:GenotypeASO") %>%
  ggplot(aes(x = mean,
             y = fct_rev(tissue),
             fill = analyte_class)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  geom_bar(size = 10^-1, stat = "identity", width = 0.6, color = "black") +
  theme_bw() +
  labs(y = NULL, x = "Analyte Class Mean Estimate",
       fill = "Class") +
  facet_wrap( ~ variable, labeller = labeller(variable = facelabs), scales = "free_x") +
  scale_fill_manual(values = analyte_class_colors) +
  simple_theme()

ggsave(
  estimate_sig_free_interact,
  filename = "data/Analyte_class_estimates/Analyte_Class_Estimate_average_significant_free_interactionOnly.svg",
  width = 8,
  height = 5
  )



#-------------------------------------------------------------------------------
#                     Tissue Specific Class Enrichment
#-------------------------------------------------------------------------------


#----------------------------------
#        Plotting Loop
#----------------------------------

for(tissue_target in unique(df.barplot.sig$tissue)){
  
  cat("Tissue: ", tissue_target, "\n")
  # tissue_rlang <- rlang::sym(tissue)
  
  plot <- df.barplot.sig %>%
    filter(tissue == tissue_target,
           variable != "(Intercept)") %>%
    ggplot(aes(x=mean, y = fct_reorder(analyte_class, desc(mean)), fill = analyte_class)) +
    theme_bw() +
    geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
    facet_wrap(~variable, nrow = 1, scales = "free_x",
               labeller = labeller(variable = facelabs)) +
    geom_segment(aes(x=0, xend=mean,
                     y=fct_reorder(analyte_class, desc(mean)),
                     yend=fct_reorder(analyte_class, desc(mean))), color = "grey") +
    geom_point(aes(size = n), shape=21, stroke = 0.2) +
    labs(y = NULL, x = "Analyte Class Mean Estimate",
         fill = "Class") +
    scale_fill_manual(values = analyte_class_colors) +
    simple_theme()
  
  plot.legend <- cowplot::plot_grid(cowplot::get_legend(plot))
  plot <- plot + theme(legend.position = "none")
  print(plot)

  ggsave(plot, filename =
           paste0("data/Analyte_class_estimates/Subplot_", tissue_target ,"_Analyte_class_enrichment.svg"),
         width = 9, height = 4)
  ggsave(plot.legend, filename =
           paste0("data/Analyte_class_estimates/Subplot_", tissue_target ,"_Analyte_class_enrichment_LEGEND.svg"),
         width = 5, height = 10)
  
  }



#----------------------------------
#        Brain Tissue
#----------------------------------

brain.summary <- 
  df.barplot.sig %>%
  filter(tissue %in% c("Brainstem", "Cortex", 
                       "Striatum", "Nigra"),
         variable != "(Intercept)") %>%
  ggplot(aes(x=mean, y = fct_reorder(analyte_class, desc(mean)), fill = tissue)) +
  theme_bw() +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  facet_wrap(~variable, nrow = 1, scales = "free_x",
             labeller = labeller(variable = facelabs)) +
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
         paste0("data/Analyte_class_estimates/Subplot_BrainSummary_Analyte_class_enrichment.svg"),
       width = 10, height = 6)

#----------------------------------
#        Gut Tissue
#----------------------------------

gut.summary <- 
  df.barplot.sig %>%
  filter(tissue %in% c("Duodenum", "Duodenum Content", 
                       "Colon", "Colon Content", "Cecum"),
         variable != "(Intercept)") %>%
  ggplot(aes(x=mean, y = fct_reorder(analyte_class, desc(mean)), fill = tissue)) +
  theme_bw() +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  facet_wrap(~variable, nrow = 1, scales = "free_x",
             labeller = labeller(variable = facelabs)) +
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
         paste0("data/Analyte_class_estimates/Subplot_GutSummary_Analyte_class_enrichment.svg"),
       width = 10, height = 6)


