
source("src/load_packages.R")

facelabs<- c("Genotype Effect (WT / ASO)",
             "Microbiota Effect (SPF / GF)",
             "Microbiota x Genotype Interaction")
names(facelabs) <- c("GenotypeWT", 
                     "MicrobiotaSPF", 
                     "MicrobiotaSPF:GenotypeWT")

colormash <- c(
  # Set 1
  '#e41a1c', '#377eb8', '#4daf4a', '#984ea3',
  '#ff7f00', '#ffff33', '#a65628', '#f781bf','#999999',
  # Accent
  '#7fc97f','#beaed4','#fdc086','#ffff99',
  '#386cb0','#f0027f','#bf5b17','#666666',
  # Dark
  '#1b9e77','#d95f02','#7570b3','#e7298a',
  '#66a61e','#e6ab02','#a6761d','#666666')


# DF prep

df.interact.A <-
  read_excel("files/caltech.PD mice.011920.xlsx", 
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
  scale_fill_manual(values = colorRampPalette(colormash)(length(unique(
    df.barplot$analyte_class)))) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white"),
    plot.title = element_text(hjust = 0.5),
    plot.margin = unit(c(1, 1, 1, 2), "cm"))

ggsave(
  estimate_all,
  filename = "data/Analyte_class_estimates/Analyte_Class_Estimate_average_all.png",
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
  scale_fill_manual(values = colorRampPalette(colormash)(length(unique(
    df.barplot$analyte_class)))) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white"),
    plot.title = element_text(hjust = 0.5),
    plot.margin = unit(c(1, 1, 1, 2), "cm"))

ggsave(
  estimate_all_free,
  filename = "data/Analyte_class_estimates/Analyte_Class_Estimate_average_all_free.png",
  width = 15,
  height = 4.5,
  dpi = 600
)


#-------------------------------------------------------------------------------
#                     Plotting Significant Estimates
#-------------------------------------------------------------------------------
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
  scale_fill_manual(values = colorRampPalette(colormash)(length(unique(
    df.barplot.sig$analyte_class)))) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white"),
    plot.title = element_text(hjust = 0.5),
    plot.margin = unit(c(1, 1, 1, 2), "cm"))

ggsave(
  estimate_sig,
  filename = "data/Analyte_class_estimates/Analyte_Class_Estimate_average_significant.png",
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
  scale_fill_manual(values = colorRampPalette(colormash)(length(unique(
    df.barplot.sig$analyte_class)))) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white"),
    plot.title = element_text(hjust = 0.5),
    plot.margin = unit(c(1, 1, 1, 2), "cm"))

ggsave(
  estimate_sig_free,
  filename = "data/Analyte_class_estimates/Analyte_Class_Estimate_average_significant_free.png",
  width = 15,
  height = 4.5,
  dpi = 600
)




