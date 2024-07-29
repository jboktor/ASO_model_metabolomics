
source("src/load_packages.R")
source("src/functions.R")
analyte_class_colors <- readRDS(glue("files/analyte_class_colors_grouped.rds"))
# load(file = "files/analyte_class_colors.RData") # analyte_class_colors

df.interact.A <-
  read_excel("files/caltech.PD mice.022223.xlsx", 
             sheet = "interaction.regcoef") %>% 
  janitor::clean_names() %>%
  dplyr::filter(imputation == "imputed") 

#' Metabolites of interest are altered between 
#' ASO and WT mice in SPF conditions, 
#' not altered in a germ-free context, 
#' and metabolites altered between GF and SPF animals in an ASO genotype, 
#' not altered in a WT background

df.genotype.conditional <-
  read_excel("files/caltech.PD mice.022223.xlsx",
             sheet = "conditional.genotype.effect") %>%
  janitor::clean_names() %>%
  dplyr::filter(imputation == "imputed")

df.microbiota.conditional <-
  read_excel("files/caltech.PD mice.022223.xlsx",
             sheet = "conditional.microbita.effect") %>%
  janitor::clean_names() %>%
  dplyr::filter(imputation == "imputed")



for (toi in unique(df.interact.A$tissue)){
  
  cat('Plotting: ' , toi, '\n')
  
  df.interact <- df.interact.A %>%
    dplyr::filter(tissue == toi) %>%
    pivot_wider(
      names_from = "term",
      names_sep = ".",
      values_from = c(
        "estimate",
        "std_error",
        "statistic",
        "p_value",
        "x2_5_percent",
        "x97_5_percent"
      )
    )
  
  genotype_sig_mets <- df.interact %>% 
    filter(p_value.GenotypeASO <= 0.05) %>% 
    pull(metabolite)
  microbiota_sig_mets <- df.interact %>% 
    filter(p_value.MicrobiotaSPF <= 0.05) %>% 
    pull(metabolite)
  
  conditional_genotype_hits <-
    df.genotype.conditional %>%
    filter(tissue == toi) %>%
    select(microbiota, metabolite, wilcox_p) %>%
    pivot_wider(names_from = 'microbiota',
                names_prefix = 'ASO_vs_WT_in_',
                values_from = 'wilcox_p') %>%
    filter(ASO_vs_WT_in_SPF <= 0.05 & ASO_vs_WT_in_GF > 0.05) %>% 
    filter(metabolite %in% genotype_sig_mets)
  
  conditional_microbiota_hits <-
    df.microbiota.conditional %>%
    filter(tissue == toi) %>%
    select(genotype, metabolite, wilcox_p) %>%
    pivot_wider(names_from = 'genotype',
                names_prefix = 'GF_vs_SPF_in_',
                values_from = 'wilcox_p') %>%
    filter(GF_vs_SPF_in_ASO <= 0.05 &
             GF_vs_SPF_in_WT > 0.05) %>% 
    filter(metabolite %in% microbiota_sig_mets)
  
  features_of_interest <- full_join(
    conditional_genotype_hits,
    conditional_microbiota_hits, 
    by = 'metabolite'
  )
  
  plot <- interact_plot(df.interact, features_of_interest) 
  print(plot)
  ggsave(paste0("figures/Interaction_plots/", Sys.Date(), "_Interaction_plot_", toi, "_v2.0.svg"), 
         #dpi = 600,
         width = 6, height = 6)
}


