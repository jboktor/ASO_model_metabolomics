
source("src/load_packages.R")
source("src/functions.R")

df.interact.A <-
  read_excel("files/caltech.PD mice.040820.xlsx", 
             sheet = "interaction.regcoef") %>% 
  janitor::clean_names() %>%
  dplyr::filter(imputation == "imputed") 

for (toi in unique(df.interact.A$tissue)){
  print(toi)
  df.interact <- df.interact.A %>%
    dplyr::filter(tissue == toi) %>%
    pivot_wider(
      names_from = "variable",
      names_sep = ".",
      values_from = c(
        "estimate",
        "std_error",
        "t_value",
        "pr_t",
        "x2_5_percent",
        "x97_5_percent"
      )
    )
  plot <- interact_plot(df.interact) 
  print(plot)
  ggsave(paste0("data/Interaction_plots/Interaction_plot_", toi, ".png"), 
         width = 6, height = 6)
}


