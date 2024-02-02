# 

source("src/load_packages.R")
source("src/functions.R")
library(magrittr)
library(glue)
library(ppcor)
library(plotly)

# load in table 
df_met <- read.csv("files/caltechdata.csv", header = TRUE)
df_met %<>% filter(imputation == "imputed") %>%
  mutate(group = glue("{Genotype}_{Microbiota}")) %>%
  mutate(
    Genotype_binary = case_when(Genotype == "ASO" ~ 1,
                                Genotype == "WT" ~ 0),
    Microbiota_binary = case_when(Microbiota == "SPF" ~ 1,
                                  Microbiota == "GF" ~ 0),
  )


tube_metadata <- c(
  "Sample.Bar.Code",
  "Sample.Identification",
  "Tube.label",
  "Material",
  "Species",
  "Plate.Number.LC",
  "Plate.Number.FIA",
  "Additional.Information",
  "Animal",
  "Microbiota",
  "Genotype",
  "BW",
  "Brainstem",
  "Nigra",
  "Striatum",
  "Cortex",
  "Duodenum",
  "Colon",
  "Colon_cont.",
  "Duo_cont.",
  "Ceacum",
  "Tissue",
  "imputation",
  "group",
  "Genotype_binary",
  "Microbiota_binary"
)


metabolites <- df_met %>% dplyr::select(-tube_metadata) %>% colnames()

calculate_body_weight_correlation <- function(df, metabolite) {
  suppressWarnings({
    stat_row_1 <- cor.test(df[['BW']], df[[metabolite]], method = "spearman") %>%
      broom::tidy() %>% 
      mutate(partial_correlation = FALSE) %>% 
      janitor::clean_names() %>% 
      dplyr::select(-c(alternative, method))
  })
  
  stat_row_2 <- pcor.test(
    x = df[['BW']],
    y = df[[metabolite]],
    z = list(df[['Microbiota_binary']],
             df[['Genotype_binary']]),
    method = "spearman"
  ) %>% 
    mutate(partial_correlation = TRUE) %>% 
    janitor::clean_names() %>% 
    dplyr::select(-c(n, gp, method))
  
  df_out <- bind_rows(stat_row_1, stat_row_2)
  return(df_out)
}


pb <-
  progress::progress_bar$new(total = length(metabolites) * length(df_met$Tissue))
corr_df <- tibble()
# For each tissue, take
for (t in unique(df_met$Tissue)) {
  df_tissue <- df_met %>% filter(Tissue == t)
  for (m in metabolites) {
    pb$tick()
    if (all(is.na(c(df_tissue[[m]]))) |
        length(unique(df_tissue[[m]])) < 2) {
      next
    }
    corr_df %<>% bind_rows(
      calculate_body_weight_correlation(df_tissue, m) %>%
        mutate(metabolite = m, tissue = t)
    )
  }
}

saveRDS(corr_df,
        glue(
          "{getwd()}/files/{Sys.Date()}_body-weight-correlations.rds"
        ))


corr_df_wide <- corr_df %>%
  mutate(partial_correlation = case_when(
    partial_correlation ~ "corrected",
    TRUE ~"uncorrected"
  )) %>% 
  pivot_wider(values_from = c("estimate", "statistic", "p_value"), 
              names_from = "partial_correlation")


p_pval <- corr_df_wide %>% 
  ggplot(aes(p_value_corrected, p_value_uncorrected)) +
  geom_point(alpha = 0.4, aes(group = metabolite)) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  facet_wrap(~tissue)

ggplotly(p_pval, tooltip = c("p_value_corrected", "p_value_uncorrected", "metabolite"))

p_rho <- corr_df_wide %>% 
  ggplot(aes(estimate_corrected, estimate_uncorrected)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(alpha = 0.4, aes(group = metabolite) ) +
  theme_bw() +
  coord_fixed() +
  facet_wrap(~tissue)

ggplotly(p_rho, tooltip = c("estimate_corrected", "estimate_uncorrected", "metabolite"))



corr_df_wide %>% 
  ggplot(aes(estimate_corrected, estimate_uncorrected)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(alpha = 0.4, aes(group = metabolite) ) +
  theme_bw() +
  coord_fixed() +
  facet_wrap(~tissue)


p_rhod <- corr_df_wide %>% 
  mutate(delta_rho = estimate_uncorrected - estimate_corrected) %>% 
  ggplot(aes(fct_reorder(metabolite, delta_rho), delta_rho)) +
  geom_point(aes(group = metabolite, color = p_value_corrected)) +
  theme_bw() +
  facet_grid(cols = vars(tissue), space = "free", scales = "free") +
  scale_colour_viridis_c(option = "G") +
  theme(axis.text.x = element_blank())

ggplotly(p_rhod, tooltip = c("estimate_corrected", "estimate_uncorrected", "metabolite"))





stats_df <- readRDS(glue("{wkdir}/files/2023-01-24_lm-comparisons-BW-BWquant.rds"))



# visualizing ranked metabolite significance and estimates using each model

plot_df <- stats_df %>% 
  mutate(term_category = case_when(
    grepl("BW", term) ~ "BodyWeight",
    TRUE ~ term
  )) %>% 
  filter(Tissue == "Striatum", 
         term != "(Intercept)")

metabolite_rank <- plot_df %>% 
  filter(term_category == "BodyWeight") %>% 
  group_by(term_category, metabolite) %>% 
  summarize(mean_est = mean(estimate)) %>% 
  arrange(mean_est) %>% 
  pull(metabolite)

bwp <- plot_df %>%
  mutate(metabolite = factor(metabolite, levels = metabolite_rank)) %>% 
  ggplot(aes(x = metabolite, y = estimate)) +
  geom_point(aes(shape = BW_data, color = p.value)) +
  facet_grid(rows = vars(term_category), cols = vars(Tissue), scales = "free") +
  scale_shape_manual(values = c("binned" = 16, "continuous" = 2)) +
  scale_color_viridis_c(option = "G") +
  theme_bw() +
  theme(axis.text.x = element_blank())

ggplotly(bwp, tooltip = c("estimate", "p.value", "metabolite"))


# estimate and p-value distributions
estimate_density <- plot_df %>% 
  ggplot(aes(estimate)) +
  geom_density(aes(color = BW_data)) +
  facet_grid(rows = vars(term_category), scales = "free") +
  theme_bw()
pval_density <- plot_df %>% 
  ggplot(aes(p.value)) +
  geom_density(aes(color = BW_data)) +
  facet_grid(rows = vars(term_category), scales = "free") +
  theme_bw()

density_plots <- (estimate_density + pval_density) + plot_layout(guides = "collect") 
  


estimate_plots <- list()
density_plots <- list()
for (t in tissues_list) {
  
  plot_df <- stats_df %>% 
    mutate(term_category = case_when(
      grepl("BW", term) ~ "BodyWeight",
      TRUE ~ term
    )) %>% 
    filter(Tissue == !!t, 
           term != "(Intercept)")
  
  metabolite_rank <- plot_df %>% 
    filter(term_category == "BodyWeight") %>% 
    group_by(term_category, metabolite) %>% 
    summarize(mean_est = mean(estimate)) %>% 
    arrange(mean_est) %>% 
    pull(metabolite)
  
  bwp <- plot_df %>%
    mutate(metabolite = factor(metabolite, levels = metabolite_rank)) %>% 
    ggplot(aes(x = metabolite, y = estimate)) +
    geom_point(aes(shape = BW_data, color = p.value)) +
    facet_grid(rows = vars(term_category), cols = vars(Tissue), scales = "free") +
    scale_shape_manual(values = c("binned" = 16, "continuous" = 2)) +
    scale_color_viridis_c(option = "G") +
    theme_bw() +
    theme(axis.text.x = element_blank())
  
  estimate_plots[[t]] <- bwp
  # ggplotly(bwp, tooltip = c("estimate", "p.value", "metabolite"))
  
  # estimate and p-value distributions
  estimate_density <- plot_df %>% 
    ggplot(aes(estimate)) +
    geom_density(aes(color = BW_data)) +
    facet_grid(rows = vars(term_category), scales = "free") +
    theme_bw()
  pval_density <- plot_df %>% 
    ggplot(aes(p.value)) +
    geom_density(aes(color = BW_data)) +
    facet_grid(rows = vars(term_category), scales = "free") +
    theme_bw()
  
  density_plots[[t]] <- (estimate_density + pval_density) + plot_layout(guides = "collect") 
}

estimate_plots$Brainstem
density_plots$Brainstem

# mouse_meta <-
#   readxl::read_xlsx(
#     "files/Metadata_Metabolomics 2019 _checked 2020_updated.xlsx",
#     sheet = 'Sheet2'
#   )
# mouse_meta %<>% select(Animal, group, body_weight) 
# df_met <- mouse_meta %<>% full_join(df_met, by = "Animal")