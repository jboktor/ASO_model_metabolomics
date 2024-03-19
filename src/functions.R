# analysis functions


#------------------------------------------------------------------------------
#                       Aesthetic Variables
#------------------------------------------------------------------------------
tissue_order <-
  c("Plasma",
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

simple_theme <- function() {
  th <- ggplot2::theme_bw() +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "white"),
      plot.title = element_text(hjust = 0.5),
      plot.margin = unit(c(1, 1, 1, 2), "cm")
    )
  return(th)
}

facelabs<- c("Genotype Effect (ASO / WT)",
             "Microbiota Effect (SPF / GF)",
             "Microbiota x Genotype Interaction")
names(facelabs) <- c("GenotypeASO", 
                     "MicrobiotaSPF", 
                     "MicrobiotaSPF:GenotypeASO")


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


# Tissue color palette
tissue_cols <- ggsci::pal_npg()(10)
names(tissue_cols) <- tissue_order

# Group Color palette
groups_cols <-
  c("ASO GF" = "#a6cee3",
    "ASO SPF" = "#1f78b4",
    "WT GF" = "#b2df8a",
    "WT SPF" =  "#33a02c")



#------------------------------------------------------------------------------
#  Functions
#------------------------------------------------------------------------------

'%ni%' <- Negate('%in%')

#------------------------------------------------------------------------------

color_loop_generator <- function(names_col){
  # Init
  accent_pal <- 
    c('#7fc97f','#beaed4','#fdc086','#ffff99',
      '#386cb0','#f0027f','#bf5b17','#666666')
  color_output <- c()
  pal_length <- length(accent_pal)
  feats <- unique(names_col)
  n_cols <- length(feats)
  color_output <- c(rep(accent_pal, floor(n_cols/pal_length)), 
                    accent_pal[1:(n_cols %% pal_length)])
  names(color_output) <- feats
  
  return(color_output)
}


#--------------------------------------------------------------------------------
#                  Interaction Effect Plotting Function 
#--------------------------------------------------------------------------------

interact_plot <- function(df.interact, features_of_interest){
  
  #' Creates XY plot of Microbiome/Genotype and labels top _25_ significant 
  #' interaction effect variables
  
  df.interact.nonsig <- df.interact %>% filter(`p_value.MicrobiotaSPF:GenotypeASO` > 0.05)
  df.interact.sig <- df.interact %>% 
    filter(metabolite %in% features_of_interest$metabolite)
  df.interact.sig2name <- df.interact.sig %>% 
    # filter(metabolite %in% features_of_interest$metabolite) %>% 
    # filter(`pr_t.MicrobiotaSPF:GenotypeASO` <= 0.05) %>%
    dplyr::arrange(`p_value.MicrobiotaSPF:GenotypeASO`) %>%
    slice_head(n = 10)
  
  interact.plot <- ggplot() +
    geom_point(data = df.interact.nonsig, 
               aes(x = estimate.MicrobiotaSPF, y = estimate.GenotypeASO), 
               alpha = 0.5, fill = "grey", shape = 21, stroke = 0.1, size = 2) +
    geom_point(data = df.interact.sig,
               aes(x = estimate.MicrobiotaSPF, y = estimate.GenotypeASO, fill = analyte_class), 
               shape = 21, stroke = 0.1, size = 2 #color = "#ffa400"
               ) +
    theme_bw() +
    geom_abline(intercept = 0, slope = 1, linetype = 3, color = "darkgrey") +
    geom_abline(intercept = 0, slope = -1, linetype = 3, color = "darkgrey") +
    geom_vline(xintercept = 0, linetype = 1, color = "grey") +
    geom_hline(yintercept = 0, linetype = 1, color = "grey") +
    ggplot2::annotate("text", x=2.25, y=-2.25, label= "SPF WT \n Enriched") + 
    ggplot2::annotate("text", x=2.25, y=2.25, label= "SPF ASO \n Enriched") + 
    ggplot2::annotate("text", x=-2.25, y=-2.25, label= "GF WT \n Enriched") + 
    ggplot2::annotate("text", x=-2.25, y=2.25, label= "GF ASO \n Enriched") + 
    scale_fill_manual(values = analyte_class_colors) +
    lims(x = c(-2.5, 2.5), y = c(-2.5, 2.5)) +
    labs(x = " Microbiota Effect Estimate (SPF / GF)",
         y = " Genotype Effect Estimate (ASO / WT)",
         title = unique(df.interact$tissue)) +
    geom_label_repel(data = df.interact.sig2name,
                     aes(x = estimate.MicrobiotaSPF, y = estimate.GenotypeASO, label = metabolite), 
                     segment.alpha = 0.5,
                     segment.size = 0.2, 
                     alpha = 0.7,
                     size = 3.5,
                     force = 50,
                     max.time = 1,
                     max.iter = Inf,
                     max.overlaps = Inf) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "None")
    return(interact.plot)
}


#--------------------------------------------------------------------------------

my_clean_theme <- function() {
  th <- ggplot2::theme_bw() +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "white"),
      plot.title = element_text(hjust = 0.5),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    )
  return(th)
}
#--------------------------------------------------------------------------------


# Creating a new grouped color palette for analyte class -----
analyte_class_regrouped <- list(
  "G1" = c("Fatty Acids", "Diglycerides", "Triglycerides", "Acylcarnitines"),
  "G2" = c("Amino Acids", "Amino Acid Related", "Amine Oxides", "Biogenic Amines"),
  "G3" = c( "Ceramides", "Dihydroceramides", "Hexosylceramides", "Dihexosylceramides", "Trihexosylceramides", "Sphingomyelins"),
  "G4" = c("Phosphatidylcholines", "Lysophosphatidylcholines"),
  "G5" = c("Bile Acids", "Cholesteryl Esters"),
  "G6" = c("Indoles and Derivatives", "Nucleobases and Related", "Cresols", "Vitamins and Cofactors", "Alkaloids", "Carboxylic Acids", "Carbohydrates and Related")
)

col_pal_list <- c(
  "Blues",
  "Greens",
  "Purples",
  "Greys",
  "PuRd",
  "Reds"
)

i <- 1
regrouped_analyte_class_colpal <- list()
for (class in names(analyte_class_regrouped)){
  class_n <- analyte_class_regrouped[class] %>% unlist() %>% length()
  # if (class_n < 3){
  #   newcols <- brewer.pal(3, col_pal_list[i])[-1] # remove the top color
  # } else {
  #   newcols <- brewer.pal(class_n, col_pal_list[i])
  # }
  newcols <- brewer.pal(class_n + 2, col_pal_list[i])[-1][-1] # remove the top color 
  regrouped_analyte_class_colpal <- c(regrouped_analyte_class_colpal, newcols)
  i %<>% + 1
}

regrouped_analyte_class_colpal %<>% unlist()
names(regrouped_analyte_class_colpal) <- unname(unlist(analyte_class_regrouped))
saveRDS(regrouped_analyte_class_colpal,
        glue("files/analyte_class_colors_grouped.rds"))



