# analysis functions


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

interact_plot <- function(df.interact){
  
  #' Creates XY plot of Microbiome/Genotype and labels top _25_ significant 
  #' interaction effect variables
  
  df.interact.nonsig <- df.interact %>% filter(`pr_t.MicrobiotaSPF:GenotypeWT` > 0.05)
  df.interact.sig <- df.interact %>% filter(`pr_t.MicrobiotaSPF:GenotypeWT` <= 0.05) %>% 
    dplyr::arrange(`pr_t.MicrobiotaSPF:GenotypeWT`, desc(abs(`estimate.MicrobiotaSPF:GenotypeWT`))) %>%
    slice_head(n = 25)
  
  interact.plot <- ggplot() +
    geom_point(data = df.interact.nonsig, 
               aes(x = estimate.MicrobiotaSPF, y = estimate.GenotypeWT), 
               alpha = 0.7, color = "grey") +
    geom_point(data = df.interact.sig,
               aes(x = estimate.MicrobiotaSPF, y = estimate.GenotypeWT, 
                   color = "red")) +
    theme_bw() +
    geom_abline(intercept = 0, slope = 1, linetype = 3, color = "darkgrey") +
    geom_abline(intercept = 0, slope = -1, linetype = 3, color = "darkgrey") +
    geom_vline(xintercept = 0, linetype = 1, color = "grey") +
    geom_hline(yintercept = 0, linetype = 1, color = "grey") +
    ggplot2::annotate("text", x=2.25, y=2.25, label= "SPF WT \n Enriched") + 
    ggplot2::annotate("text", x=2.25, y=-2.25, label= "SPF ASO \n Enriched") + 
    ggplot2::annotate("text", x=-2.25, y=2.25, label= "GF WT \n Enriched") + 
    ggplot2::annotate("text", x=-2.25, y=-2.25, label= "GF ASO \n Enriched") + 
    lims(x = c(-2.5, 2.5), y = c(-2.5, 2.5)) +
    labs(x = " Microbiota Effect Estimate (SPF / GF)",
         y = " Genotype Effect Estimate (WT / ASO)",
         title = df.interact.sig$tissue) +
    geom_label_repel(data = df.interact.sig,
                     aes(x = estimate.MicrobiotaSPF, y = estimate.GenotypeWT, label = metabolite), 
                     segment.alpha = 0.5,
                     segment.size = 0.2, 
                     size = 2.75,
                     force = 1,
                     max.time = 1,
                     max.iter = Inf,
                     max.overlaps = Inf) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "None")
  interact.plot + 
    return(interact.plot)
}

#--------------------------------------------------------------------------------


