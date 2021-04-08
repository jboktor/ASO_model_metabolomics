# analysis functions


#------------------------------------------------------------------------------
#  Functions
#------------------------------------------------------------------------------
GO_bubble_plot <- function(df){
  
  df %>% 
    top_n(15, wt=-q.value.FDR.B.H) %>% 
    mutate(hitsPerc=Hit.Count.in.Query.List*100/Hit.Count.in.Genome) %>%
    arrange(hitsPerc) %>%
    mutate(Name=factor(Name, levels=Name)) %>%
    ggplot(aes(x=hitsPerc, 
               y=Name, 
               colour=q.value.FDR.B.H, 
               size=Hit.Count.in.Query.List)) +
    geom_point() +
    theme_bw()+
    scale_color_viridis_c(option = "cividis") +
    expand_limits(x=1) +
    labs(x="Hits (%)", y="", colour="FDR", size="Count")
}

#------------------------------------------------------------------------------

GO_barplot <- function(df, fill){
  
  df %>% 
    drop_na("q.value.FDR.B.H") %>%
    mutate(Hit.Count.in.Query.List = as.integer(Hit.Count.in.Query.List)) %>%
    dplyr::mutate(hitsPerc=Hit.Count.in.Query.List*100/Hit.Count.in.Genome) %>%
    dplyr::arrange(q.value.FDR.B.H, desc(hitsPerc)) %>%
    slice_head(n = 15) %>% 
    dplyr::mutate(yvar = paste(Name, ID)) %>% 
    ggplot(aes(x= -log10(q.value.FDR.B.H), 
               y=reorder(yvar, -q.value.FDR.B.H))) +
    geom_bar(stat = "identity", fill = fill, width = 0.5) + 
    labs(x = expression(paste(log[10], " (Adjusted P-value)")), 
         y = "") +
    theme_classic() +
    scale_x_continuous(expand = c(0, 0)) 
}

#------------------------------------------------------------------------------


xysummary <- function(df, genekeydf, siglist,
                      x, y, xgroup, ygroup, title, fill, force){
  
  cols <- c("nodiff" = "#808080", "diff" = fill)
  col.rims <- c("nodiff" = "#ffffff", "diff" = fill)
  
  # remerge genename
  # select df of only significant genes
  sigfilter <- siglist %>% 
    dplyr::arrange(q_value, desc(abs(log2.fold_change.))) %>%
    slice_head(n = 25)
  
  cat("Top 20 most significant genes: \n\n")
  print(sigfilter$gene)
  
  df %>% 
    dplyr::mutate(genelabels = if_else(gene %in% sigfilter$gene, gene, "")) %>% 
    dplyr::mutate(goi_col = if_else(gene %in% siglist$gene, "diff", "nodiff")) %>% 
    dplyr::mutate(log2_value_1 = log2(value_1 + 1)) %>% 
    dplyr::mutate(log2_value_2 = log2(value_2 + 1)) %>%
    dplyr::arrange(desc(goi_col)) %>% 
    ggplot(aes(x=log2_value_1, y=log2_value_2)) +
    geom_point(aes(fill = goi_col, color = goi_col),
               shape=21, size=0.7, alpha = 0.7) +
    geom_text_repel(aes(label = genelabels), segment.alpha = 0.2, segment.size = 0.2, size =1.75, force = force, max.iter=6000) +
    theme_classic() +
    geom_abline(intercept = 0, slope = 1) +
    labs(x = paste0("expression (log2 FPKM + 1): ", xgroup),
         y = paste0("expression (log2 FPKM + 1): ", ygroup),
         title = title) +
    scale_color_manual(values = col.rims, name ="Group") +
    scale_fill_manual(values = cols, name ="Group") +
    theme(legend.position = "none", 
          plot.title = element_text(hjust = 0.5))
}



#------------------------------------------------------------------------------
boxplot_transformed <- function(goi, cols = group.cols, 
                    title = blank.title, ylabel = blank.ylabel){
  
  blank.title = " "; blank.ylabel = " "
  group.cols = c("ASO_GF"= "#FF6000", "ASO_SPF" = "#077E97", 
                 "WT_GF" = "#A00000", "WT_SPF" = "808080")
  
  set.seed(123)
  
  df.long %>% 
    filter(Genes == goi) %>% 
    ggplot(aes(x=group, y= log2(count + 1))) +
    geom_boxplot(aes(fill = group), alpha = 0.75, outlier.alpha = 0, width = 0.9) +
    geom_point(aes(fill = group), position = position_jitterdodge(jitter.width = 0.5), 
               shape=21, size=1.5, alpha = 1) +
    theme_classic() +
    ggtitle(goi) +
    labs(y = expression(paste(log[2], " (FPKM + 1)"))) +
    scale_color_manual(values = cols, name ="Group") +
    scale_fill_manual(values = cols, name ="Group") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(),
          panel.grid.major.y = element_blank(),
          legend.position = "none")
  
}

#------------------------------------------------------------------------------
boxplot_fpkm <- function(goi, cols = group.cols, 
                    title = blank.title, ylabel = blank.ylabel){
  
  blank.title = " "; blank.ylabel = " "
  group.cols = c("ASO_GF"= "#FF6000", "ASO_SPF" = "#077E97", 
                 "WT_GF" = "#A00000", "WT_SPF" = "808080")
  
  set.seed(123)
  
  df.long %>% 
    filter(Genes == goi) %>% 
    ggplot(aes(x=group, y= count)) +
    geom_boxplot(aes(fill = group), alpha = 0.75, outlier.alpha = 0, width = 0.9) +
    geom_point(aes(fill = group), position = position_jitterdodge(jitter.width = 0.5), 
               shape=21, size=1.5, alpha = 1) +
    theme_classic() +
    ggtitle(goi) +
    labs(y = "FPKM") +
    scale_color_manual(values = cols, name ="Group") +
    scale_fill_manual(values = cols, name ="Group") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(),
          panel.grid.major.y = element_blank(),
          legend.position = "none")
  
}

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



