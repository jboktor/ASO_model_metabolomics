#####################################################################################
### Script for preprocessing Cal Tech PD Mouse Model - Q500
# Author: Sia Dehkordi - September 2020     
# Notes: Data was normlized for batch effect by Biocrates
# QC Steps:
#     1. Read in metabolomics data
#     2. Exclude Mets with >40% <LOD
#     3. impute LOD/2
#####################################################################################

library(openxlsx)
library(plyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(Hmisc)
library(ggrepel)
library(reshape2)
library(MASS)
library(M3C)
library(stringr)
library(data.table)
library(tsne)
library(coin)
library(ggforce)


setwd("/users/sm667/Desktop/Duke/Cal Tech PD Mice/")
rm(list=ls())
set.seed(20)   #set random seed for reproducibility purposed

source("./scripts/impute.knn.obs.sel.R")

set.seed(20)   #set random seed for reproducibility purposed

#####################################################################################
###### Read-in metabolomics data, LOD/ save them in one excell spreadsheet
#####################################################################################
mets.metadata <- read.xlsx("./upload to box/5370_Cal Tech PD Mouse Model_Gut brain axis study_09112020.xlsx","Q500 Metabolites Info")
rownames(mets.metadata) <- mets.metadata$Abbreviation

dta <-read.xlsx("./upload to box/5370_Cal Tech PD Mouse Model_Gut brain axis study_09112020.xlsx","q500 data")



meta <- read.xlsx("./upload to box/Metadata_Metabolomics 2019 _checked 2020_updated.xlsx")
meta <- meta[1:24, c(1,3:14)]
colnames(meta) <- meta[1,]
meta = meta[-1,]
meta [-c(1:3)] <- apply(meta [-c(1:3)] , 2, as.numeric)

LOD <- read.xlsx("./upload to box/5370_Cal Tech PD Mouse Model_Gut brain axis study_09112020.xlsx","LOD Thresholds") 
LOD[-c(1:2)] <- apply(LOD[-c(1:2)] , 2, as.numeric)
#LLOD <-read.xlsx("./upload to box/5370_Cal Tech PD Mouse Model_Gut brain axis study_09112020.xlsx","%<LOD")

colnames(dta) <- gsub("\\."," ",colnames(dta))
colnames(LOD) <- gsub("\\."," ",colnames(LOD))
#  colnames(LLOD) <- gsub("\\."," ",colnames(LLOD))
mets <- colnames(dta)[10:ncol(dta)] 

dta <- left_join(dta , meta , by = "Animal")

dta$Material = stringr::str_to_title(dta$Material)
dta$Tissue = str_extract(dta$`Sample Identification` , "[^_[0-9]+$]*") %>% str_to_title
#LLOD$Material = stringr::str_to_title(LLOD$Material)

#table( dta$Material)
#brain tissue             Gut content intestinal tract tissue                  plasma 
#88                      69                      46                            23   

# rr = dlply(dta ,.(Tissue) , function(x){
#   ddd = ddply(x , .(Microbiota), function(xx){
#     percent.b.LOD <- apply(xx %>% dplyr::select(mets), 2, function(xx){
#       mean(!is.na(xx))
#     })
#     
#   })
#   
#   return(ddd)
# })

#Include metabolites that have >%30 in SPF
mets.included = dlply(dta ,.(Material) , function(x){
  percent.b.LOD <- apply(x %>% subset(Microbiota=="SPF") %>% dplyr::select(mets), 2, function(xx){
    mean(!is.na(xx))
  })
  return (names(percent.b.LOD)[which(percent.b.LOD>0.3)])
})
lapply(mets.included, length)
#$`Brain Tissue`  : 297     Gut Content: 503 Intestinal Tract Tissue=459 Plasma =539

mets.metadata[mets.included$Plasma,]

mets.metadata[mets.included$Plasma,] %>% dplyr::select(Analyte.Class) %>% table(., useNA = "always")
#Acylcarnitines                 Alkaloids              Amine Oxides        Amino Acid Related               Amino Acids 
#17                         1                         1                        23                        20 
#Bile Acids           Biogenic Amines Carbohydrates and Related          Carboxylic Acids                 Ceramides 
#6                         7                         1                         5                        22 
#Cholesteryl Esters                   Cresols              Diglycerides        Dihexosylceramides          Dihydroceramides 
#22                         1                        24                         9                         3 
#Fatty Acids          Hexosylceramides   Indoles and Derivatives  Lysophosphatidylcholines   Nucleobases and Related 
#7                        18                         3                        13                         2 
#Phosphatidylcholines            Sphingomyelins             Triglycerides       Trihexosylceramides    Vitamins and Cofactors 
#74                        15                       239                         5                         1 


#Find metabolites that survived <LOD
common_mets <-  Reduce( intersect , mets.included )



dta$Genotype <- factor(dta$Genotype , levels = c("WT","ASO"))
dta$Microbiota <- factor(dta$Microbiota , levels = c("GF","SPF"))

#####################################################################################
#####################################################################################
#################impute data using LOD /2
extracols <- setdiff(colnames(dta) , mets)
dta.imputed <- left_join(dta , 
                         LOD %>% dplyr::select(-`Plate Number FIA`),
                         by = "Plate Number LC",
                         suffix = c("", ".LOD"))
##impute missing values with LOD/2
for ( met in mets){
  idx <- which(is.na(dta.imputed[,met]))
  dta.imputed[idx,met] <- dta.imputed[idx,paste0(met,".LOD")] / 2
}
dta.imputed <- dta.imputed %>% dplyr::select(-contains("LOD"))
###End of Imputation
####################################################################################
####################################################################################




anal.dta <- rbind(dta %>% mutate(imputation = "unimputed"),
                  dta.imputed %>%mutate(imputation = "imputed"))

####################################################################################
##########################tSNE Analysis
####################################################################################
  set.seed(40)  
  pdf(file="./figures/August 2021/tsne_caltech_jan2021_lod2Imputed_August 2021.pdf" ,onefile=T  , height = 6, width = 6.5)
  #Run tSNE on all samples on metabolites survived in all tissues
  tsne.res <-  tsne(dta.imputed[,common_mets] %>% log2  %>% scale , perplexity=50)
  tsne.res <- data.frame(tsne.res)
  colnames(tsne.res) <- c("t-SNE1","t-SNE2")
  tsne.res$Material <- dta$Material
  tsne.res$Genotype <- dta$Genotype
  tsne.res$Microbiota <- dta$Microbiota
  tsne.res$Tissue <- dta$Tissue
  tsne.res$G_by_M <- paste(tsne.res$Genotype , tsne.res$Microbiota)
  
  
  geno_cols <-c("ASO" = "#1f78b4", "WT" = "#33a02c")
  microb_cols <-c("GF" = "#1f78b4","SPF" = "#33a02c")
  G_by_M_cols <-c("ASO GF" = "#a6cee3","ASO SPF" = "#1f78b4","WT GF" = "#b2df8a","WT SPF" =  "#33a02c")
  material_order <-
    c("Plasma",
      "Brain Tissue",
      "Intestinal Tract Tissue",
      "Gut Content"
    )
  # Tissue color palette
  material_cols <- ggsci::pal_npg()(10)
  names(material_cols) <- material_order
  
  tisse_order <-c("Plasma","Brainstem","Cortex","Nigra","Striatum","Duodenum","Duodenum Content","Cecum","Colon","Colon Content")
  # Tissue color palette
  tissue_cols <- ggsci::pal_npg()(10)
  names(tissue_cols) <- tisse_order
  
  gg <- ggplot(tsne.res, aes(x=`t-SNE1`, y=`t-SNE2`,col=Material)) + 
    geom_point( size=4)+ggtitle("All Samples") +
    scale_color_manual(values =material_cols , aesthetics = c("colour", "fill"))+
    scale_x_continuous(limits = c(-18, 18))+
    scale_y_continuous(limits = c(-18, 18))+
    theme_bw()+
    theme( axis.text.x = element_text( face = "bold"),
           strip.text = element_text(face = "bold"),
           plot.title = element_text(hjust = 0.5, size=20),
           legend.position="bottom")  
  grid.arrange(gg)
  
  gg <- ggplot(tsne.res, aes(x=`t-SNE1`, y=`t-SNE2`,col=Tissue)) + 
    geom_point( size=4)+ggtitle("All Samples") +
    scale_color_manual(values =material_cols , aesthetics = c("colour", "fill"))+
    scale_x_continuous(limits = c(-18, 18))+
    scale_y_continuous(limits = c(-18, 18))+
    theme_bw()+
    theme( axis.text.x = element_text( face = "bold"),
           strip.text = element_text(face = "bold"),
           plot.title = element_text(hjust = 0.5, size=20),
           legend.position="bottom")  
  grid.arrange(gg)
  
  gg <- ggplot(tsne.res, aes(x=`t-SNE1`, y=`t-SNE2`,col=Genotype)) + 
    geom_point( size=4)+ggtitle("All Samples") +
    scale_color_manual(values =geno_cols , aesthetics = c("colour", "fill"))+
    scale_x_continuous(limits = c(-18, 18))+
    scale_y_continuous(limits = c(-18, 18))+
    theme_bw()+
    theme( axis.text.x = element_text( face = "bold"),
           strip.text = element_text(face = "bold"),
           plot.title = element_text(hjust = 0.5, size=20),
           legend.position="bottom")  
  grid.arrange(gg)
  
  gg <- ggplot(tsne.res, aes(x=`t-SNE1`, y=`t-SNE2`,col=Microbiota)) + 
    geom_point( size=4)+ggtitle("All Samples") +
    scale_color_manual(values =microb_cols , aesthetics = c("colour", "fill"))+
    scale_x_continuous(limits = c(-18, 18))+
    scale_y_continuous(limits = c(-18, 18))+
    theme_bw()+
    theme( axis.text.x = element_text( face = "bold"),
           strip.text = element_text(face = "bold"),
           plot.title = element_text(hjust = 0.5, size=20),
           legend.position="bottom")  
  grid.arrange(gg)
  #Plotting tSNE - all samples together , for each Material coloring tissues
  temp <- tsne.res
  for ( x in unique(tsne.res$Material) ){
    temp$color <- ifelse(temp$Material == x , temp$Tissue, paste0("Non-",x))
    gg <- ggplot(temp, aes(x=`t-SNE1`, y=`t-SNE2`,col=color)) + 
      geom_point( size=4)+ggtitle(x) +
      scale_color_manual(breaks = c(names(tissue_cols),paste0("Non-",x)) , 
                         values =c (tissue_cols,"lightgrey") %>% setNames(.,c(names(tissue_cols),paste0("Non-",x))) , 
                         aesthetics = c("colour", "fill") )+
      scale_x_continuous(limits = c(-18, 18))+
      scale_y_continuous(limits = c(-18, 18))+
      labs(color='Tissue') +
      theme_bw()+
      theme( axis.text.x = element_text( face = "bold"),
             strip.text = element_text(face = "bold"),
             plot.title = element_text(hjust = 0.5, size=20),
             legend.position="bottom")  
    grid.arrange(gg)
  }
  
  temp <- tsne.res
  #Plotting tSNE - all samples together , for each tissue coloring is done using genotype, microbiota or their interactions
  for ( x in unique(tsne.res$Tissue)){
    temp$color <- ifelse(temp$Tissue == x , temp$Microbiota %>% as.character, paste0("Non-",x))
    gg <- ggplot(temp, aes(x=`t-SNE1`, y=`t-SNE2`,col=color)) + 
      geom_point( size=4)+ggtitle(x) +
      scale_color_manual(breaks = c(names(microb_cols),paste0("Non-",x)) , 
                         values =c (microb_cols,"lightgrey") %>% setNames(.,c(names(microb_cols),paste0("Non-",x))) , 
                         aesthetics = c("colour", "fill") )+
      scale_x_continuous(limits = c(-18, 18))+
      scale_y_continuous(limits = c(-18, 18))+
      labs(color='Microbiota') +
      theme_bw()+
      theme( axis.text.x = element_text( face = "bold"),
             strip.text = element_text(face = "bold"),
             plot.title = element_text(hjust = 0.5, size=20),
             legend.position="bottom")  
    grid.arrange(gg)
    
    temp$color <- ifelse(temp$Tissue == x , temp$Genotype %>% as.character, paste0("Non-",x))
    gg <- ggplot(temp, aes(x=`t-SNE1`, y=`t-SNE2`,col=color),
                 legend.position="bottom") + 
      geom_point( size=4)+ggtitle(x) +
      scale_color_manual(breaks = c(names(geno_cols),paste0("Non-",x)) , 
                         values =c (geno_cols,"lightgrey") %>% setNames(.,c(names(geno_cols),paste0("Non-",x))) , 
                         aesthetics = c("colour", "fill") )+
      scale_x_continuous(limits = c(-18, 18))+
      scale_y_continuous(limits = c(-18, 18))+
      labs(color='Genotype') +
      theme_bw()+
      theme( axis.text.x = element_text( face = "bold"),
             strip.text = element_text(face = "bold"),
             plot.title = element_text(hjust = 0.5, size=20),
             legend.position="bottom")  
    grid.arrange(gg)
    
    temp$color <- ifelse(temp$Tissue == x , paste(temp$Genotype,temp$Microbiota) %>% as.character, paste0("Non-",x))
    gg <- ggplot(temp, aes(x=`t-SNE1`, y=`t-SNE2`,col=color)) + 
      geom_point( size=4)+
      scale_color_manual(breaks = c(names(G_by_M_cols),paste0("Non-",x)) , 
                         values =c (G_by_M_cols,"lightgrey") %>% setNames(.,c(names(G_by_M_cols),paste0("Non-",x))) , 
                         aesthetics = c("colour", "fill") )+
      
      ggtitle(x) +
      scale_x_continuous(limits = c(-18, 18))+
      scale_y_continuous(limits = c(-18, 18))+
      labs(color='Genotype:Microbiota') +
      theme_bw()+
      theme( axis.text.x = element_text( face = "bold"),
             strip.text = element_text(face = "bold"),
             plot.title = element_text(hjust = 0.5, size=20),
             legend.position="bottom")  
    grid.arrange(gg)
  }
  #####End tissue specific plotting of tSNE results  ###################################
  dev.off()
  ###################################################################################
  ##################End tSNE plots
  ###################################################################################
  
  ######################################################################################################################################################################
  #################Conditional Effects For Each Tissue , 
  ################# ASO vs WT/ SPF, ASO vs WT/ GF, GF vs SPF/ ASO, GF vs SPF/ WT
  ######################################################################################################################################################################
  pdf(file="./figures/August 2021/manhattan_conditional_August21.pdf" ,onefile=T  , height = 16, width = 20)
  
  microbiota.res <- read.xlsx("./results/caltech.PD mice.040820.xlsx","conditional.microbita.effect") %>%
    subset(imputation == "imputed" )
  genotype.res <- read.xlsx("./results/caltech.PD mice.040820.xlsx","conditional.genotype.effect") %>%
    subset(imputation == "imputed" )
  ldply(microbiota.res$Tissue %>% unique %>% as.list %>% setNames(.,microbiota.res$Tissue %>% unique) , function(x){
    if ( x%in% c("Brainstem","Cortex","Nigra","Striatum")){
      material ="Brain Tissue"
    }else if ( x%in% c("Cecum","Colon Content","Duodenum Content")){
      material ="Gut Content"
    }else if ( x%in% c("Colon","Duodenum")){
      material ="Intestinal Tract Tissue"
    }else{
      material <- "Plasma"
    }
    print (c(x,":",material))
    
    axisdf <- data.frame (Metabolite= mets.included[[material]] , mets.metadata[mets.included[[material]],c("Analyte.Class","Metabolite.Name")]) %>%
      # Compute number of analyte for each class
      group_by(Analyte.Class) %>% 
      dplyr::summarise(n=n()) %>% 
      # Calculate cumulative position of each analyte class
      mutate(gap=ifelse(n>6,0,3),tot=cumsum(n)-n+cumsum(gap),  start = tot+1, end=tot+n , center = round((start+end)/2) ) %>%
      dplyr::select(-gap)
    axisdf$col=rep(c("grey", "white"), len = nrow(axisdf))
    
    plotdf <- data.frame (Metabolite= mets.included[[material]] , mets.metadata[mets.included[[material]],c("Analyte.Class","Metabolite.Name")]) %>%
      ddply(., .(Analyte.Class), function(x){  
        x$MetNum = 1:nrow(x)
        return (x)
      }) %>% 
      right_join(., 
                 rbind.fill(
                   microbiota.res %>% subset(Tissue == x) %>% dplyr::select(Metabolite , Genotype , t,ttest.p ,wilcox.p) %>% mutate(Effect=paste0("GF vs. SPF|",Genotype)),
                   genotype.res   %>% subset(Tissue == x) %>% dplyr::select(Metabolite , Microbiota, t,ttest.p ,wilcox.p) %>% mutate(Effect=paste0("WT vs. ASO|",Microbiota))
                 )
                 , by = "Metabolite") %>%
      left_join(., axisdf , by = "Analyte.Class") %>%
      mutate(x = MetNum + tot) %>%
      # Add highlight and annotation information
      mutate( is_highlight=ifelse(wilcox.p<0.05, "yes", "no")) 
      #mutate( is_annotate=ifelse(wilcox.p<0.05 , "yes", "no"))       
    
    #For each effect type find top 25 significant metabolites having largest changes and smallest p-values
    plotdf <- ddply(plotdf , .(Effect) , function(X){
      ord <- with(X , order(abs(t)/wilcox.p, decreasing = T))
      topN = X[ord , ] %>% subset(wilcox.p <0.05) %>% head(25) %>% dplyr::select(Metabolite) %>% unlist
      
      X$is_annotate=ifelse(X$Metabolite %in% topN , "yes", "no")
      return(X)
    })
    
    p<- ggplot(plotdf ,  aes(x=x, y=-sign(t)*log10(wilcox.p))) +
      facet_wrap(~Effect , ncol=2)+
      # Show all points
      #geom_point( aes(color=as.factor(Analyte.Class)), alpha=0.8, size=1.3) +
      geom_point(size=1.3)+
      #scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
      
      # custom X axis:
      scale_x_continuous( label = axisdf$Analyte.Class, breaks= axisdf$center ) +
      #scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
      #add rectangles to distinguish classes
      geom_rect(data = axisdf,
                aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill=col),
                alpha = 0.5  , inherit.aes = FALSE)+    
      scale_fill_manual(values = alpha(c("grey", "white"), .3))+
      #add line to show significance
      geom_hline(yintercept=c(log10(0.05), -log10(0.05)), linetype="dashed", color = "red")+
      #geom_hline(yintercept=log10(0.05), linetype="dashed", color = "red")+
      
      # Add highlighted points
      geom_point(data=subset(plotdf, is_highlight=="yes"), color="orange", size=2) +
      
      # Add label using ggrepel to avoid overlapping
      #geom_label_repel( data=subset(plotdf, is_annotate=="yes"), aes(label=Metabolite), size=2) +
      ggtitle(x )+  xlab("") +
      geom_label_repel( data=subset(plotdf, is_annotate=="yes"), 
                        aes(label=Metabolite),box.padding = unit(0.5, "lines"), size=6)+

      # Custom the theme:
      theme_bw()+
      theme( 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text( face = "bold",size = 12,  angle = 90, vjust = 0.5, hjust=1),
        strip.text = element_text(face = "bold", size = 18),
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        legend.position = "none")  
    grid.arrange(p )
    
    
    return (1)
  })
  
  dev.off()


  ################################################################################################################
  #################Marginal Plots-   For Each Tissue , Plot GF vs SPF ; ASO vs WT
  ################################################################################################################
  pdf(file="./figures/August 2021/manhattan_marginal_August21.pdf" ,onefile=T  , height = 16, width = 14)
  microbiota.res <- read.xlsx("./results/caltech.PD mice.040820.xlsx","marginal.microbita.effect") %>%
    subset(imputation == "imputed" )
  genotype.res <- read.xlsx("./results/caltech.PD mice.040820.xlsx","marginal.genotype.effect") %>%
    subset(imputation == "imputed" )
  ldply(microbiota.res$Tissue %>% unique %>% as.list %>% setNames(.,microbiota.res$Tissue %>% unique) , function(x){
    if ( x%in% c("Brainstem","Cortex","Nigra","Striatum")){
      material ="Brain Tissue"
    }else if ( x%in% c("Cecum","Colon Content","Duodenum Content")){
      material ="Gut Content"
    }else if ( x%in% c("Colon","Duodenum")){
      material ="Intestinal Tract Tissue"
    }else{
      material <- "Plasma"
    }
    print (c(x,":",material))
    
    axisdf <- data.frame (Metabolite= mets.included[[material]] , mets.metadata[mets.included[[material]],c("Analyte.Class","Metabolite.Name")]) %>%
      # Compute number of analyte for each class
      group_by(Analyte.Class) %>% 
      dplyr::summarise(n=n()) %>% 
      # Calculate cumulative position of each analyte class
      mutate(gap=ifelse(n>6,0,3),tot=cumsum(n)-n+cumsum(gap),  start = tot+1, end=tot+n , center = round((start+end)/2) ) %>%
      dplyr::select(-gap)
    axisdf$col=rep(c("grey", "white"), len = nrow(axisdf))
    
    plotdf <- data.frame (Metabolite= mets.included[[material]] , mets.metadata[mets.included[[material]],c("Analyte.Class","Metabolite.Name")]) %>%
      ddply(., .(Analyte.Class), function(x){  
        x$MetNum = 1:nrow(x)
        return (x)
      }) %>% 
      right_join(., 
                 rbind(
                   microbiota.res %>% subset(Tissue == x) %>% dplyr::select(Metabolite , t ,wilcox.p) %>% mutate(Effect="Microbiota (GF vs. SPF)"),
                   genotype.res   %>% subset(Tissue == x) %>% dplyr::select(Metabolite , t ,wilcox.p) %>% mutate(Effect="Genotype (WT vs. ASO)")
                 )
                 , by = "Metabolite") %>%
      left_join(., axisdf , by = "Analyte.Class") %>%
      mutate(x = MetNum + tot) %>%
      # Add highlight and annotation information
      mutate( is_highlight=ifelse(wilcox.p<0.05, "yes", "no"))
      #mutate( is_annotate=ifelse(wilcox.p<0.05 , "yes", "no"))       
    #For each effect type, find top 25 significant metabolites having largest changes and smallest p-values
    plotdf <- ddply(plotdf , .(Effect) , function(X){
      ord <- with(X , order(abs(t)/wilcox.p, decreasing = T))
      topN = X[ord , ] %>% subset(wilcox.p <0.05) %>% head(25) %>% dplyr::select(Metabolite) %>% unlist
      
      X$is_annotate=ifelse(X$Metabolite %in% topN , "yes", "no")
      return(X)
    })
    
    p<- ggplot(plotdf ,  aes(x=x, y=-sign(t)*log10(wilcox.p))) +
      facet_wrap(~Effect , ncol=1)+
      # Show all points
      #geom_point( aes(color=as.factor(Analyte.Class)), alpha=0.8, size=1.3) +
      geom_point(size=1.3)+
      #scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
      
      # custom X axis:
      scale_x_continuous( label = axisdf$Analyte.Class, breaks= axisdf$center ) +
      #scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
      #add rectangles to distinguish classes
      geom_rect(data = axisdf,
                aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill=col),
                alpha = 0.5  , inherit.aes = FALSE)+    
      scale_fill_manual(values = alpha(c("grey", "white"), .3))+
      #add line to show significance
      geom_hline(yintercept=c(log10(0.05), -log10(0.05)), linetype="dashed", color = "red")+
      #geom_hline(yintercept=log10(0.05), linetype="dashed", color = "red")+
      
      # Add highlighted points
      geom_point(data=subset(plotdf, is_highlight=="yes"), color="orange", size=2) +
      
      # Add label using ggrepel to avoid overlapping
      #geom_label_repel( data=subset(plotdf, is_annotate=="yes"), aes(label=Metabolite), size=2) +
      ggtitle(x )+  xlab("") +
      geom_label_repel( data=subset(plotdf, is_annotate=="yes"), 
                        aes(label=Metabolite),box.padding = unit(0.5, "lines"), size=6)+
      # Custom the theme:
      theme_bw()+
      theme( 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text( face = "bold", angle = 90, vjust = 0.5, hjust=1 , size = 12),
        strip.text = element_text(face = "bold", size = 18),
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        legend.position = "none")  
    grid.arrange(p )
    
    
    return (1)
  })
  
  dev.off()
  

  ################################################################################################################
  ##### Volcano plots
  ################################################################################################################
  pdf(file="./figures/August 2021/volcano_conditional_August21.pdf" ,onefile=T  , height = 12, width = 14)
  microbiota.res <- read.xlsx("./results/caltech.PD mice.040820.xlsx","conditional.microbita.effect") %>%
    subset(imputation == "imputed" )
  genotype.res <- read.xlsx("./results/caltech.PD mice.040820.xlsx","conditional.genotype.effect") %>%
    subset(imputation == "imputed" )
  ldply(microbiota.res$Tissue %>% unique %>% as.list %>% setNames(.,microbiota.res$Tissue %>% unique) , function(x){
    if ( x%in% c("Brainstem","Cortex","Nigra","Striatum")){
      material ="Brain Tissue"
    }else if ( x%in% c("Cecum","Colon Content","Duodenum Content")){
      material ="Gut Content"
    }else if ( x%in% c("Colon","Duodenum")){
      material ="Intestinal Tract Tissue"
    }else{
      material <- "Plasma"
    }
    print (c(x,":",material))
    
    
    plotdf <-  rbind.fill(
      microbiota.res %>% subset(Tissue == x) %>% dplyr::select(Metabolite , Genotype , t,ttest.p ,wilcox.p) %>% mutate(Effect=paste0("GF vs. SPF|",Genotype)),
      genotype.res   %>% subset(Tissue == x) %>% dplyr::select(Metabolite , Microbiota, t,ttest.p ,wilcox.p) %>% mutate(Effect=paste0("WT vs. ASO|",Microbiota))
      ) %>% 
      left_join(., mets.metadata %>% dplyr::select(Analyte.Class , Abbreviation), 
                by = c("Metabolite"="Abbreviation") )
    
    tCcutoff = mean(plotdf$t)+2*sd(plotdf$t)   #highlight those metabolites deviate more than 2sd from mean change
    
    plotdf <- plotdf %>%
      # Add highlight and annotation information
      mutate(Sig = case_when(
        abs(t)>tCcutoff  & wilcox.p <0.05~ "p-value and t",
        abs(t)>tCcutoff  ~ "t",
        wilcox.p <0.05  ~ "p-value",
        TRUE ~ "NS"
      ) %>% factor(. , levels = c("NS","t","p-value" , "p-value and t")))
    #For each effect type find top 25 significant metabolites having largest changes and smallest p-values
    plotdf <- ddply(plotdf , .(Effect) , function(X){
      ord <- with(X , order(abs(t)/wilcox.p, decreasing = T))
      topN = X[ord , ] %>% subset(wilcox.p <0.05) %>% head(20) %>% dplyr::select(Metabolite) %>% unlist
      
      X$is_annotate=ifelse(X$Metabolite %in% topN , "yes", "no")
      return(X)
    }) %>% 
      mutate(#replace metabolites with abs(t) > mean(t)+ 2sd with mean(t)+ 2sd
        Metabolite = ifelse(abs(t)>(mean(t)+3*sd(t)) ,  paste0(Metabolite,"*") , Metabolite),
        t = ifelse(abs(t)>(mean(t)+3*sd(t)) ,  mean(t)+sign(t)*3*sd(t) , t)
      )  
    
    
    print(range(plotdf$t))
    l = max(abs(plotdf$t)) %>% ceil  #l shoulbd around 3sd
    x_breaks <- c(-l , -tCcutoff,0,tCcutoff , l  ) %>% round(2)
    print(x_breaks)
    p<- ggplot(plotdf ,  aes(x=t, y=-log10(wilcox.p) , color = Sig)) +
      facet_wrap(~Effect , ncol=2)+
      geom_point(size=4)+
      scale_y_continuous(limits = c(0,max(-log10(plotdf$wilcox.p))+2) ) +     
      scale_x_continuous(limits = c(-l,l) ,breaks=x_breaks, expand = c(0, 0) ) +     
      #add rectangles to distinguish classes   
      scale_color_manual(values = c("NS"= "grey30","t"= "forestgreen","p-value"= "royalblue","p-value and t"= "red2"))+
      #add line to show significance
      geom_hline(yintercept= -log10(0.05), linetype="longdash", color = "black")+
      geom_vline(xintercept=c(-tCcutoff , tCcutoff), linetype="longdash", color = "black")+
      
      # Add highlighted points
      geom_point(alpha = 0.5 ,  size=2) +
      
      # Add label using ggrepel to avoid overlapping
      #geom_label_repel( data=subset(plotdf, is_annotate=="yes"), aes(label=Metabolite), size=2) +
      ggtitle(x)+  xlab("t") + ylab(expression("log"  [10] * "Wilcoxon P")  ) +
      geom_label_repel( data=subset(plotdf, is_annotate=="yes"), 
                        aes(label=Metabolite),color="black",  box.padding = unit(0.5, "lines"), size=6)+
      
      # Custom the theme:
      theme_bw()+
      theme( 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text( face = "bold"),
        strip.text = element_text(face = "bold" , size = 18),
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))  
    grid.arrange(p )
    
    
    return (1)
  })
  dev.off()
  
  ######Volcano plots for marginal effects
  ################################################################################################################
  pdf(file="./figures/August 2021/volcano_marginal_August21.pdf" ,onefile=T  , height = 12, width = 14)
  microbiota.res <- read.xlsx("./results/caltech.PD mice.040820.xlsx","marginal.microbita.effect") %>%
    subset(imputation == "imputed" )
  genotype.res <- read.xlsx("./results/caltech.PD mice.040820.xlsx","marginal.genotype.effect") %>%
    subset(imputation == "imputed" )
  ldply(microbiota.res$Tissue %>% unique %>% as.list %>% setNames(.,microbiota.res$Tissue %>% unique) , function(x){
    if ( x%in% c("Brainstem","Cortex","Nigra","Striatum")){
      material ="Brain Tissue"
    }else if ( x%in% c("Cecum","Colon Content","Duodenum Content")){
      material ="Gut Content"
    }else if ( x%in% c("Colon","Duodenum")){
      material ="Intestinal Tract Tissue"
    }else{
      material <- "Plasma"
    }
    print (c(x,":",material))
    
    plotdf <-  rbind.fill(
      microbiota.res %>% subset(Tissue == x) %>% dplyr::select(Metabolite , t ,wilcox.p) %>% mutate(Effect="Microbiota (GF vs. SPF)"),
      genotype.res   %>% subset(Tissue == x) %>% dplyr::select(Metabolite , t ,wilcox.p) %>% mutate(Effect="Genotype (WT vs. ASO)")
    ) %>% 
      left_join(., mets.metadata %>% dplyr::select(Analyte.Class , Abbreviation), 
                by = c("Metabolite"="Abbreviation") )
    tCcutoff = mean(plotdf$t)+2*sd(plotdf$t)   #highlight those metabolites deviate more than 2sd from mean change
    
    plotdf <- plotdf %>%
      # Add highlight and annotation information
      mutate(Sig = case_when(
        abs(t)>tCcutoff  & wilcox.p <0.05~ "p-value and t",
        abs(t)>tCcutoff  ~ "t",
        wilcox.p <0.05  ~ "p-value",
        TRUE ~ "NS"
      ) %>% factor(. , levels = c("NS","t","p-value" , "p-value and t")))
    #For each effect type find top 25 significant metabolites having largest changes and smallest p-values
    plotdf <- ddply(plotdf , .(Effect) , function(X){
      ord <- with(X , order(abs(t)/wilcox.p, decreasing = T))
      topN = X[ord , ] %>% subset(wilcox.p <0.05) %>% head(20) %>% dplyr::select(Metabolite) %>% unlist
      
      X$is_annotate=ifelse(X$Metabolite %in% topN , "yes", "no")
      return(X)
    }) %>% 
      mutate(#replace metabolites with abs(t) > mean(t)+ 3sd with mean(t)+ 2sd
        Metabolite = ifelse(abs(t)>(mean(t)+3*sd(t)) ,  paste0(Metabolite,"*") , Metabolite),
        t = ifelse(abs(t)>(mean(t)+3*sd(t)) ,  mean(t)+sign(t)*3*sd(t) , t)
      )  
    
    
    print(range(plotdf$t))
    l = max(abs(plotdf$t)) %>% ceil  #l shoulbd around 3sd
    x_breaks <- c(-l , -tCcutoff,0,tCcutoff , l  ) %>% round( 2)
    print(x_breaks)
    p<- ggplot(plotdf ,  aes(x=t, y=-log10(wilcox.p) , color = Sig)) +
      facet_wrap(~Effect , ncol=2)+
      geom_point(size=4)+
      scale_y_continuous(limits = c(0,max(-log10(plotdf$wilcox.p))+2) ) +     
      scale_x_continuous(limits = c(-l,l) ,breaks=x_breaks, expand = c(0, 0) ) +     
      #add rectangles to distinguish classes
      scale_color_manual(values = c("NS"= "grey30","t"= "forestgreen","p-value"= "royalblue","p-value and t"= "red2"))+
      #add line to show significance
      geom_hline(yintercept= -log10(0.05), linetype="longdash", color = "black")+
      geom_vline(xintercept=c(-tCcutoff , tCcutoff), linetype="longdash", color = "black")+
      
      # Add highlighted points
      geom_point(alpha = 0.5 ,  size=2) +
      
      # Add label using ggrepel to avoid overlapping
      #geom_label_repel( data=subset(plotdf, is_annotate=="yes"), aes(label=Metabolite), size=2) +
      ggtitle(x)+  xlab("t") + ylab(expression("log"  [10] * "Wilcoxon P")  ) +
      geom_label_repel( data=subset(plotdf, is_annotate=="yes"), 
                        aes(label=Metabolite),color="black",  box.padding = unit(0.5, "lines"), size=6)+
      
      # Custom the theme:
      theme_bw()+
      theme( 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text( face = "bold"),
        strip.text = element_text(face = "bold" , size = 18),
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))  
    grid.arrange(p )
    
    
    return (1)
  })
  dev.off()
  
  