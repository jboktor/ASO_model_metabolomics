#####################################################################################
### Script for regression analysis of Cal Tech PD Mouse Model - Q500
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

setwd("/users/siamakdehkordi/Desktop/Duke/Cal Tech PD Mice/")
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
  
  #Find Gut Metabolites that (those that  expressed only in SPF but not in GF)
  mets.SPFonly = dlply(dta ,.(Tissue) , function(x){
    percent.b.LOD <- ddply(x , .(Microbiota), function(xx){
      apply(xx %>% dplyr::select(mets), 2, function(xxx){
        mean(!is.na(xxx))
      })      
    })
    idx <- percent.b.LOD [percent.b.LOD$Microbiota=="GF" , -1] <.3 & percent.b.LOD [percent.b.LOD$Microbiota=="SPF" , -1] >=.5
    
    return (colnames(percent.b.LOD[-1])[idx])
  })
  
  #Find metabolites that survived <LOD
  common_mets <-  Reduce( intersect , mets.included )



  dta$Genotype <- factor(dta$Genotype , levels = c("WT","ASO"))
  dta$Microbiota <- factor(dta$Microbiota , levels = c("GF","SPF"))
  #################impute data using LOD /2
  extracols <- setdiff(colnames(dta) , mets)
  dta.imputed <- left_join(dta , 
                           LOD %>% dplyr::select(-`Plate Number FIA`),
                           by = "Plate Number LC",
                           suffix = c("", ".LOD"))
  #impute missing values with LOD/2
  for ( met in mets){
    idx <- which(is.na(dta.imputed[,met]))
    dta.imputed[idx,met] <- dta.imputed[idx,paste0(met,".LOD")] / 2
  }
  dta.imputed <- dta.imputed %>% dplyr::select(-contains("LOD"))
  ###End of Imputation
  
  
  
  anal.dta <- rbind(dta %>% mutate(imputation = "unimputed"),
                    dta.imputed %>%mutate(imputation = "imputed"))
  
  runTSNE = F
  if (runTSNE){
    set.seed(40)  
    #pdf(file="./figures/tsne_caltech_jan2021_lod2Imputed.pdf" ,onefile=T  , height = 6, width = 6.5)
    #pdf(file="./figures/tsne_caltech_jan2021_lod2Imputed_april 2021.pdf" ,onefile=T  , height = 6, width = 6.5)
    pdf(file="./figures/tsne_caltech_jan2021_lod2Imputed_August 2021.pdf" ,onefile=T  , height = 6, width = 6.5)
    #Run tSNE on all samples 
    #tsne_caltech_mice <-  tsne(dta[,common_mets] %>% impute.knn.obs.sel(.)  %>% scale , perplexity=50)
    tsne_caltech_mice <-  tsne(dta.imputed[,common_mets] %>% log2  %>% scale , perplexity=50)
    tsne_caltech_mice <- data.frame(tsne_caltech_mice)
    colnames(tsne_caltech_mice) <- c("t-SNE1","t-SNE2")
    tsne_caltech_mice$Material <- dta$Material
    tsne_caltech_mice$Genotype <- dta$Genotype
    tsne_caltech_mice$Microbiota <- dta$Microbiota
    tsne_caltech_mice$Tissue <- dta$Tissue
    tsne_caltech_mice$G_by_M <- paste(tsne_caltech_mice$Genotype , tsne_caltech_mice$Microbiota)
    
    
    geno_cols <-c("ASO" = "#1f78b4", "WT" = "#33a02c")
    micro_cols <-c("GF" = "#1f78b4","SPF" = "#33a02c")
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

    gg <- ggplot(tsne_caltech_mice, aes(x=`t-SNE1`, y=`t-SNE2`,col=Material)) + 
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
    
    
    gg <- ggplot(tsne_caltech_mice, aes(x=`t-SNE1`, y=`t-SNE2`,col=Genotype)) + 
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
    
    gg <- ggplot(tsne_caltech_mice, aes(x=`t-SNE1`, y=`t-SNE2`,col=Microbiota)) + 
      geom_point( size=4)+ggtitle("All Samples") +
      scale_color_manual(values =micro_cols , aesthetics = c("colour", "fill"))+
      scale_x_continuous(limits = c(-18, 18))+
      scale_y_continuous(limits = c(-18, 18))+
      theme_bw()+
      theme( axis.text.x = element_text( face = "bold"),
             strip.text = element_text(face = "bold"),
             plot.title = element_text(hjust = 0.5, size=20),
             legend.position="bottom")  
    grid.arrange(gg)
    
    library(ggforce)
    temp <- tsne_caltech_mice
    for ( x in unique(tsne_caltech_mice$Material)){
      temp$color <- ifelse(temp$Material == x , temp$Microbiota %>% as.character, NA)
      gg <- ggplot(temp, aes(x=`t-SNE1`, y=`t-SNE2`,col=color)) + 
        geom_point( size=4)+ggtitle(x) +
        scale_color_manual(values =micro_cols , aesthetics = c("colour", "fill") , na.value="lightgrey")+
        scale_x_continuous(limits = c(-18, 18))+
        scale_y_continuous(limits = c(-18, 18))+
        geom_mark_ellipse(data = temp %>% subset(Material == x), aes(color = Tissue,
                              label=Tissue),
                          expand = unit(0.5,"mm"),
                          label.buffer = unit(-5, 'mm'))+        
        theme_bw()+
        theme( axis.text.x = element_text( face = "bold"),
               strip.text = element_text(face = "bold"),
               plot.title = element_text(hjust = 0.5, size=20),
               legend.position="bottom")  
      grid.arrange(gg)
      
      temp$color <- ifelse(temp$Tissue == x , temp$Genotype %>% as.character, NA)
      gg <- ggplot(temp, aes(x=`t-SNE1`, y=`t-SNE2`,col=color),
                   legend.position="bottom") + 
        geom_point( size=4)+ggtitle(x) +
        scale_color_manual(values =geno_cols , aesthetics = c("colour", "fill"), na.value="lightgrey")+
        scale_x_continuous(limits = c(-18, 18))+
        scale_y_continuous(limits = c(-18, 18))+
        theme_bw()+
        theme( axis.text.x = element_text( face = "bold"),
               strip.text = element_text(face = "bold"),
               plot.title = element_text(hjust = 0.5, size=20),
               legend.position="bottom")  
      grid.arrange(gg)
      
      temp$color <- ifelse(temp$Tissue == x , paste(temp$Genotype,temp$Microbiota) %>% as.character, NA)
      gg <- ggplot(temp, aes(x=`t-SNE1`, y=`t-SNE2`,col=color)) + 
        geom_point( size=4)+
        scale_color_manual(values =G_by_M_cols , aesthetics = c("colour", "fill"), na.value="lightgrey")+
        ggtitle(x) +
        scale_x_continuous(limits = c(-18, 18))+
        scale_y_continuous(limits = c(-18, 18))+
        theme_bw()+
        theme( axis.text.x = element_text( face = "bold"),
               strip.text = element_text(face = "bold"),
               plot.title = element_text(hjust = 0.5, size=20),
               legend.position="bottom")  
      grid.arrange(gg)
    }
    dev.off()
    
    
    #RUN PCA
    library(ggbiplot)
  }
  
  ############Run anova to see which metabolite is significant
  #1 SPF vs GF
  #2- SPF ASO x SPF WT 2)SPF ASO x GF ASO 3) GF WT x SPF WT 4) GF WT x GF ASO
  #1- Blood and Brain
  
  runregression=F
  if (runregression){
    regression.interaction = dlply(anal.dta, .(Tissue,imputation), function(dat){
      material <- unique(dat$Material)[1]
      tissue <- unique ( dat$Tissue)
      
      reg.res <- c()
      anova.res <- c()
      for ( met in mets.included[[material]]){
        print(sprintf("G x M Iteraction- %s:%s:%s",tissue,dat$imputation[1],met))
  
        tryCatch({
          f <- as.formula(sprintf("log(`%s`)~Microbiota*Genotype" , met))
          fit1 <- lm(f ,  dat)
          res.met = data.frame(Metabolite = met,
                               summary(fit1)$coef %>% data.frame(., check.names = F) %>% 
                                 mutate(Variable=rownames(.)) %>% cbind(.,confint(fit1)),check.names = F) 
          reg.res <- rbind(reg.res,res.met)
          
          res.met = data.frame(Metabolite = met ,N= fit1$model%>% nrow,anova(fit1) %>% 
                                 data.frame(., check.names = F) %>% mutate(Variable=rownames(.)),check.names = F)
          anova.res <- rbind(anova.res,res.met)
          
          
        },
        error = function(e){
          message(paste(material,":",paste0(tissue, collapse = ":"),"----",met,"\n"), e)})
      }
      reg.res <- left_join(reg.res , mets.metadata %>% dplyr::rename(Metabolite=Abbreviation) ,by = "Metabolite")
      anova.res <- left_join(anova.res , mets.metadata %>% dplyr::rename(Metabolite=Abbreviation) ,by = "Metabolite")
      return (list(regression = reg.res, anova= anova.res))
    })
    
    #Compare Microbiota effect ASO vs WT for different tissue types and Microbiota
    conditional.genotype.effect = ddply(anal.dta, .(imputation,Microbiota , Tissue ), function(dat){
      material <- dat$Material[1]
      tissue <- dat$Tissue[1]
      microbiota <- dat$Microbiota[1]
      
      diff.abund.test <- c()
      for ( met in mets.included[[material]]){
        print(sprintf("ASO vs WT- %s:%s:%s:%s",microbiota, tissue  ,dat$imputation[1],met))
        dat$Metabolite = dat[,met]
        tryCatch({
  
          idx <- !is.na(dat[,met])
          wilcox <- wilcox.test(dat[idx,met]~ dat[idx,"Genotype"])
          ttest <- t.test(log(dat[idx,met])~ dat[idx,"Genotype"])
          
          fisher.pitman = coin::oneway_test(Metabolite~as.factor(Genotype),data=dat[idx,],
                            distribution=approximate(nresample=20000))
          
          diff.abund.test <- rbind(diff.abund.test , 
                                   data.frame(Metabolite = met, Item = "ASO vs WT" , Description = "ASO - WT", 
                                              ttest.estimate= paste(names(ttest$estimate) , ":", round(ttest$estimate,3) , collapse = ";") , 
                                              t=ttest$statistic, ttest.df=ttest$parameter,
                                              ttest.confint=paste(c("Lower","Upper"),":", ttest$conf.int, collapse = ";"),
                                              ttest.p = ttest$p.value, 
                                              wilcox.p=wilcox$p.value , wilcox.W=wilcox$statistic,
                                              fisher.pitman.p=pvalue(fisher.pitman) , fisher.pitman.Z=statistic(fisher.pitman),
                                              check.names = F))
          
        },
        error = function(e){
          message(paste(material,":",paste0(tissue, collapse = ":"),"----",met,"\n"), e)})
      }
      diff.abund.test <- left_join(diff.abund.test , mets.metadata %>% dplyr::rename(Metabolite=Abbreviation) ,by = "Metabolite")
      return (diff.abund.test)
    })
    
    #Compare Microbiota effect : GF vs SPF for different tissue types and genotypes
    conditional.Microbiota.effect = ddply(anal.dta, .(imputation,Genotype , Tissue ), function(dat){
      material <- dat$Material[1]
      tissue <- dat$Tissue[1]
      microbiota <- dat$Microbiota[1]
      
      diff.abund.test <- c()
      for ( met in mets.included[[material]]){
        print(sprintf("GF vs SPF- %s:%s:%s:%s",microbiota, tissue  ,dat$imputation[1],met))
        dat$Metabolite = dat[,met]
        tryCatch({
          
          idx <- !is.na(dat[,met])
          wilcox <- wilcox.test(dat[idx,met]~ dat[idx,"Microbiota"])
          ttest <- t.test(log(dat[idx,met])~ dat[idx,"Microbiota"])
          
          fisher.pitman = coin::oneway_test(Metabolite~as.factor(Microbiota),data=dat[idx,],
                                            distribution=approximate(nresample=20000))
          
          diff.abund.test <- rbind(diff.abund.test , 
                                   data.frame(Metabolite = met, Item = "GF vs SPF" , Description = "GF - SPF", 
                                              ttest.estimate= paste(names(ttest$estimate) , ":", round(ttest$estimate,3) , collapse = ";") , 
                                              t=ttest$statistic, ttest.df=ttest$parameter,
                                              ttest.confint=paste(c("Lower","Upper"),":", ttest$conf.int, collapse = ";"),
                                              ttest.p = ttest$p.value, 
                                              wilcox.p=wilcox$p.value , wilcox.W=wilcox$statistic,
                                              fisher.pitman.p=pvalue(fisher.pitman) , fisher.pitman.Z=statistic(fisher.pitman),
                                              check.names = F))
          
        },
        error = function(e){
          message(paste(material,":",paste0(tissue, collapse = ":"),"----",met,"\n"), e)})
      }
      diff.abund.test <- left_join(diff.abund.test , mets.metadata %>% dplyr::rename(Metabolite=Abbreviation) ,by = "Metabolite")
      return (diff.abund.test)
    })
    
    #######################
    #Compare Microbiota effect ASO vs WT for different tissue types
    marginal.genotype.effect = ddply(anal.dta, .(imputation , Tissue ), function(dat){
      material <- dat$Material[1]
      tissue <- dat$Tissue[1]

      diff.abund.test <- c()
      for ( met in mets.included[[material]]){
        print(sprintf("ASO vs WT- %s:%s:%s", tissue  ,dat$imputation[1],met))
        dat$Metabolite = dat[,met]
        tryCatch({
          
          idx <- !is.na(dat[,met])
          wilcox <- wilcox.test(dat[idx,met]~ dat[idx,"Genotype"])
          ttest <- t.test(log(dat[idx,met])~ dat[idx,"Genotype"])
          
          fisher.pitman = coin::oneway_test(Metabolite~as.factor(Genotype),data=dat[idx,],
                                            distribution=approximate(nresample=20000))
          
          diff.abund.test <- rbind(diff.abund.test , 
                                   data.frame(Metabolite = met, Item = "ASO vs WT" , Description = "ASO - WT", 
                                              ttest.estimate= paste(names(ttest$estimate) , ":", round(ttest$estimate,3) , collapse = ";") , 
                                              t=ttest$statistic, ttest.df=ttest$parameter,
                                              ttest.confint=paste(c("Lower","Upper"),":", ttest$conf.int, collapse = ";"),
                                              ttest.p = ttest$p.value, 
                                              wilcox.p=wilcox$p.value , wilcox.W=wilcox$statistic,
                                              fisher.pitman.p=pvalue(fisher.pitman) , fisher.pitman.Z=statistic(fisher.pitman),
                                              check.names = F))
          
        },
        error = function(e){
          message(paste(material,":",paste0(tissue, collapse = ":"),"----",met,"\n"), e)})
      }
      diff.abund.test <- left_join(diff.abund.test , mets.metadata %>% dplyr::rename(Metabolite=Abbreviation) ,by = "Metabolite")
      return (diff.abund.test)
    })
    
    #Compare Microbiota effect : GF vs SPF for different tissue types and genotypes
    marginal.microbiota.effect = ddply(anal.dta, .(imputation , Tissue ), function(dat){
      material <- dat$Material[1]
      tissue <- dat$Tissue[1]

      diff.abund.test <- c()
      for ( met in mets.included[[material]]){
        print(sprintf("GF vs SPF- %s:%s:%s", tissue  ,dat$imputation[1],met))
        dat$Metabolite = dat[,met]
        tryCatch({
          
          idx <- !is.na(dat[,met])
          wilcox <- wilcox.test(dat[idx,met]~ dat[idx,"Microbiota"])
          ttest <- t.test(log(dat[idx,met])~ dat[idx,"Microbiota"])
          
          fisher.pitman = coin::oneway_test(Metabolite~as.factor(Microbiota),data=dat[idx,],
                                            distribution=approximate(nresample=20000))
          
          diff.abund.test <- rbind(diff.abund.test , 
                                   data.frame(Metabolite = met, Item = "GF vs SPF" , Description = "GF - SPF", 
                                              ttest.estimate= paste(names(ttest$estimate) , ":", round(ttest$estimate,3) , collapse = ";") , 
                                              t=ttest$statistic, ttest.df=ttest$parameter,
                                              ttest.confint=paste(c("Lower","Upper"),":", ttest$conf.int, collapse = ";"),
                                              ttest.p = ttest$p.value, 
                                              wilcox.p=wilcox$p.value , wilcox.W=wilcox$statistic,
                                              fisher.pitman.p=pvalue(fisher.pitman) , fisher.pitman.Z=statistic(fisher.pitman),
                                              check.names = F))
          
        },
        error = function(e){
          message(paste(material,":",paste0(tissue, collapse = ":"),"----",met,"\n"), e)})
      }
      diff.abund.test <- left_join(diff.abund.test , mets.metadata %>% dplyr::rename(Metabolite=Abbreviation) ,by = "Metabolite")
      return (diff.abund.test)
    })
    
    
    #######Run log ratio test to test the effect of weight on metabolites:
        #Full model : log(metabolite)~ BW + G * M
        #Reduced Model:  log(metabolite)~ G * M
    lrtest_BW_res <- ddply(anal.dta, .(Tissue,imputation), function(dat){
      material <- unique(dat$Material)[1]
      tissue <- unique ( dat$Tissue)
      res <- c()
      for ( met in mets.included[[material]]){
        print(sprintf("G x M Iteraction- %s:%s:%s",tissue,dat$imputation[1],met))
        
        tryCatch({
          f1 <- as.formula(sprintf("log(`%s`)~BW+Microbiota*Genotype" , met))
          f2 <- as.formula(sprintf("log(`%s`)~Microbiota*Genotype" , met))
          fit1 <- lm(f1 ,  dat)
          fit2 <- lm(f2 ,  dat)
          
          res.met = rbind(
            tidy(fit1) %>% mutate(Model  = "Full"),
            tidy(fit2) %>% mutate(Model  = "Reduced")
          ) %>% filter(term != "(Intercept)") %>% mutate(Type = "anova")
          lrtest_res <- lrtest(fit1, fit2)
          res.met <- rbind.fill(
            res.met , 
            data.frame(Type="lrtest" , `Pr(>Chisq)` = lrtest_res$`Pr(>Chisq)`[2] , check.names  = F)) %>% 
            mutate(Metabolite = met)
          res <- rbind(res , res.met)
        },
        error = function(e){
          message(paste(material,":",paste0(tissue, collapse = ":"),"----",met,"\n"), e)})
      }
      res <- left_join(res , mets.metadata %>% dplyr::rename(Metabolite=Abbreviation) ,by = "Metabolite")
      return (res)
    })
    
    #####save lrtest results
    left_join(
      lrtest_BW_res %>% subset(Tissue == "Plasma" & imputation=="imputed" & `Pr(>Chisq)`<0.05) %>% dplyr::select(Metabolite , `Pr(>Chisq)` , Analyte.Class)  , 
      regression.interaction$Plasma.imputed$anova %>% filter(`Pr(>F)`<0.05 & Variable=="Microbiota") %>% dplyr::select(Metabolite , `Pr(>F)`) %>% dplyr::rename(`p(Microbiota)_reduced`=`Pr(>F)`),
      by = "Metabolite"
    ) %>%
      left_join(
        .  , 
        regression.interaction$Plasma.imputed$anova %>% filter(`Pr(>F)`<0.05 & Variable=="Genotype") %>% dplyr::select(Metabolite , `Pr(>F)`) %>% dplyr::rename(`p(Genotype)_reduced`=`Pr(>F)`),
        by = "Metabolite"
      ) %>%
      left_join(
        .  , 
        regression.interaction$Plasma.imputed$anova %>% filter(`Pr(>F)`<0.05 & Variable=="Microbiota:Genotype") %>% dplyr::select(Metabolite , `Pr(>F)`) %>% dplyr::rename(`p(Microbiota:Genotype)_reduced`=`Pr(>F)`),
        by = "Metabolite"
      ) %>%
      left_join(
        .  , 
        regression.interaction$Plasma.imputed$anova %>% filter(`Pr(>F)`<0.05 & Variable=="Microbiota") %>% dplyr::select(Metabolite , `Pr(>F)`) %>% dplyr::rename(`p(Microbiota)_full`=`Pr(>F)`),
        by = "Metabolite"
      ) %>%
      left_join(
        .  , 
        regression.interaction$Plasma.imputed$anova %>% filter(`Pr(>F)`<0.05 & Variable=="Genotype") %>% dplyr::select(Metabolite , `Pr(>F)`) %>% dplyr::rename(`p(Genotype)_full`=`Pr(>F)`),
        by = "Metabolite"
      ) %>%
      left_join(
        .  , 
        regression.interaction$Plasma.imputed$anova %>% filter(`Pr(>F)`<0.05 & Variable=="Microbiota:Genotype") %>% dplyr::select(Metabolite , `Pr(>F)`) %>% dplyr::rename(`p(Microbiota:Genotype)full`=`Pr(>F)`),
        by = "Metabolite"
      )  %>% write.csv(. , "lrtest BW.csv")    
    ##############################################################################################################
    ##############################################################################################################
    
    
    ##############################################################################################################
    ##############################################################################################################
    ###Save Results
    ##############################################################################################################
    ##############################################################################################################
    wb <- createWorkbook()
    
    temp = ldply(regression.interaction , function(x){
      return(x$regression)
    })
    addWorksheet(wb,"interaction.regcoef")
    writeDataTable(wb, "interaction.regcoef" ,temp)
  
    temp = ldply(regression.interaction , function(x){
      return(x$anova)
    })
    addWorksheet(wb,"interaction.anova")
    writeDataTable(wb, "interaction.anova" ,temp)
    
    addWorksheet(wb,"conditional.microbita.effect")
    writeDataTable(wb, "conditional.microbita.effect" ,conditional.microbiota.effect)
  
    addWorksheet(wb,"conditional.genotype.effect")
    writeDataTable(wb, "conditional.genotype.effect" ,conditional.genotype.effect)

    addWorksheet(wb,"marginal.microbita.effect")
    writeDataTable(wb, "marginal.microbita.effect" ,marginal.microbiota.effect)
    
    addWorksheet(wb,"marginal.genotype.effect")
    writeDataTable(wb, "marginal.genotype.effect" , marginal.genotype.effect)
    
    addWorksheet(wb,"lrtest_BW__res")
    writeDataTable(wb, "lrtest_BW__res" , lrtest_BW_res)
    
    ##Change reference group for genotype (WT as the reference)
          #saveWorkbook(wb, "./results/caltech.PD mice.040820.xlsx", overwrite = TRUE) 
    #####Added loglikelihood ratio tests for significance of BW
    saveWorkbook(wb, "./results/caltech.PD mice.122622.xlsx", overwrite = TRUE) 
  }
  
  if (FALSE){
    pdf(file="./figures/Sep 22/boxplot_metabolite_by_genotype_micro.pdf" ,onefile=T , height=10 ,width=10)
      
    my_comparisons <- list( c("WT.GF", "ASO.GF"), c("WT.GF", "WT.SPF"), c("WT.GF", "ASO.SPF"),c("ASO.GF","WT.SPF"), c("ASO.GF","ASO.SPF"),c("WT.SPF","ASO.SPF")) 
    rr <-     ddply(dta.imputed , .(Tissue) , function(X){
      tissue <- X$Tissue[1]
      material <- X$Material[1]
      X$Grouping = interaction(X$Genotype , X$Microbiota)
      print(paste(material , "/",tissue))
      ldply(split(mets.included[[material]], ceiling(seq_along(mets.included[[material]])/9)), function(m){
        print ( "\n")
        print (m)
        XX <- X %>% dplyr::select(Grouping,Genotype,Microbiota , all_of(m)) %>% 
          pivot_longer(4:ncol(.) , names_to = "Metabolite")

        pp <- ggboxplot(XX, x = "Grouping", y = "value",palette="jco", 
                  color = "Microbiota", width = 0.8 ) + 
          facet_wrap(~Metabolite , scales = "free_y" , ncol = 3)+
          #geom_jitter(aes(shape = Genotype))+
          stat_compare_means(comparisons = my_comparisons)+
          stat_compare_means(method = "anova")+
          labs(title=ifelse(tissue==material , material , paste0(material,"/",tissue)),x="Genotype X Microbiota", y = "Concetration(micro M)")+
          theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
                panel.spacing=unit(1, "lines"),
                strip.background = element_rect(color = "black", size = 1),
                strip.text  = element_text(face = "bold" , size = 14),
                axis.text.x = element_text(face = "bold"),
                plot.title = element_text(hjust = 0.5 , size = 14 , face="bold") , 
                plot.subtitle = element_text(hjust = 0.5, size = 14, face = "bold" , color= "red")
          )
        grid.arrange(pp)
        return(1)
      })
      return(1)
      
      
    })
    dev.off()
    
    
    pdf(file="./figures/Sep 22/boxplot_metabolite_by_genotype_micro_2.pdf" ,onefile=T , height=10 ,width=10)
    
    my_comparisons <- list( c("WT.GF", "ASO.GF"), c("WT.GF", "WT.SPF"), c("WT.GF", "ASO.SPF"),c("ASO.GF","WT.SPF"), c("ASO.GF","ASO.SPF"),c("WT.SPF","ASO.SPF")) 
    for ( met in mets.metadata[order(mets.metadata$Analyte.Class),]$Abbreviation){
      print (met)
      
      XX <-  dta.imputed %>% dplyr::select(Tissue , Genotype,Microbiota , all_of(met)) %>% 
        pivot_longer(4:ncol(.) , names_to = "Metabolite") %>%
        mutate(Grouping = interaction(Genotype , Microbiota)) %>%
        mutate(Tissue = factor(Tissue , levels = c("Plasma" , "Brainstem" , "Cortex" , "Nigra" , "Striatum" , "Colon" , "Duodenum" , "Cecum" , "Colon Content" , "Duodenum Content")))
      
      pp <- ggboxplot(XX, x = "Grouping", y = "value",palette="jco", 
                      color = "Microbiota", width = 0.8 ) + 
        facet_wrap(~Tissue , scales = "free_y" , ncol = 3)+
        #geom_jitter(aes(shape = Genotype))+
        stat_compare_means(comparisons = my_comparisons)+
        labs(title=met,x="Genotype X Microbiota", y = "Concetration(micro M)")+
        theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
              panel.spacing=unit(1, "lines"),
              strip.background = element_rect(color = "black", size = 1),
              strip.text  = element_text(face = "bold" , size = 14),
              axis.text.x = element_text(face = "bold"),
              plot.title = element_text(hjust = 0.5 , size = 14 , face="bold") , 
              plot.subtitle = element_text(hjust = 0.5, size = 14, face = "bold" , color= "red")
        )
      grid.arrange(pp)
      
    }
    dev.off()
    

  }
  if (FALSE){  #PERMANOVA analysis
    library(vegan)
    set.seed(1360)
    perm_anova_res <-     ddply(dta.imputed , .(Tissue) , function(X){
      profile <- X %>% dplyr::select(all_of(mets.included[[X$Material[1]]])) %>% 
        scale %>%  data.frame %>%
        select_if(~sum(!is.na(.)) > 0)
      
      
      animals <- X %>% dplyr::select(Microbiota,Genotype)
      profile.dist <- vegdist(profile, method="euclidean")
      
      dispersion <- betadisper(profile.dist, group=interaction(animals$Genotype, animals$Microbiota))
      dispersion_test <- (permutest(dispersion))$tab %>% data.frame
      
      profile.res <- adonis2(profile ~ Microbiota*Genotype, data = animals, permutations = 10000, method="euclidean")
      overall.res <- adonis2(profile ~ Microbiota*Genotype, data = animals, permutations = 10000, method="euclidean", by = NULL)
      
      return(
        rbind.fill(
          dispersion_test %>% tibble::rownames_to_column (var = "term") %>% mutate(Type = "dispersion test"),
          profile.res %>% as.data.frame %>% tibble::rownames_to_column (var = "term") %>% mutate(Type = "PERMANOVA"),
          overall.res %>% as.data.frame %>% tibble::rownames_to_column (var = "term") %>% mutate(Type = "Overall")
          
        )
      )
    })
    
    wb <- createWorkbook()
    addWorksheet(wb,"perm_anova")
    writeDataTable(wb, "perm_anova" , perm_anova_res)
    
    saveWorkbook(wb, "./results/caltech.PD mice.PERMANOVA.xlsx", overwrite = TRUE) #added hampsy, hamsa, hamsom to analysis
    

    d1 <- 
      perm_anova_res %>% subset(Type == "Overall" & term %in% c("Model")) %>% 
      dplyr::select(Tissue ,term , F , R2, `Pr(>F)`) %>%
      mutate(F =sprintf("%0.02f", F) ,R2 =sprintf("%0.02f", R2) , `Pr(>F)` =sprintf("%0.02E", `Pr(>F)`) ) %>%
      tidyr::pivot_wider(id_cols=Tissue,
                         names_from = term,
                         values_from = c("F", "R2","Pr(>F)")
      )
    
    d2 <-
      perm_anova_res %>% subset(Type == "PERMANOVA" & term %in% c("Microbiota" , "Genotype" , "Microbiota:Genotype")) %>% 
      dplyr::select(Tissue , term , F , R2, `Pr(>F)`) %>%
      mutate(F =sprintf("%0.02f", F) ,R2 =sprintf("%0.02f", R2) , `Pr(>F)` =sprintf("%0.02E", `Pr(>F)`) ) %>%
      tidyr::pivot_wider(id_cols=Tissue,
                         names_from = term,
                         values_from = c("F", "R2","Pr(>F)")
                         )
    permanova_suptable <- left_join(d1 , d2 ) %>%
      mutate(Tissue = factor(Tissue , levels = c("Plasma" , "Brainstem" , "Cortex" , "Nigra" , "Striatum" , "Colon" , "Duodenum" , "Cecum" , "Colon Content" , "Duodenum Content")))
    permanova_suptable <- permanova_suptable[order(permanova_suptable$Tissue) , ]
    write.csv(permanova_suptable , "./results/caltech.PD mice.PERMANOVA suptable.csv")
  }
  
  if (FALSE){
    library(ggpubr)
    library(grid)
    library(ggpattern)
    
    library(ComplexHeatmap)
    library(circlize)
    library(tidyverse)
    
    #For each Tissue type, print number of metabolites/Class survived
    write.csv( 
      ldply(mets.included , function(X){
        temp <- ddply(mets.metadata , .(Analyte.Class), function(XX){
          rr =intersect( X ,  XX$Abbreviation)
          return(   paste0(length(rr),"/",nrow(XX))   )
        })
        return (temp)
      }),
      "./log/plasma_mets_survived.csv" , row.names = F
    )
    
    pdf(file="./figures/plasma_cluster_Jan21.pdf" ,onefile=T  , height = 12, width = 16)
    set.seed(40)
    robust_dist = function(x, y) {
      qx = quantile(x, c(0.1, 0.9))
      qy = quantile(y, c(0.1, 0.9))
      l = x > qx[1] & x < qx[2] & y > qy[1] & y < qy[2]
      x = x[l]
      y = y[l]
      sqrt(sum((x - y)^2))
    }
    
    dta.plasma = dta.imputed %>% subset(Material=="Plasma")
    rownames(dta.plasma) = dta.plasma$Animal
    
    mets_to_feedin = mets.included$Plasma
    class_colors <-  colorspace::rainbow_hcl(length(unique(mets.metadata[mets_to_feedin,]$Analyte.Class))) 
    names(class_colors) <- unique(mets.metadata[mets_to_feedin,]$Analyte.Class)
    
    dta.plasma%>% dplyr::select(mets_to_feedin) %>% scale %>%
    Heatmap(., name = "scaled concentrations", cluster_columns = FALSE,# turn off row clustering
            show_column_names = F,
              clustering_distance_rows = robust_dist,
            right_annotation = rowAnnotation(Genotype = dta.plasma$Genotype, Microbiota = dta.plasma$Microbiota , col = list(Genotype = c(2:3) %>% setNames(., c("ASO","WT")), Microbiota=4:5 %>% setNames(., c("GF","SPF")))),
            top_annotation = HeatmapAnnotation(Class = mets.metadata[mets_to_feedin,]$Analyte.Class , col = list(Class=class_colors))
              ) 
    
    mets_to_feedin =  setdiff( mets.included$Plasma,(mets.metadata %>% subset(Analyte.Class=="Triglycerides"))$Abbreviation)
    class_colors <-  colorspace::rainbow_hcl(length(unique(mets.metadata[mets_to_feedin,]$Analyte.Class))) 
    names(class_colors) <- unique(mets.metadata[mets_to_feedin,]$Analyte.Class)
    dta.plasma%>% dplyr::select(mets_to_feedin) %>% scale %>%
      Heatmap(., name = "scaled concentrations", cluster_columns = FALSE,# turn off row clustering
              show_column_names = F,
              clustering_distance_rows = robust_dist,
              right_annotation = rowAnnotation(Genotype = dta.plasma$Genotype, Microbiota = dta.plasma$Microbiota , col = list(Genotype = c(2:3) %>% setNames(., c("ASO","WT")), Microbiota=4:5 %>% setNames(., c("GF","SPF")))),
              top_annotation = HeatmapAnnotation(Class = mets.metadata[mets_to_feedin,]$Analyte.Class , col = list(Class=class_colors))
      ) 
    dev.off()
    
    pdf(file="./figures/plasma_manhattan_anova_Jan21.pdf" ,onefile=T  , height = 18, width = 14)

    anova.res <- read.xlsx("./results/caltech.PD mice.011920.xlsx","interaction.anova") %>%
      subset(Tissue =="Plasma" & imputation == "imputed" & Variable %in% c("Microbiota","Genotype","Microbiota:Genotype"))
      
    axisdf <- data.frame (Metabolite= mets.included$Plasma , mets.metadata[mets.included$Plasma,c("Analyte.Class","Metabolite.Name")]) %>%
      # Compute number of analyte for each class
      group_by(Analyte.Class) %>% 
      dplyr::summarise(n=n()) %>% 
      # Calculate cumulative position of each analyte class
      mutate(gap=ifelse(n>6,0,3),tot=cumsum(n)-n+cumsum(gap),  start = tot+1, end=tot+n , center = round((start+end)/2), col=c(rep(c("grey", "white"), 12 ),"grey")) %>%
      dplyr::select(-gap)

    plotdf <- data.frame (Metabolite= mets.included$Plasma , mets.metadata[mets.included$Plasma,c("Analyte.Class","Metabolite.Name")]) %>%
      ddply(., .(Analyte.Class), function(x){  
        x$MetNum = 1:nrow(x)
        return (x)
        #return (data.frame(Metabolite = x$Metabolite, MetNum = 1:nrow(x)))
      }) %>% 
      left_join(., 
                anova.res %>% 
                  dplyr::select(Metabolite ,Variable, `Pr(>F)`)
                , by = "Metabolite") %>%
      left_join(., axisdf , by = "Analyte.Class") %>%
      mutate(x = MetNum + tot) %>%
      # Add highlight and annotation information
      mutate( is_highlight=ifelse(`Pr(>F)`<0.05, "yes", "no")) %>%
      mutate( is_annotate=ifelse(`Pr(>F)`<0.01 , "yes", "no"))       

    p<- ggplot(plotdf ,  aes(x=x, y=-log10(`Pr(>F)`))) +
      facet_wrap(~Variable, ncol = 1 ,  scales='free')+
      # Show all points
      #geom_point( aes(color=as.factor(Analyte.Class)), alpha=0.8, size=1.3) +
      geom_point(size=1.3)+
      #scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +

      # custom X axis:
      scale_x_continuous( label = axisdf$Analyte.Class, breaks= axisdf$center ) +
      scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
      #add rectangles to distinguish classes
      geom_rect(data = axisdf,
                aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill=col),
                alpha = 0.5  , inherit.aes = FALSE)+    
      scale_fill_manual(values = alpha(c("grey", "white"), .3))+
      #add line to show significance
      geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red")+

      # Add highlighted points
      geom_point(data=subset(plotdf, is_highlight=="yes"), color="orange", size=2) +
      
      # Add label using ggrepel to avoid overlapping
      geom_label_repel( data=subset(plotdf, is_annotate=="yes"), aes(label=Metabolite), size=2) +
      ggtitle("Plasma - Anova Genotype x Microbiota" )+  xlab("") +
      # Custom the theme:
      theme_bw()+
      theme( 
             panel.grid.major.x = element_blank(),
             panel.grid.minor.x = element_blank(),
             axis.text.x = element_text( face = "bold", angle = 90, vjust = 0.5, hjust=1),
             strip.text = element_text(face = "bold"),
             plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
             legend.position = "none")  
    grid.arrange(p)
    dev.off()
    ##############################################################################################################
    ###############    ###############    Plot Microbiota Effect    ###############    ###############    ########
    ##############################################################################################################
    pdf(file="./figures/plasma_manhattan_microb_marginal_Jan21.pdf" ,onefile=T  , height = 8, width = 14)

    microbiota.res <- read.xlsx("./results/caltech.PD mice.011920.xlsx","marginal.microbita.effect") %>%
      subset(Tissue =="Plasma" & imputation == "imputed" )
    
    axisdf <- data.frame (Metabolite= mets.included$Plasma , mets.metadata[mets.included$Plasma,c("Analyte.Class","Metabolite.Name")]) %>%
      # Compute number of analyte for each class
      group_by(Analyte.Class) %>% 
      dplyr::summarise(n=n()) %>% 
      # Calculate cumulative position of each analyte class
      mutate(gap=ifelse(n>6,0,3),tot=cumsum(n)-n+cumsum(gap),  start = tot+1, end=tot+n , center = round((start+end)/2), col=c(rep(c("grey", "white"), 12 ),"grey")) %>%
      dplyr::select(-gap)
    
    plotdf <- data.frame (Metabolite= mets.included$Plasma , mets.metadata[mets.included$Plasma,c("Analyte.Class","Metabolite.Name")]) %>%
      ddply(., .(Analyte.Class), function(x){  
        x$MetNum = 1:nrow(x)
        return (x)
        #return (data.frame(Metabolite = x$Metabolite, MetNum = 1:nrow(x)))
      }) %>% 
      right_join(., 
                 microbiota.res %>% dplyr::select(-Analyte.Class,-Metabolite.Name)
                 , by = "Metabolite") %>%
      left_join(., axisdf , by = "Analyte.Class") %>%
      mutate(x = MetNum + tot) %>%
      # Add highlight and annotation information
      mutate( is_highlight=ifelse(wilcox.p<0.05, "yes", "no")) %>%
      mutate( is_annotate=ifelse(wilcox.p<0.05 , "yes", "no"))       
    
    p<- ggplot(plotdf ,  aes(x=x, y=-sign(t)*log10(wilcox.p))) +
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
      ggtitle("Plasma - GF vs. SPF" )+  xlab("") +
      # Custom the theme:
      theme_bw()+
      theme( 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text( face = "bold", angle = 90, vjust = 0.5, hjust=1),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        legend.position = "none")  
    grid.arrange(p)
    grid.arrange(p+geom_label_repel( data=subset(plotdf, is_annotate=="yes"), aes(label=Metabolite), size=2) )
    dev.off()
    
    pdf(file="./figures/plasma_manhattan_microb_conditional_Jan21.pdf" ,onefile=T  , height = 16, width = 14)
    microbiota.res <- read.xlsx("./results/caltech.PD mice.011920.xlsx","conditional.microbita.effect") %>%
      subset(Tissue =="Plasma" & imputation == "imputed" )
    
    axisdf <- data.frame (Metabolite= mets.included$Plasma , mets.metadata[mets.included$Plasma,c("Analyte.Class","Metabolite.Name")]) %>%
      # Compute number of analyte for each class
      group_by(Analyte.Class) %>% 
      dplyr::summarise(n=n()) %>% 
      # Calculate cumulative position of each analyte class
      mutate(gap=ifelse(n>6,0,3),tot=cumsum(n)-n+cumsum(gap),  start = tot+1, end=tot+n , center = round((start+end)/2), col=c(rep(c("grey", "white"), 12 ),"grey")) %>%
      dplyr::select(-gap)
    
    plotdf <- data.frame (Metabolite= mets.included$Plasma , mets.metadata[mets.included$Plasma,c("Analyte.Class","Metabolite.Name")]) %>%
      ddply(., .(Analyte.Class), function(x){  
        x$MetNum = 1:nrow(x)
        return (x)
        #return (data.frame(Metabolite = x$Metabolite, MetNum = 1:nrow(x)))
      }) %>% 
      right_join(., 
                microbiota.res %>% dplyr::select(-Analyte.Class,-Metabolite.Name)
                , by = "Metabolite") %>%
      left_join(., axisdf , by = "Analyte.Class") %>%
      mutate(x = MetNum + tot) %>%
      # Add highlight and annotation information
      mutate( is_highlight=ifelse(wilcox.p<0.05, "yes", "no")) %>%
      mutate( is_annotate=ifelse(wilcox.p<0.05 , "yes", "no"))       
    
    p<- ggplot(plotdf ,  aes(x=x, y=-sign(t)*log10(wilcox.p))) +
      facet_wrap(~Genotype, ncol = 1 ,  scales='free')+
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
      ggtitle("Plasma - GF vs. SPF" )+  xlab("") +
      # Custom the theme:
      theme_bw()+
      theme( 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text( face = "bold", angle = 90, vjust = 0.5, hjust=1),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        legend.position = "none")  
    grid.arrange(p)
    grid.arrange(p+geom_label_repel( data=subset(plotdf, is_annotate=="yes"), aes(label=Metabolite), size=2) )
    dev.off()
    
    ##############################################################################################################
    ###############    ###############    Plot Genotype Effect    ###############    ###############    ########
    ##############################################################################################################
    pdf(file="./figures/plasma_manhattan_genotype_marginal_Jan21.pdf" ,onefile=T  , height = 8, width = 14)
    
    genotype.res <- read.xlsx("./results/caltech.PD mice.011920.xlsx","marginal.genotype.effect") %>%
      subset(Tissue =="Plasma" & imputation == "imputed" )
    
    axisdf <- data.frame (Metabolite= mets.included$Plasma , mets.metadata[mets.included$Plasma,c("Analyte.Class","Metabolite.Name")]) %>%
      # Compute number of analyte for each class
      group_by(Analyte.Class) %>% 
      dplyr::summarise(n=n()) %>% 
      # Calculate cumulative position of each analyte class
      mutate(gap=ifelse(n>6,0,3),tot=cumsum(n)-n+cumsum(gap),  start = tot+1, end=tot+n , center = round((start+end)/2), col=c(rep(c("grey", "white"), 12 ),"grey")) %>%
      dplyr::select(-gap)
    
    plotdf <- data.frame (Metabolite= mets.included$Plasma , mets.metadata[mets.included$Plasma,c("Analyte.Class","Metabolite.Name")]) %>%
      ddply(., .(Analyte.Class), function(x){  
        x$MetNum = 1:nrow(x)
        return (x)
        #return (data.frame(Metabolite = x$Metabolite, MetNum = 1:nrow(x)))
      }) %>% 
      right_join(., 
                 genotype.res %>% dplyr::select(-Analyte.Class,-Metabolite.Name)
                 , by = "Metabolite") %>%
      left_join(., axisdf , by = "Analyte.Class") %>%
      mutate(x = MetNum + tot) %>%
      # Add highlight and annotation information
      mutate( is_highlight=ifelse(wilcox.p<0.05, "yes", "no")) %>%
      mutate( is_annotate=ifelse(wilcox.p<0.01 , "yes", "no"))       
    
    p<- ggplot(plotdf ,  aes(x=x, y=-sign(t)*log10(wilcox.p))) +
      # Show all points
      #geom_point( aes(color=as.factor(Analyte.Class)), alpha=0.8, size=1.3) +
      geom_point(size=1.3)+
      # custom X axis:
      scale_x_continuous( label = axisdf$Analyte.Class, breaks= axisdf$center ) +
      #add rectangles to distinguish classes
      geom_rect(data = axisdf,
                aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill=col),
                alpha = 0.5  , inherit.aes = FALSE)+    
      scale_fill_manual(values = alpha(c("grey", "white"), .3))+
      #add line to show significance
      geom_hline(yintercept=c(log10(0.05), -log10(0.05)), linetype="dashed", color = "red")+

      # Add highlighted points
      geom_point(data=subset(plotdf, is_highlight=="yes"), color="orange", size=2) +
      
      # Add label using ggrepel to avoid overlapping
      #geom_label_repel( data=subset(plotdf, is_annotate=="yes"), aes(label=Metabolite), size=2) +
      ggtitle("Plasma - ASO vs WT" )+  xlab("") +
      # Custom the theme:
      theme_bw()+
      theme( 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text( face = "bold", angle = 90, vjust = 0.5, hjust=1),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        legend.position = "none")  
    grid.arrange(p)
    
    grid.arrange(p+geom_label_repel( data=subset(plotdf, is_annotate=="yes"), aes(label=Metabolite), size=2) )
    
    dev.off()
    pdf(file="./figures/plasma_manhattan_genotype_conditional_Jan21.pdf" ,onefile=T  , height = 16, width = 14)
    
    genotype.res <- read.xlsx("./results/caltech.PD mice.011920.xlsx","conditional.genotype.effect") %>%
      subset(Tissue =="Plasma" & imputation == "imputed" )
    
    axisdf <- data.frame (Metabolite= mets.included$Plasma , mets.metadata[mets.included$Plasma,c("Analyte.Class","Metabolite.Name")]) %>%
      # Compute number of analyte for each class
      group_by(Analyte.Class) %>% 
      dplyr::summarise(n=n()) %>% 
      # Calculate cumulative position of each analyte class
      mutate(gap=ifelse(n>6,0,3),tot=cumsum(n)-n+cumsum(gap),  start = tot+1, end=tot+n , center = round((start+end)/2), col=c(rep(c("grey", "white"), 12 ),"grey")) %>%
      dplyr::select(-gap)
    
    plotdf <- data.frame (Metabolite= mets.included$Plasma , mets.metadata[mets.included$Plasma,c("Analyte.Class","Metabolite.Name")]) %>%
      ddply(., .(Analyte.Class), function(x){  
        x$MetNum = 1:nrow(x)
        return (x)
        #return (data.frame(Metabolite = x$Metabolite, MetNum = 1:nrow(x)))
      }) %>% 
      right_join(., 
                 genotype.res %>% dplyr::select(-Analyte.Class,-Metabolite.Name)
                 , by = "Metabolite") %>%
      left_join(., axisdf , by = "Analyte.Class") %>%
      mutate(x = MetNum + tot) %>%
      # Add highlight and annotation information
      mutate( is_highlight=ifelse(wilcox.p<0.05, "yes", "no")) %>%
      mutate( is_annotate=ifelse(wilcox.p<0.01 , "yes", "no"))       
    
    p<- ggplot(plotdf ,  aes(x=x, y=-sign(t)*log10(wilcox.p))) +
      facet_wrap(~Microbiota, ncol = 1 ,  scales='free')+
      # Show all points
      #geom_point( aes(color=as.factor(Analyte.Class)), alpha=0.8, size=1.3) +
      geom_point(size=1.3)+
      # custom X axis:
      scale_x_continuous( label = axisdf$Analyte.Class, breaks= axisdf$center ) +
      #add rectangles to distinguish classes
      geom_rect(data = axisdf,
                aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill=col),
                alpha = 0.5  , inherit.aes = FALSE)+    
      scale_fill_manual(values = alpha(c("grey", "white"), .3))+
      #add line to show significance
      geom_hline(yintercept=c(log10(0.05), -log10(0.05)), linetype="dashed", color = "red")+
      
      # Add highlighted points
      geom_point(data=subset(plotdf, is_highlight=="yes"), color="orange", size=2) +
      
      # Add label using ggrepel to avoid overlapping
      #geom_label_repel( data=subset(plotdf, is_annotate=="yes"), aes(label=Metabolite), size=2) +
      ggtitle("Plasma - ASO vs WT" )+  xlab("") +
      # Custom the theme:
      theme_bw()+
      theme( 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text( face = "bold", angle = 90, vjust = 0.5, hjust=1),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        legend.position = "none")  
    grid.arrange(p)
    
    grid.arrange(p+geom_label_repel( data=subset(plotdf, is_annotate=="yes"), aes(label=Metabolite), size=2) )
    
    dev.off()
    
    #######PLOT BOXPLOTS
    pdf(file="./figures/boxplot_caltech_Jan21.pdf" ,onefile=T  , height = 6, width = 6)
    res_microbiota_marginal <- read.xlsx("./results/caltech.PD mice.011920.xlsx","marginal.microbita.effect") %>%
      subset( imputation == "imputed" )
    res_microbiota_cond <- read.xlsx("./results/caltech.PD mice.011920.xlsx","conditional.microbita.effect") %>%
      subset( imputation == "imputed" )
    res_genotype_marginal <- read.xlsx("./results/caltech.PD mice.011920.xlsx","marginal.genotype.effect") %>%
      subset( imputation == "imputed" )
    res_genotype_cond <- read.xlsx("./results/caltech.PD mice.011920.xlsx","conditional.genotype.effect") %>%
      subset( imputation == "imputed" )
    
    th <- theme_bw()+
      theme( axis.text.x = element_text( face = "bold", angle = 0, vjust = 0.5, hjust=1),
             strip.text = element_text(face = "bold"),
             plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
             legend.position = "bottom")  
    
    ddply(dta.imputed,.(Tissue), function(dat){
      material <- unique(dat$Material)[1]
      tissue <- unique ( dat$Tissue)[1]
      
      for ( met in mets.included[[material]]){
        print(paste0(material , ":",paste0(tissue, collapse = ";"), ":",met))
        
        tryCatch({
          print(paste0(material , ":",paste0(tissue, collapse = ";"), ":",met))
          plottitle = tissue
          dat$Metabolite = dat[,met]
          
          temp = ifelse(material == "Plasma" ,"plasma", ifelse(length(tissue)>1,"brain",paste0("brain.",tolower(tissue))))
          
          p_ASO_vs_WT <- (res_genotype_marginal %>% subset(Metabolite==met & Tissue == tissue & imputation=="imputed" ))["wilcox.p"]
          p_GF_vs_SPF <- (res_microbiota_marginal %>% subset(Metabolite==met & Tissue == tissue & imputation=="imputed" ))["wilcox.p"]
          
          p_ASO_GF_vs_SPF <- (res_microbiota_cond %>% subset(Metabolite==met & Tissue == tissue & imputation=="imputed" & Genotype=="ASO"))["wilcox.p"]
          p_WT_GF_vs_SPF <- (res_microbiota_cond %>% subset(Metabolite==met & Tissue == tissue & imputation=="imputed" & Genotype=="WT"))["wilcox.p"]
          
          p_GF_ASO_vs_WT <- (res_genotype_cond %>% subset(Metabolite==met & Tissue == tissue & imputation=="imputed" & Microbiota=="GF"))["wilcox.p"]
          p_SPF_ASO_vs_WT <- (res_genotype_cond %>% subset(Metabolite==met & Tissue == tissue & imputation=="imputed" & Microbiota=="SPF"))["wilcox.p"]
          
          annot_tex = sprintf("Wilcoxon rank-sum test:\nGF vs SPF:%0.2E    ASO vs WT:%0.2E\nASO(GF vs SPF):%0.2E    WT(GF vs SPF):%0.2E\nGF(ASO vs WT):%0.2E   SPF(ASO vs WT):%0.2E",p_GF_vs_SPF,p_ASO_vs_WT, p_ASO_GF_vs_SPF,p_WT_GF_vs_SPF,p_GF_ASO_vs_WT,p_SPF_ASO_vs_WT)
          grob <- grid::grobTree(textGrob(annot_tex, x=0.1,  y=0.9, hjust=0,
                                          gp=gpar(col="red", fontsize=10, fontface="italic")))
          p1 =  ggviolin (dat, x = "Genotype", y = "Metabolite",
                          color = "Microbiota", 
                          add = "jitter", shape = "Microbiota",
                          draw_quantiles = c(0.25,0.5,0.75))+
            ggtitle(sprintf("%s Genotype x Microbiota\n%s",met , plottitle ))+ylab(expression("Concentration("* mu*"M)") ) +
            annotation_custom(grob) +th
          
          grid.arrange(p1)
          
        },
        error = function(e){
          message(paste(material,":",paste0(tissue, collapse = ":"),"----",met,"\n"), e)})
      }
      return (1)
      
    })
    
    dev.off()

    ####BOX PLOT FIRST LOOKED AT THE MICROBIOTA
    pdf(file="./figures/boxplot_caltech_Jan21_2.pdf" ,onefile=T  , height = 6, width = 6)
    res_microbiota_marginal <- read.xlsx("./results/caltech.PD mice.011920.xlsx","marginal.microbita.effect") %>%
      subset( imputation == "imputed" )
    res_microbiota_cond <- read.xlsx("./results/caltech.PD mice.011920.xlsx","conditional.microbita.effect") %>%
      subset( imputation == "imputed" )
    res_genotype_marginal <- read.xlsx("./results/caltech.PD mice.011920.xlsx","marginal.genotype.effect") %>%
      subset( imputation == "imputed" )
    res_genotype_cond <- read.xlsx("./results/caltech.PD mice.011920.xlsx","conditional.genotype.effect") %>%
      subset( imputation == "imputed" )
    
    th <- theme_bw()+
      theme( axis.text.x = element_text( face = "bold", angle = 0, vjust = 0.5, hjust=1),
             strip.text = element_text(face = "bold"),
             plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
             legend.position = "bottom")  
    
    ddply(dta.imputed,.(Tissue), function(dat){
      material <- unique(dat$Material)[1]
      tissue <- unique ( dat$Tissue)[1]
      
      for ( met in mets.included[[material]]){
        print(paste0(material , ":",paste0(tissue, collapse = ";"), ":",met))
        
        tryCatch({
          print(paste0(material , ":",paste0(tissue, collapse = ";"), ":",met))
          plottitle = tissue
          dat$Metabolite = dat[,met]
          
          temp = ifelse(material == "Plasma" ,"plasma", ifelse(length(tissue)>1,"brain",paste0("brain.",tolower(tissue))))
          
          p_ASO_vs_WT <- (res_genotype_marginal %>% subset(Metabolite==met & Tissue == tissue & imputation=="imputed" ))["wilcox.p"]
          p_GF_vs_SPF <- (res_microbiota_marginal %>% subset(Metabolite==met & Tissue == tissue & imputation=="imputed" ))["wilcox.p"]
          
          p_ASO_GF_vs_SPF <- (res_microbiota_cond %>% subset(Metabolite==met & Tissue == tissue & imputation=="imputed" & Genotype=="ASO"))["wilcox.p"]
          p_WT_GF_vs_SPF <- (res_microbiota_cond %>% subset(Metabolite==met & Tissue == tissue & imputation=="imputed" & Genotype=="WT"))["wilcox.p"]
          
          p_GF_ASO_vs_WT <- (res_genotype_cond %>% subset(Metabolite==met & Tissue == tissue & imputation=="imputed" & Microbiota=="GF"))["wilcox.p"]
          p_SPF_ASO_vs_WT <- (res_genotype_cond %>% subset(Metabolite==met & Tissue == tissue & imputation=="imputed" & Microbiota=="SPF"))["wilcox.p"]
          
          annot_tex = sprintf("Wilcoxon rank-sum test:\nGF vs SPF:%0.2E    ASO vs WT:%0.2E\nASO(GF vs SPF):%0.2E    WT(GF vs SPF):%0.2E\nGF(ASO vs WT):%0.2E   SPF(ASO vs WT):%0.2E",p_GF_vs_SPF,p_ASO_vs_WT, p_ASO_GF_vs_SPF,p_WT_GF_vs_SPF,p_GF_ASO_vs_WT,p_SPF_ASO_vs_WT)
          grob <- grid::grobTree(textGrob(annot_tex, x=0.1,  y=0.9, hjust=0,
                                          gp=gpar(col="red", fontsize=10, fontface="italic")))
          p1 =  ggviolin (dat, x = "Microbiota", y = "Metabolite",
                          color = "Genotype", 
                          add = "jitter", shape = "Microbiota",
                          draw_quantiles = c(0.25,0.5,0.75))+
            ggtitle(sprintf("%s Genotype x Microbiota\n%s",met , plottitle ))+ylab(expression("Concentration("* mu*"M)") ) +
            annotation_custom(grob) +th
          
          grid.arrange(p1)
          
        },
        error = function(e){
          message(paste(material,":",paste0(tissue, collapse = ":"),"----",met,"\n"), e)})
      }
      return (1)
      
    })
    
    dev.off()
    
##################################    
    ################################################################################################################
    ################################################################################################################
    ######April
    #################Marginal-   For Each Tissue , Plot GF vs SPF ; ASO vs WT
    #################Conditional For Each Tissue , Plot ASO vs WT/ SPF, ASO vs WT/ GF, GF vs SPF/ ASO, GF vs SPF/ WT
    ################################################################################################################
    pdf(file="./figures/April 2021/manhattan_marginal_April21.pdf" ,onefile=T  , height = 16, width = 14)
    
    microbiota.res <- read.xlsx("./results/caltech.PD mice.011920.xlsx","marginal.microbita.effect") %>%
      subset(imputation == "imputed" )
    genotype.res <- read.xlsx("./results/caltech.PD mice.011920.xlsx","marginal.genotype.effect") %>%
      subset(imputation == "imputed" )
    ldply(microbiota.res$Tissue %>% unique %>% as.list %>% setNames(.,microbiota.res$Tissue %>% unique) , function(x){
      print (x)
      if ( x%in% c("Brainstem","Cortex","Nigra","Striatum")){
        material ="Brain Tissue"
      }else if ( x%in% c("Cecum","Colon Content","Duodenum Content")){
        material ="Gut Content"
      }else if ( x%in% c("Colon","Duodenum")){
        material ="Gut Content"
      }else{
        material <- "Plasma"
      }
      
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
                    genotype.res   %>% subset(Tissue == x) %>% dplyr::select(Metabolite , t ,wilcox.p) %>% mutate(Effect="Genotype (ASO vs. WT)")
                   )
                   , by = "Metabolite") %>%
        left_join(., axisdf , by = "Analyte.Class") %>%
        mutate(x = MetNum + tot) %>%
        # Add highlight and annotation information
        mutate( is_highlight=ifelse(wilcox.p<0.05, "yes", "no")) %>%
        mutate( is_annotate=ifelse(wilcox.p<0.05 , "yes", "no"))       
      
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
        geom_label_repel( data=subset(plotdf, is_annotate=="yes"), aes(label=Metabolite), size=2)+
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
    
    pdf(file="./figures/April 2021/manhattan_conditional_April21.pdf" ,onefile=T  , height = 16, width = 20)
    
    microbiota.res <- read.xlsx("./results/caltech.PD mice.011920.xlsx","conditional.microbita.effect") %>%
      subset(imputation == "imputed" )
    genotype.res <- read.xlsx("./results/caltech.PD mice.011920.xlsx","conditional.genotype.effect") %>%
      subset(imputation == "imputed" )
    ldply(microbiota.res$Tissue %>% unique %>% as.list %>% setNames(.,microbiota.res$Tissue %>% unique) , function(x){
      print (x)
      if ( x%in% c("Brainstem","Cortex","Nigra","Striatum")){
        material ="Brain Tissue"
      }else if ( x%in% c("Cecum","Colon Content","Duodenum Content")){
        material ="Gut Content"
      }else if ( x%in% c("Colon","Duodenum")){
        material ="Gut Content"
      }else{
        material <- "Plasma"
      }
      
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
                     microbiota.res %>% subset(Tissue == x) %>% dplyr::select(Metabolite , Genotype , t ,wilcox.p) %>% mutate(Effect=paste0("GF vs. SPF|",Genotype)),
                     genotype.res   %>% subset(Tissue == x) %>% dplyr::select(Metabolite , Microbiota, t ,wilcox.p) %>% mutate(Effect=paste0("ASO vs. WT|",Microbiota))
                   )
                   , by = "Metabolite") %>%
        left_join(., axisdf , by = "Analyte.Class") %>%
        mutate(x = MetNum + tot) %>%
        # Add highlight and annotation information
        mutate( is_highlight=ifelse(wilcox.p<0.05, "yes", "no")) %>%
        mutate( is_annotate=ifelse(wilcox.p<0.05 , "yes", "no"))       
      
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
        geom_label_repel( data=subset(plotdf, is_annotate=="yes"), aes(label=Metabolite), size=2)+
        
        
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
    ############################################################################################################################################
    ############################################################################################################################################
    ############For each combination of tissue, geno, micro compute average analyte class then correlate between tissue correlations
    #average.class.concentration <- ddply(dta.imputed %>% subset(Material %in% c("Brain Tissue","Plasma","Intestinal Tract Tissue")), .(Tissue ) , function(dat){
    average.class.concentration <- ddply(dta.imputed , .(Tissue ) , function(dat){
      print(dat$Tissue[1])
      df <- data.frame(
        Metabolite = mets.included[[dat$Material[1]]],
        stringsAsFactors = F
      ) %>% dplyr::mutate(class = mets.metadata[Metabolite,]$Analyte.Class )
      class.mean <- ddply(df ,   .(class), function(xx){
        y = dat %>% dplyr::select(Animal ,Genotype,Microbiota, all_of(xx$Metabolite)) 
        y$mean <- apply(y[-(1:3)] %>% scale, 1 , mean ,na.rm=T)
        return(y %>% dplyr::select(Animal,Genotype,Microbiota , mean))
      })
      return( dcast( class.mean , Animal+Genotype+Microbiota ~ class, value.var = "mean" ) )
    })
    dat <- average.class.concentration %>% reshape2::melt(.,id.vars = c("Tissue","Animal","Genotype","Microbiota")) %>% 
      dcast(., Animal+Genotype+Microbiota ~ Tissue+variable , value.var = "value" )
    get_upper_tri <- function(cormat){
      cormat[lower.tri(cormat)]<- NA
      return(cormat)
    }
    plot_corr <- function(dat ){
      grid.arrange( 
        ggplot() +# Draw ggplot2 plot with text only
          annotate("text",x = 1,y = 1,size = 8,label = dat$facetingtext[1]) + theme_void()
      )
      
      cormat.tissues <- unique(dta$Tissue)
      mypalette <- ggpubr::get_palette("jco",length(cormat.tissues)) %>% setNames(., cormat.tissues )
      
      slicing <- c("All Samples","ASO","WT","GF","SPF","ASO GF","ASO SPF","WT GF","WT SPF")
      mypalette_splicing <- ggpubr::get_palette("jco",length(slicing)) %>% setNames(., slicing )
      
      cor.metabs_r <- ddply(dat , .(faceting) , function(xx){
        dd = psych::corr.test(xx[-(1:5)] , method="spearman" , ci=F)$r
        dd <- data.frame(metab = rownames(dd) , stringsAsFactors = F) %>% cbind(dd )
        return(dd )
      })
      cor.metabs_p <- ddply(dat , .(faceting) , function(xx){
        dd = psych::corr.test(xx[-(1:5)] , method="spearman" , ci=F)$p
        dd = apply(dd, 2, function(xxx){
          xxx = cut (xxx , c(0,.001,.01,.05 , 1) , labels = c("***","**","*","") , include.lowest = T)
          xxx[is.na(xxx)]=""
          
          return(xxx)
        })
        dd <- data.frame(metab = colnames(dd) , stringsAsFactors = F ) %>% cbind(dd )
      })
      
      #Plot correlation of Plasma with other tissues
      for ( tissue in cormat.tissues){
        ncol=1
        nrow=1
        if (length(unique(dat$faceting))==2){
          ncol=2
        }else if (length(unique(dat$faceting))==4){
          nrow=2
          ncol=2
        }
        i=1
        j=1
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nr = nrow, nc =ncol)))
        
        for (facet in unique(dat$faceting)){
          r <- cor.metabs_r  %>%
            subset(grepl("^Plasma",metab) & faceting==facet) %>% dplyr::select(metab, grep(paste0(tissue,"_"), colnames(cor.metabs_r)))
          p <- cor.metabs_p  %>%
            subset(grepl("^Plasma",metab) & faceting==facet) %>% dplyr::select(metab, grep(paste0(tissue,"_"), colnames(cor.metabs_r)))
          class(r) = class(p) = "data.frame"
          
          #Define column and wo splits
          if (dat$facetingtext[1]=="All Samples"){
            csplit_main = rep("All Samples" , ncol(r)-1) %>% factor
          }else if (dat$facetingtext[1]=="Genotype"){
            csplit_main = rep(facet , ncol(r)-1) %>% factor
          }else if (dat$facetingtext[1]=="Microbiota"){
            csplit_main = rep(facet , ncol(r)-1) %>% factor
          }else if (dat$facetingtext[1]=="G_x_M"){
            csplit_main = rep(facet , ncol(r)-1) %>% factor
          }
          
          csplit  = ldply(strsplit(colnames(r)[-1],"_") , function(x){return(x[[1]])}) %>% unlist %>% factor(., levels = cormat.tissues )
          rsplit  = rep("Plasma", nrow(r)) %>% factor(., levels = cormat.tissues )
          
          colnames(r)[-1] <- ldply(strsplit(colnames(r)[-1],"_") , function(x){return(x[[2]])}) %>% unlist
          rownames(r) <- ldply(strsplit(r$metab,"_") , function(x){return(x[[2]])}) %>% unlist
          r= r[-1]
          p = p[-1]
          
          if (tissue == "Plasma"){
            r <- r %>% get_upper_tri
            p <- p %>% get_upper_tri
            p[is.na(p)]=""
          }
          topanno = HeatmapAnnotation(
            Tissue = anno_block(
              gp = gpar(fill = mypalette[tissue]),
              labels = tissue, 
              labels_gp = gpar(col = "white", fontsize = 10)
            ),
            Split = anno_block(
              gp = gpar(fill = mypalette_splicing[as.character(csplit_main)]),
              labels = facet, 
              labels_gp = gpar(col = "white", fontsize = 10)
            )
          )
          leftanno = rowAnnotation(
            Tissue = anno_block(
              gp = gpar(fill = mypalette["Plasma"]),
              labels = "Plasma", 
              labels_gp = gpar(col = "white", fontsize = 10)
            )
          )
          
          ht<-Heatmap(r ,
                      rect_gp = gpar(col = "white"),
                      col=colorRamp2(c(-1,0,1), c("red" , "white", "blue")), na_col = "lightgrey",
                      cell_fun = function(i, j, x, y, width, height, fill) {
                        grid.text(p[j, i], x = x, y = y)}, 
                      cluster_columns = F, cluster_rows  = F,
                      column_split  = csplit_main,
                      row_split = rsplit,
                      left_annotation = leftanno,
                      top_annotation = topanno,
                      column_title_side = "top",
                      #width = unit(20, "cm"),
                      show_heatmap_legend = F,
                      name="Spearman r")
          print(c(i,j))
          pushViewport(viewport(layout.pos.row = i, layout.pos.col = j))
          draw(ht, newpage = FALSE)
          upViewport()
          j = j+1
          if ( j>ncol){
            j=1
            i=i+1
          }
        }
      }
    
  }
  
  m <- colnames(dat)[-c(1:3)]
  dat_cor <- rbind(
    dat %>% mutate(facetingtext = "All Samples" , faceting = "All Samples") %>% dplyr::select(Animal,Genotype,Microbiota ,facetingtext, faceting  , all_of(m)) ,
    dat %>% mutate(facetingtext = "Genotype", faceting =Genotype) %>% dplyr::select(Animal,Genotype,Microbiota ,facetingtext, faceting  , all_of(m)) ,
    dat %>% mutate(facetingtext = "Microbiota", faceting =Microbiota) %>% dplyr::select(Animal,Genotype,Microbiota ,facetingtext, faceting  , all_of(m)) ,
    dat %>% mutate(facetingtext = "G_x_M" , faceting =paste(Genotype,Microbiota)) %>% dplyr::select(Animal,Genotype,Microbiota ,facetingtext, faceting  , all_of(m)) 
  )
  library(Cairo)
  CairoPDF(file="./figures/April 2021/cross_tissue_cor_all_samples_050621.pdf" ,onefile=T , height=8 ,width=8)
  plot_corr(dat_cor %>% subset(facetingtext == "All Samples"))
  dev.off()
  CairoPDF(file="./figures/April 2021/cross_tissue_cor_G_x_M_050621.pdf" ,onefile=T , height=16 ,width=16)
  plot_corr(dat_cor %>% subset(facetingtext == "G_x_M"))
  dev.off()
  
  CairoPDF(file="./figures/April 2021/cross_tissue_cor_Geno Micro_050621.pdf" ,onefile=T , height=8 ,width=16)
  ddply(dat_cor %>% subset(facetingtext %in% c("Genotype","Microbiota")), .(facetingtext), function(x){
    plot_corr(x)
    return(1)
  })
  dev.off()
  
  anal.dta <- read.csv("caltechdata.csv" , stringsAsFactors = F , check.names = F)
  mets <- colnames(anal.dta) [10:643]
  corr <- dlply(anal.dta , .(imputation), function(X){
    print(unique(X$imputation))
    dat <-dcast(X %>% as.data.table,Animal~Tissue , value.var = mets) %>% dplyr::select(-1)
    return( psych::corr.test(dat , method="spearman" , ci=F) )
  })
  wb <- createWorkbook()
  
  addWorksheet(wb,"imputed data - Spearman r")
  writeDataTable(wb, "imputed data - Spearman r" ,corr$imputed$r %>% data.frame (., check.names = F) , rowNames = T)
  
  addWorksheet(wb,"imputed data - Spearman r-p")
  writeDataTable(wb, "imputed data - Spearman r-p" ,corr$imputed$p%>% data.frame (., check.names = F), rowNames = T)
  
  addWorksheet(wb,"unimputed data - Spearman r")
  writeDataTable(wb, "unimputed data - Spearman r" ,corr$unimputed$r%>% data.frame (., check.names = F), rowNames = T)
  
  addWorksheet(wb,"unimputed data - Spearman r-p")
  writeDataTable(wb, "unimputed data - Spearman r-p" ,corr$unimputed$p%>% data.frame (., check.names = F), rowNames = T)
  saveWorkbook(wb, "./results/caltech PD mice model Cor.060921.xlsx", overwrite = TRUE) #added hampsy, hamsa, hamsom to analysis
  
  #######For Blood select metabolites that are significant with 
  if (FALSE){
    regression.res$plasma.unimputed %>% subset( `Pr(>|t|)` <0.2 & Variable=="MicrobiotaSPF") %>% dplyr::select(Analyte.Class) %>% table
    # Acylcarnitines       Amino Acid Related              Amino Acids          Biogenic Amines         Carboxylic Acids 
    # 1                       10                       10                        2                        1 
    # Ceramides       Cholesteryl Esters             Diglycerides       Dihexosylceramides              Fatty Acids 
    # 7                        8                        4                        4                        3 
    # Hexosylceramides  Indoles and Derivatives Lysophosphatidylcholines     Phosphatidylcholines           Sphingomyelins 
    # 8                        2                        3                       15                        6 
    # Triglycerides 
    # 2   
    mets2tsne.microbiota = regression.res$plasma.unimputed %>% subset( `Pr(>|t|)` <0.2 & Variable=="MicrobiotaSPF") %>% dplyr::select(Metabolite) %>% unlist
    pdf(file="./figures/tsne_blood_oct_10.pdf" ,onefile=T  , height = 6, width = 6.5)
    #Run tSNE on all samples 
    set.seed(40)  
    tsne_caltech_mice <-  tsne(dta%>% subset(Material=="Plasma" & Genotype=="WT") %>% dplyr::select(mets2tsne.microbiota) %>% impute.knn.obs.sel(.)  %>% scale , perplexity=50)
    tsne_caltech_mice <- data.frame(tsne_caltech_mice)
    colnames(tsne_caltech_mice) <- c("t-SNE1","t-SNE2")
    tsne_caltech_mice$Microbiota <- dta%>% subset(Material=="Plasma" & Genotype=="WT") %>% dplyr::select(Microbiota) %>% unlist

    gg <- ggplot(tsne_caltech_mice, aes(x=`t-SNE1`, y=`t-SNE2`,col=Microbiota)) + 
      geom_point( size=4)+ggtitle("Plasma- WT") +
#      scale_x_continuous(limits = c(-18, 18))+
#      scale_y_continuous(limits = c(-18, 18))+
      theme_bw()+
      theme( axis.text.x = element_text( face = "bold"),
             strip.text = element_text(face = "bold"),
             plot.title = element_text(hjust = 0.5, size=20),
             legend.position="bottom")  
    grid.arrange(gg)
    
    set.seed(40)  
    tsne_caltech_mice <-  tsne(dta%>% subset(Material=="Plasma" & Genotype=="ASO") %>% dplyr::select(mets2tsne.microbiota) %>% impute.knn.obs.sel(.)  %>% scale , perplexity=50)
    tsne_caltech_mice <- data.frame(tsne_caltech_mice)
    colnames(tsne_caltech_mice) <- c("t-SNE1","t-SNE2")
    tsne_caltech_mice$Microbiota <- dta%>% subset(Material=="Plasma" & Genotype=="ASO") %>% dplyr::select(Microbiota) %>% unlist
    
    gg <- ggplot(tsne_caltech_mice, aes(x=`t-SNE1`, y=`t-SNE2`,col=Microbiota)) + 
      geom_point( size=4)+ggtitle("Plasma- ASO") +
      #      scale_x_continuous(limits = c(-18, 18))+
      #      scale_y_continuous(limits = c(-18, 18))+
      theme_bw()+
      theme( axis.text.x = element_text( face = "bold"),
             strip.text = element_text(face = "bold"),
             plot.title = element_text(hjust = 0.5, size=20),
             legend.position="bottom")  
    grid.arrange(gg)
    
    
    dev.off()
  }
  
  microb.res_marginal <- read.xlsx("./results/caltech.PD mice.011920.xlsx","marginal.microbita.effect") %>%
    subset(Tissue =="Plasma"  ) %>%
    as.data.table %>% 
    data.table::dcast(.,Metabolite+Analyte.Class~imputation , value.var= c("ttest.estimate","t","ttest.df","ttest.p","wilcox.p","wilcox.W","fisher.pitman.p","fisher.pitman.Z") ) %>%
    mutate(wilcox.q_imputed = p.adjust(wilcox.p_imputed , "BH"), wilcox.q_unimputed = p.adjust(wilcox.p_unimputed , "BH")) %>% 
    subset(wilcox.p_imputed<0.05 | wilcox.p_unimputed <0.05) %>%
    mutate(name_class = paste0(Metabolite,"(",Analyte.Class,")"), 
           t_p_q_imputed=sprintf("t=%0.2f\np=%0.2e(q=%0.02e)",t_imputed, wilcox.p_imputed,wilcox.q_imputed),
           t_p_q_unimputed=sprintf("t=%0.2f\np=%0.2e(q=%0.02e)",t_unimputed, wilcox.p_unimputed,wilcox.q_unimputed))
  
  ord = order(microb.res_marginal$wilcox.p_imputed)
  microb.res_marginal <- microb.res_marginal[ord,]
  
  mets.included <- read.xlsx("./results/caltech.PD mice.011920.xlsx","conditional.microbita.effect") %>%
    subset(Tissue =="Plasma"  ) %>% subset(wilcox.p<0.05) %>% dplyr::select(Metabolite) %>% unlist
  microb.res_cond <- read.xlsx("./results/caltech.PD mice.011920.xlsx","conditional.microbita.effect") %>%
    subset(Tissue =="Plasma" & Metabolite %in% mets.included ) %>%
    as.data.table %>% 
    data.table::dcast(.,Metabolite+Analyte.Class~Genotype+imputation , value.var= c("ttest.estimate","t","ttest.df","ttest.p","wilcox.p","wilcox.W","fisher.pitman.p","fisher.pitman.Z") ) %>%
    mutate(wilcox.q_ASO_imputed = p.adjust(wilcox.p_ASO_imputed , "BH"), wilcox.q_ASO_unimputed = p.adjust(wilcox.p_ASO_unimputed , "BH"),
           wilcox.q_WT_imputed = p.adjust(wilcox.p_WT_imputed , "BH"), wilcox.q_WT_unimputed = p.adjust(wilcox.p_WT_unimputed , "BH")) %>% 
    mutate(name_class = paste0(Metabolite,"(",Analyte.Class,")"), 
           t_p_q_ASO.imputed=sprintf("t=%0.2f\np=%0.2e(q=%0.02e)",t_ASO_imputed, wilcox.p_ASO_imputed,wilcox.q_ASO_imputed),t_p_q_ASO_unimputed=sprintf("t=%0.2f\np=%0.2e(q=%0.02e)",t_ASO_unimputed, wilcox.p_ASO_unimputed,wilcox.q_ASO_unimputed),
           t_p_q_WT.imputed=sprintf("t=%0.2f\np=%0.2e(q=%0.02e)",t_WT_imputed, wilcox.p_WT_imputed,wilcox.q_WT_imputed),t_p_q_WT_unimputed=sprintf("t=%0.2f\np=%0.2e(q=%0.02e)",t_WT_unimputed, wilcox.p_WT_unimputed,wilcox.q_WT_unimputed))
  ord = order(microb.res_cond$wilcox.p_WT_imputed)
  microb.res_cond <- microb.res_cond[ord,]
  
  
  geno.res_marginal <- read.xlsx("./results/caltech.PD mice.011920.xlsx","marginal.genotype.effect") %>%
    subset(Tissue =="Plasma"  ) %>%
    as.data.table %>% 
    data.table::dcast(.,Metabolite+Analyte.Class~imputation , value.var= c("ttest.estimate","t","ttest.df","ttest.p","wilcox.p","wilcox.W","fisher.pitman.p","fisher.pitman.Z") ) %>%
    mutate(wilcox.q_imputed = p.adjust(wilcox.p_imputed , "BH"), wilcox.q_unimputed = p.adjust(wilcox.p_unimputed , "BH")) %>% 
    subset(wilcox.p_imputed<0.05 | wilcox.p_unimputed <0.05) %>%
    mutate(name_class = paste0(Metabolite,"(",Analyte.Class,")"), 
           t_p_q_imputed=sprintf("t=%0.2f\np=%0.2e(q=%0.02e)",t_imputed, wilcox.p_imputed,wilcox.q_imputed),
           t_p_q_unimputed=sprintf("t=%0.2f\np=%0.2e(q=%0.02e)",t_unimputed, wilcox.p_unimputed,wilcox.q_unimputed))
  
  ord = order(geno.res_marginal$wilcox.p_imputed)
  geno.res_marginal <- geno.res_marginal[ord,]
  
  mets.included <- read.xlsx("./results/caltech.PD mice.011920.xlsx","conditional.genotype.effect") %>%
    subset(Tissue =="Plasma"  ) %>% subset(wilcox.p<0.05) %>% dplyr::select(Metabolite) %>% unlist
  geno.res_cond <- read.xlsx("./results/caltech.PD mice.011920.xlsx","conditional.genotype.effect") %>%
    subset(Tissue =="Plasma" & Metabolite %in% mets.included ) %>%
    as.data.table %>% 
    data.table::dcast(.,Metabolite+Analyte.Class~Microbiota+imputation , value.var= c("ttest.estimate","t","ttest.df","ttest.p","wilcox.p","wilcox.W","fisher.pitman.p","fisher.pitman.Z") ) %>%
    mutate(wilcox.q_GF_imputed = p.adjust(wilcox.p_GF_imputed , "BH"), wilcox.q_GF_unimputed = p.adjust(wilcox.p_GF_unimputed , "BH"),
           wilcox.q_SPF_imputed = p.adjust(wilcox.p_SPF_imputed , "BH"), wilcox.q_SPF_unimputed = p.adjust(wilcox.p_SPF_unimputed , "BH")) %>% 
    mutate(name_class = paste0(Metabolite,"(",Analyte.Class,")"), 
           t_p_q_GF.imputed=sprintf("t=%0.2f\np=%0.2e(q=%0.02e)",t_GF_imputed, wilcox.p_GF_imputed,wilcox.q_GF_imputed),t_p_q_GF_unimputed=sprintf("t=%0.2f\np=%0.2e(q=%0.02e)",t_GF_unimputed, wilcox.p_GF_unimputed,wilcox.q_GF_unimputed),
           t_p_q_SPF.imputed=sprintf("t=%0.2f\np=%0.2e(q=%0.02e)",t_SPF_imputed, wilcox.p_SPF_imputed,wilcox.q_SPF_imputed),t_p_q_SPF_unimputed=sprintf("t=%0.2f\np=%0.2e(q=%0.02e)",t_SPF_unimputed, wilcox.p_SPF_unimputed,wilcox.q_SPF_unimputed))
  ord = order(geno.res_cond$wilcox.p_SPF_imputed)
  geno.res_cond <- geno.res_cond[ord,]
  
  wb <- createWorkbook()
  
  addWorksheet(wb,"microb_marginal")
  writeDataTable(wb, "microb_marginal" ,microb.res_marginal)
  
  addWorksheet(wb,"microb_conditional")
  writeDataTable(wb, "microb_conditional" ,microb.res_cond)
  
  addWorksheet(wb,"genotype_marginal")
  writeDataTable(wb, "genotype_marginal" ,geno.res_marginal)
  
  addWorksheet(wb,"genotype_conditional")
  writeDataTable(wb, "genotype_conditional" ,geno.res_cond)
  
  saveWorkbook(wb, "./results/forppt.011920.xlsx", overwrite = TRUE) #added hampsy, hamsa, hamsom to analysis
  
  
  
  
  L_geno_microb <- list()
  microb.res_cond <- read.xlsx("./results/caltech.PD mice.011920.xlsx","conditional.microbita.effect") 
  L_microb <- dlply(microb.res_cond , .(Tissue), function(x){return (x)})
  geno.res_cond <- read.xlsx("./results/caltech.PD mice.011920.xlsx","conditional.genotype.effect")
  L_geno <- dlply(geno.res_cond , .(Tissue), function(x){return (x)})
  L_geno_microb <- dlply(data.frame(Tissue = names(L_geno) , stringsAsFactors = F)   ,.(Tissue), function(x){
    tissue  = x$Tissue[1]
    
    d1 <- L_microb[[tissue]] %>% subset(imputation == "imputed") %>%
      as.data.table %>% 
      dplyr::rename(t_GF_vs_SPF = t ,wilcox.p_GF_vs_SPF = wilcox.p ) %>% 
      data.table::dcast(.,Metabolite+Analyte.Class~Genotype , value.var= c("t_GF_vs_SPF","wilcox.p_GF_vs_SPF") ) %>%
      mutate(wilcox.q_GF_vs_SPF_ASO = p.adjust(wilcox.p_GF_vs_SPF_ASO , "BH"),
             wilcox.q_GF_vs_SPF_WT = p.adjust(wilcox.p_GF_vs_SPF_WT , "BH")) %>% 
      mutate(name_class = paste0(Metabolite,"(",Analyte.Class,")"), 
             t_p_q_ASO=sprintf("t=%0.2f\np=%0.2e(q=%0.02e)",t_GF_vs_SPF_ASO, wilcox.p_GF_vs_SPF_ASO,wilcox.q_GF_vs_SPF_ASO),
             t_p_q_WT=sprintf("t=%0.2f\np=%0.2e(q=%0.02e)",t_GF_vs_SPF_WT, wilcox.p_GF_vs_SPF_WT,wilcox.q_GF_vs_SPF_WT))
    d2 <- L_geno[[tissue]] %>% subset(imputation == "imputed") %>%
      as.data.table %>% 
      dplyr::rename(t_ASO_vs_WT = t ,wilcox.p_ASO_vs_WT = wilcox.p ) %>% 
      data.table::dcast(.,Metabolite+Analyte.Class~Microbiota , value.var= c("t_ASO_vs_WT","wilcox.p_ASO_vs_WT") ) %>%
      mutate(wilcox.q_ASO_vs_WT_GF = p.adjust(wilcox.p_ASO_vs_WT_GF , "BH"),
             wilcox.q_ASO_vs_WT_SPF = p.adjust(wilcox.p_ASO_vs_WT_SPF , "BH")) %>% 
      mutate(
             t_p_q_GF=sprintf("t=%0.2f\np=%0.2e(q=%0.02e)",t_ASO_vs_WT_GF, wilcox.p_ASO_vs_WT_GF,wilcox.q_ASO_vs_WT_GF),
             t_p_q_SPF=sprintf("t=%0.2f\np=%0.2e(q=%0.02e)",t_ASO_vs_WT_SPF, wilcox.p_ASO_vs_WT_SPF,wilcox.q_ASO_vs_WT_SPF))
    
    merged <- full_join (
      d1,
      d2 %>% dplyr::select(-Analyte.Class),
      by = "Metabolite"
      ) %>% mutate(imputation = "imputed")
    return (merged)
  })
  
  ###Save New formatted Results
  temp = lapply(1:length(L_geno_microb) , function(x){
    print(x)
    return(1)
  })
  
  wb <- createWorkbook()
  temp = lapply(1:length(L_geno_microb) , function(x){
    addWorksheet(wb,names(L_geno_microb)[x])
    writeDataTable(wb, names(L_geno_microb)[x] ,L_geno_microb[[x]])
    return(1)
  })
  saveWorkbook(wb, "./results/caltech.PD mice_G_by_M_results.xlsx", overwrite = TRUE) #added hampsy, hamsa, hamsom to analysis
  
  with( anal.dta %>% subset(imputation =="imputed" &Tissue=="Plasma") , 
        interaction.plot(Microbiota, Genotype , log(`TG(16:1_34:2)`) , fixed=T)
  )
  with( anal.dta %>% subset(imputation =="imputed" &Tissue=="Plasma") , 
        lm(log(`TG(16:1_34:2)`)~Microbiota*Genotype ) %>% anova
  )
  L_geno_microb$Plasma %>% subset(wilcox.p_ASO_vs_WT_GF <0.05) %>% dplyr::select(Metabolite) %>% unlist
  
  LOD.by.genotype <- ddply(dta , .(Material) , function(x){
    percent.b.LOD.SPF <- apply(x %>% subset(Microbiota=="SPF") %>% dplyr::select(mets), 2, function(xx){
      mean(!is.na(xx))
    })
    percent.b.LOD.GF <- apply(x %>% subset(Microbiota=="GF") %>% dplyr::select(mets), 2, function(xx){
      mean(!is.na(xx))
    })
    percent.LOD <- data.frame(Metabolite = mets , SPF = percent.b.LOD.SPF*100 , GF = percent.b.LOD.GF*100 , stringsAsFactors = F)
    percent.LOD$SPF.specific <- percent.LOD$SPF>30 & percent.LOD$GF<30
    percent.LOD$GF.specific <- percent.LOD$SPF<30 & percent.LOD$GF>30
    
    return (percent.LOD)
    
  })
  
  LOD.by.genotype %>% write.csv(. , "./results/LOD.by.genotype.csv")
  
  #Read-in regression results
  geno.res <- read.xlsx("./results/caltech.PD mice.011920.xlsx","conditional.genotype.effect") 
  geno.res <- ddply(geno.res , .(imputation) , function(x){
    x$wilcox.p.FDR <- p.adjust(x$wilcox.p , "BH")
    return(x)
  })
  geno.res_SPF <- geno.res %>%  subset(Microbiota=="SPF" & Tissue=="Plasma")
  geno.res_SPF %>% subset( Metabolite %in% (LOD.by.genotype%>% subset(Material=="Plasma" & SPF.specific))$Metabolite & wilcox.p<0.4)
  
  temp <- geno.res %>% subset(imputation=="unimputed") %>% 
    dplyr::select(Tissue , Metabolite ,Microbiota ,  t , wilcox.p , wilcox.p.FDR) %>%
    as.data.table %>%
    dcast(Tissue+Metabolite~Microbiota , value.var= c("t","wilcox.p","wilcox.p.FDR"))

  temp %>% subset(wilcox.p_SPF <0.05 & wilcox.p_GF>0.05 ) %>% arrange(Tissue, wilcox.p_SPF) %>%
    write.csv(. , "./results/sig WT vs. ASO_SPF not GF.csv")
  

  microb.res <- read.xlsx("./results/caltech.PD mice.011920.xlsx","conditional.microbita.effect") 
  microb.res <- ddply(microb.res , .(imputation,Genotype,Tissue) , function(x){
    x$wilcox.p.FDR <- p.adjust(x$wilcox.p , "BH")
    return(x)
  })
  microb.res %>% dplyr::select(imputation , Genotype , Tissue , Metabolite , Description , t, wilcox.p , wilcox.p.FDR) %>%
    subset(wilcox.p<0.05) %>%
    order(Tissue, wilcox.p) %>%
    write.csv(. , "./results/ASO-GF vs ASO-SPF.csv")
  
  
  
  rr <- ddply(dta.imputed , .(Tissue) , function(X){
    material = X$Material[1]
    res <- c()
    for ( met in mets.included[[material]]){
      corr = cor.test(
        sprintf("~BW + `%s`" , met) %>% as.formula , 
        X ,
        method = "spearman"
      )
      res <- rbind(res , 
                   data.frame(metabolite = met , p = corr$p.value , cor = corr$estimate , stringsAsFactors = F))
    }
    res$p.BH = p.adjust(res$p ,  "BH")
    res$p.Bonf = p.adjust(res$p ,  "bonferroni")
    
    return(res)
  })
  