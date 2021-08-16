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

############ANOVA Analysis; test G, M , G : M Iteraction 
#1 SPF vs GF
#2- SPF ASO x SPF WT 2)SPF ASO x GF ASO 3) GF WT x SPF WT 4) GF WT x GF ASO
#1- Blood and Brain

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
  
  ###Save Results
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
  
  ##Change reference group for genotype (WT as the reference)
  saveWorkbook(wb, "./results/caltech.PD mice.040820.xlsx", overwrite = TRUE) 
