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
setwd("/users/siamakdehkordi/Desktop/Duke/Cal Tech PD Mice/")
rm(list=ls())
set.seed(20)   #set random seed for reproducibility purposed

source("./scripts/impute.knn.obs.sel.R")

#Function to identify multivariate outliers based on Mahalanobis distance
getOutliers.Mahalanobis <- function(dat , chi2.q=NULL , robust=F , dataset.name , out.fname , saveplot = F){
  #dat : n by p data frame
  #chi2.q: quantile of chi-squared distribution used to compute cutoff value
  #robust : whether use robust covariance estimation methods or classical
  dat <- as.matrix(scale(dat))
  nsamples <- nrow(dat)
  nmets <- ncol(dat)
  
  if ( is.null(chi2.q)){
    chi2.q = 1 - 0.01/nsamples
  }
  
  if ( robust){
    covv <- cov.mcd(dat)$cov
  }else{
    covv <- cov(dat)
  }
  hh = diag(dat %*% solve(covv) %*% t(dat))
  dist2 = mahalanobis( dat ,center=F, cov =covv)
  cutoff =  qchisq(chi2.q, nmets)  ## dist^2 ~ chi-square */
  
  #prepare filepath for the Mahalanobis Q-Q plot  
  df <- data.frame(sample = rownames(dat) , mahalanobis.distance=dist2 , outlier = ifelse(dist2>cutoff ,T , F))
  df$label = apply(df,1 , function(x) {
    ifelse(gsub(" ","", x["outlier"]) , x["sample"] , "")
  } )
  df <- df[order(df$mahalanobis.distance),]
  df$qt = sort(qchisq(ppoints(nsamples), df = nmets))
  
  df$sample <- as.character(df$sample)
  p<- ggplot(df , aes(x=qt , y = mahalanobis.distance , color=outlier)) +
    geom_point() +  scale_color_manual(values=c("black", "red")) +
    geom_text_repel(aes(label=label)) +
    ggtitle(paste('Q-Q plot \n',  dataset.name) )+
    xlab (substitute(paste(chi[v]^2, "-quantile "), list(v=nmets))) +
    ylab ('Squared Mahalanobis Distances') + theme(legend.position="none")
  if ( saveplot){
    ggsave(out.fname , p)
  }
  #return(which(dist2 > cutoff))
  return(df[,1:3])
}

impute.LOD <- function (dat , lod , mets, imputed.idx.out.fname,mets.percent.missing,dataset.name, saveplot=F , plot.out.fname){
  #The function imputes missing values using LOD/2
  dat.imputed <- sapply(1:length(mets), function(j){
    d <- dat[,c("Plate Bar Code",mets[j])]
    d$rownames <- rownames(d)
    d <- merge(d , lod[,c("Plate Bar Code",mets[j])] , by = "Plate Bar Code" , suffixes = c("",".lod"), sort=F )
    d <- d[match(rownames(dat),d$rownames),]   #add this line to make sure the original d has the exact same row order as dat
    
    idx <- which(is.na(d[,mets[j]]))
    
    if ( length(idx)>0 ){
      d[idx,mets[j]] <-d[idx,paste0(mets[j],".lod")] /2
    }
    return (d[,mets[j]])
  })
  colnames(dat.imputed) <- mets
  dat.imputed <- data.frame(dat.imputed , check.names = F)
  rownames(dat.imputed) <- rownames(dat)
  
  #Save missing index of metabolites   
  dat.status <- melt(t(is.na(dat[,mets])))
  colnames(dat.status) <- c("col" , "row" , "Imputed" )
  dat.status$ORIG.ID <- dat[dat.status$row , ]$ORIG.ID
  write.csv(dat.status , imputed.idx.out.fname , row.names = F)
  if ( saveplot ){
    ###################  
    #produce boxplots for the metabolites whith high %<LOD
    ###################
    dat.long <- melt(t(dat.imputed))
    colnames(dat.long) <- c("col" , "row" , "value" )
    
    dat.long <- merge(dat.long , dat.status )
    dat.long$Material <-  dat[dat.long$row,]$Material
    dat.long <- merge(dat.long , data.frame(met= names(mets.percent.missing) , percent.missing = mets.percent.missing) ,
                      by.x = "col" , by.y = "met")
    dat.long$col2 <- sprintf("%s(%0.2f%%)", dat.long$col ,dat.long$percent.missing *100)
    
    p <- ggplot(subset(dat.long, col %in% names(which(mets.percent.missing>.2 & mets.percent.missing<.4))) , 
                aes(x=col2, y=log2(value))) + 
      geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
      geom_jitter(aes(color=Imputed), position=position_jitter(width=.1, height=0  ) ) +
      xlab("Metabolite")+ylab("log2(value)")  + 
      ggtitle( paste0("Imputed metabolites with >20% & < 40% <LOD\nUsing LOD/2\n" ,dataset.name ))+
      theme(axis.text.x = element_text(angle = 90,vjust=0.5 , hjust = 1, face = "bold"),
            plot.title = element_text(hjust = 0.5))
    ggsave( plot.out.fname, p)
    
  }
  return ( dat.imputed )
}

#####################################################################################
###### Read-in metabolomics data, LOD/ save them in one excell spreadsheet
#####################################################################################
  Q500Mets <- read.xlsx("./data from biocrates/5370_Gut brain axis study_Sub1_Pilot study_Data with plate numbers_2020-09-09.xlsx","Metabolite information")
  LOD <- read.xlsx("./data from biocrates/5370_Gut brain axis study_Sub1_Pilot study_LOD thresholds_2020-09-09.xlsx")
  colnames(LOD) = LOD[4,]
  LOD = LOD[-(1:4),-3]

  PlasmaTissue <-read.xlsx("./data from biocrates/5370_Gut brain axis study_Sub1_Pilot study_Data with plate numbers_2020-09-09.xlsx","Plasma & tissue (µM)") 
  colnames(PlasmaTissue) <- PlasmaTissue[4,]
  PlasmaTissue  = PlasmaTissue[-(1:4),]
  
  Gut  <-   read.xlsx("./data from biocrates/5370_Gut brain axis study_Sub1_Pilot study_Data with plate numbers_2020-09-09.xlsx","Gut content (µM)") 
  colnames(Gut) <- Gut[4,]
  Gut  <- Gut[-(1:4),]
  
  dta.raw <- rbind.fill(PlasmaTissue , Gut) %>% subset(!is.na(`Sample Identification`))
  dta.raw[9:ncol(dta.raw)] <- apply(dta.raw[9:ncol(dta.raw)] , 2 , as.numeric)
  mets = colnames(dta.raw)[9:ncol(dta.raw)]
  #table( dta.raw$Material)
  #brain tissue             Gut content intestinal tract tissue                  plasma 
  #88                      69                      46                            23   
  QC <- dlply(dta.raw , .(Material) , function(x){
    report <- list()

    report$`%LOD` <- apply(x[mets], 2 , function(y){mean(is.na(y))} )
    report$included.mets <- names(report$`%LOD`)[report$`%LOD`<.4]
    
    return (report)
  })
  dta.raw$Animal = str_extract(dta.raw$`Sample Identification` , "[0-9]+$")
  
  dta.raw <- dta.raw %>% dplyr::select(c( setdiff(colnames(dta.raw),mets),mets)) 


  wb <- createWorkbook()
  addWorksheet(wb, "q500 data")
  writeDataTable(wb, "q500 data",dta.raw)
  
  addWorksheet(wb, "%<LOD")
  writeDataTable(wb, "%<LOD",
                 ldply(QC , function(x){x$`%LOD`}))
  
  addWorksheet(wb, "LOD Thresholds")
  writeDataTable(wb, "LOD Thresholds",LOD)
  
  addWorksheet(wb, "Q500 Metabolites Info")
  writeDataTable(wb, "Q500 Metabolites Info",Q500Mets)
  
  saveWorkbook(wb, "./upload to box/5370_Cal Tech PD Mouse Model_Gut brain axis study_09112020.xlsx", overwrite = TRUE)
  
  
