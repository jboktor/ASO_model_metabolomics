# cor scrap

anal.dta <- read.csv("files/caltechdata.csv" , stringsAsFactors = F , check.names = F)
mets <- colnames(anal.dta) [10:643]

corr <- plyr::dlply(anal.dta , .(imputation), function(X){
  print(unique(X$imputation))
  dat <-dcast(X %>% as.data.table,Animal~Tissue , value.var = mets, sep = "___") %>% dplyr::select(-1)
  return( psych::corr.test(dat , method="spearman" , ci=F) )
})

library(openxlsx)
wb <- createWorkbook()

addWorksheet(wb,"imputed data - Spearman r")
writeDataTable(wb, "imputed data - Spearman r" ,corr$imputed$r %>% data.frame (., check.names = F) , rowNames = T)

addWorksheet(wb,"imputed data - Spearman r-p")
writeDataTable(wb, "imputed data - Spearman r-p" ,corr$imputed$p%>% data.frame (., check.names = F), rowNames = T)

addWorksheet(wb,"unimputed data - Spearman r")
writeDataTable(wb, "unimputed data - Spearman r" ,corr$unimputed$r%>% data.frame (., check.names = F), rowNames = T)

addWorksheet(wb,"unimputed data - Spearman r-p")
writeDataTable(wb, "unimputed data - Spearman r-p" ,corr$unimputed$p%>% data.frame (., check.names = F), rowNames = T)
saveWorkbook(wb, "figures/Correlations/caltech PD mice model Cor.060921.xlsx", overwrite = TRUE) 

