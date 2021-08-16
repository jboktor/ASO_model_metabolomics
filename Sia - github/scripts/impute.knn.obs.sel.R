# Function to impute data points in metabolomics dataset using KNN 
# using both KNN in sample as well as in metabolite dimensions
impute.knn.obs.sel <- function(dat, K=10) {
  
  results <- list()
  cor.cutoff <- 0.2     # use only variables with cor>0.2 for distance computation
  
  da1 <- dat 
  da1list <- da2list <- rep(list(dat),length(K)) 
  
  incom.vars <- which(apply(dat,2,function(x) any(is.na(x))))
  incom.obs <- which(apply(da1,1,function(x) any(is.na(x))))
  
  Cor <- cor(da1,use="p")
  
  D2list <- lapply(incom.vars, function(j) {
    varsel <- which(abs(Cor[j,])>cor.cutoff)  
    if(length(varsel)>10) varsel <- order(abs(Cor[j,]),decreasing=T)[1:11]
    if(length(varsel)<5) varsel <- order(abs(Cor[j,]),decreasing=T)[1:6]
    D2 <- as.matrix(dist(scale(da1[,varsel])),upper=T,diag=T) 
    if(any(is.na(D2))) {
      D2a <- as.matrix(dist(scale(da1)),upper=T,diag=T)*sqrt(length(varsel)/ncol(da1)) 
      D2[is.na(D2)] <- D2a[is.na(D2)] 
    }
    diag(D2) <- NA
    D2})
  names(D2list) <- incom.vars
  
  for (i in incom.obs){
    comvars <-  complete.cases(as.numeric(da1[i,]))
    for (j in which(!comvars)) {
      D2 <- D2list[[as.character(j)]]                                 
      if(any(!is.na(D2[i,]))) {
        KNNids <- order(D2[i,],na.last=NA)
        KNNids_naomit <- KNNids[sapply(KNNids,function(ii) any(!is.na(da1[ii,j])))] 
      } else {
        KNNids  <- NULL
      }
      da1list <- lapply(1:length(da1list),function(ii) {
        k <- K[ii] 
        da <-  da1list[[ii]]
        if(!is.null(KNNids)) {
          KNNids_sel <- intersect(KNNids[1:min(k,length(KNNids))],KNNids_naomit)
        }
        if(length(KNNids_sel)<1) {
          KNNids_sel <- KNNids_naomit[1:min(floor(k/2),length(KNNids_naomit))]
        } else if (length(which(sapply(KNNids_sel,function(ii) !is.na(da1[ii,j])))) < floor(k/2) )  KNNids_sel <- KNNids_naomit[1:min(floor(k/2),length(KNNids_naomit))] 
        if(any(!is.na(D2[i,])) & length(KNNids)>=1) {
          da_sel <- da[KNNids_sel,j]
          da[i,j] <- sum(da_sel*exp(-D2[i,KNNids_sel]),na.rm=T)/sum(exp(-D2[i,KNNids_sel])[!is.na(da_sel)],na.rm=T) }
        da}) 
    }
  }
  da1list <- lapply(da1list, function(da) {
    da <- apply(da,2, function(x) {
      if(any(is.na(x))) x[is.na(x)] <- mean(x,na.rm=T)
      x}) 
    da})
  
  results <- c(results,list(da1list))
  names(results)[length(results)] <- "knn.sample.euc.sel" 
  rm(da1list,da2list)  
  return(results$knn.sample.euc.sel[[1]])
}