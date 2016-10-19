#' @title Find mutual ranks between gene expression data
#'
#' @import qtl
#' @export
mutual.rank<-function(x, verbose = TRUE, format = "clusterone", corThresh = NULL, expdecay = 25){

  if(verbose) cat("... calculating gene-wise ranks\n")
  y<-data.table(x)
  y<-y[,lapply(.SD,function(i) rank(-i)),.SDcols=1:ncol(y)]
  rmat<-data.matrix(y)

  # set ranks with r value < designated threshold to NA
  if(!is.null(corThresh)){
    if(verbose) cat("... setting values below threshold to NA\n")
    cmat<-corFast(rmat)
    rmat[x < corThresh]<-NA
  }

  # create transposed matrices so that mutual ranks can be easily calculated
  if(verbose) cat("... transposing matrices\n")
  lmat<-rmat-1
  umat<-rmat-1
  umat[lower.tri(umat, diag = TRUE)]<-NA
  lmat[upper.tri(lmat, diag = TRUE)]<-NA
  lmat<-t(lmat)

  #calculate mutual ranks
  if(verbose) cat("... calculating mutual ranks\n")
  out<-sqrt(umat * lmat)

  if(format == "clusterone"){
    if(verbose) cat("... generating clusterone - formatted output\n")
    rownames(out)<-colnames(out)
    out<-melt(out)
    out<-out[complete.cases(out),]
    colnames(out)<-c("gene1","gene2","mutual.rank")
    if(!is.null(expdecay)){
      if(verbose) cat("... calculating clusterone edge weights\n")
      # exponential decay function
      out$edge.weight<- exp(-(out$mutual.rank-1)/expdecay)
    }
  }
  return(out)
}

