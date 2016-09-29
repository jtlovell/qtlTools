#' @import qtl
#' @export
GWERk<-function(cross.in, phe, k=1, nperms=1000, permstrata = NULL, batchit=TRUE,...){
  phenos<-pull.pheno(cross.in, pheno.col=phe)
  n.ind<-nind(cross.in)
  if(is.null(permstrata)){
    permstrata<-rep(1, n.ind)
  }

  u <- unique(perm.strata)
  if(batchit){
    ord <- matrix(0, n.ind, nperms)
    for(i in u){
      for(j in 1:nperms){
        ord[permstrata==i,j]<-sample(phenos[permstrata==i])
      }
    }
    cross.in$pheno <- cbind(matrix(cross.in$pheno[,phe][ord], nrow=n.ind), cross.in$pheno)
    tem<-scanone(cross=cross.in, pheno.col = 1:nperms)
  }else{
    temp<-list()
    for(i in 1:nperms){
      ord <- 1:n.ind
      for(j in u) {
        wh <- perm.strata==j
        if(sum(wh)>1) ord[wh] <- sample(ord[wh])
      }
      ord<-matrix(ord, ncol=1)
      cross.in$pheno <- cross.in$pheno[o,phe,drop=FALSE]
      temp[[i]] <- scanone(cross=cross.in, pheno.col = phe, ...)
    }
    tem<-do.call(cbind, temp)
  }
  res <- matrix(apply(tem[,-(1:2),drop=FALSE],2, function(x) {
    m<-tapply(x,tem$chr,max, na.rm=TRUE)
    m[order(-m)][k+1]
  }), ncol=1)
  colnames(res) <- "lod"
  return(res)
}
