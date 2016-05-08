#' @title Make a cis eQTL using permutation p-values
#'
#' @import qtl
#' @export
makeCisQTL<-function(annot, cross, phe, covar, outputType="qtl", perms){

  perm.add<-as.numeric(perms[["add"]])
  perm.int<-as.numeric(perms[["int"]])
  perm.diff<-as.numeric(perms[["diff"]])

  cis.add<- find.ciseqtl(cross=cross, phe=phe, pos.gene=annot$predictedCM, chr.gene=annot$lg,
                         method="hk", covar=covar, addcovar=covar, intcovar=NULL, outputType="qtl")
  cis.int <- find.ciseqtl(cross=cross, phe=phe, pos.gene=annot$predictedCM, chr.gene=annot$lg,
                          method="hk", covar=covar, addcovar=covar, intcovar=covar, outputType="qtl")
  mods<-list(add = cis.add[[1]], int = cis.int[[1]])
  lods<-c(add = cis.add[[2]]$lod, int = cis.int[[2]]$lod)
  pvals<-c(sum(perm.add>lods[1])/length(perm.add),
           sum(perm.int>lods[2])/length(perm.int),
           sum(perm.diff>(lods[2]-lods[1]))/length(perm.diff))
  if(pvals[2]<pvals[1] | pvals[2] == pvals[1] & pvals[3]<=0.05){
    model.out<-mods[[2]]
    model.out$formula = "y ~ Q1 + Q1:covar + covar"
  }else{
    model.out<-mods[[1]]
    model.out$formula = "y ~ Q1 + covar"
  }
  return(model.out)
}
