#' @import qtl
#' @export
GWERk<-function(cross.in, phe, k=1, nperms=1000, printUpdate=TRUE, sameSeed = TRUE, ...){

  out<-lapply(phe, function(i){
    if(sameSeed) set.seed(42)
    if(printUpdate) cat(i, "\n")

    phe.in<-pull.pheno(cross.in, pheno.col=i)

    make.perm<-lapply(1:nperms, function(x)  sample(phe.in, size=length(phe.in)))

    permd.phe<-data.frame(do.call(cbind,make.perm))
    colnames(permd.phe)<-paste("perm",1:nperms,sep="_")

    cross2<-cross.in
    cross2$pheno <- cbind(cross2$pheno, permd.phe)

    perm<-scanone(cross2, pheno.col=colnames(permd.phe), ...)

    sperm<-summary(perm)[,colnames(permd.phe)]
    sum.perm<-as.numeric(apply(sperm, 2, function(x) {
      for (i in 1:k) x<-x[-which(x==max(x))]
      return(max(x))
    }))
    return(sum.perm)
  })
  names(out)<-phe
  return(do.call(cbind, out))
}
