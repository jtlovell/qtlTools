#' @title Add a qtl to a model by permutation p-values
#'
#' @import qtl
#' @export
add1QTL<-function(cross, phe, covar, perms, startModel){
  perm.add<-as.numeric(perms[["add"]])
  perm.int<-as.numeric(perms[["int"]])
  perm.diff<-as.numeric(perms[["diff"]])

  covname<-colnames(covar)[1]
  nq<-nqtl(startModel)
  form<-startModel$formula
  newforms<-paste(form, paste("Q",nq+1, sep=""), sep = " + ")
  newforms<-c(newforms, paste(newforms, paste(paste("Q",nq+1, sep=""),":",covname, sep=""), sep = " + "))

  addscans<-lapply(newforms, function(x) {
    aq<-addqtl(cross, qtl=startModel, formula=x, model="normal",
               method="hk", covar=covar, pheno.col=phe)
    for(i in 1:nq){
      aq$lod[aq$chr==startModel$chr[i] & abs(aq$pos-startModel$pos[i])<40]<-0
    }
    aq
  })
  maxs<-lapply(addscans, function(x) {
    x<-data.frame(x, stringsAsFactors=F)
    x[x$lod == max(x$lod),]
  })
  max.add<-maxs[[1]][c("chr","pos")]
  max.int<-maxs[[2]][c("chr","pos")]
  add.scan<-data.frame(addscans[[1]], stringsAsFactors=F)
  add.scan$chr = as.character(add.scan$chr)
  add.scan$pos = as.numeric(add.scan$pos)
  max.add.int<-add.scan[add.scan$pos == as.numeric(max.int[2]) & add.scan$chr == as.character(max.int[1]),]
  lods<-c(sapply(maxs, function(x) x$lod),max.add.int$lod)
  pvals<-c(sum(perm.add>lods[1])/length(perm.add),
           sum(perm.int>lods[2])/length(perm.int),
           sum(perm.diff>(lods[2]-lods[3]))/length(perm.diff))
  if(pvals[2] < pvals[1] | pvals[2] <= 0.05 & pvals[3] <= 0.05){
    model.out<-addtoqtl(cross, qtl = startModel, chr = max.int$chr, pos = max.int$pos)
    model.out$formula = newforms[2]
    pout<-pvals[2]
  }else{
    model.out<-addtoqtl(cross, qtl = startModel, chr = max.add$chr, pos = max.add$pos)
    model.out$formula = newforms[1]
    pout<-pvals[1]
  }
  return(list(model = model.out, newQTL.pvalue = pout))
}
