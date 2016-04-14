#' @title A set of functions needed to run the main qtlTools functions
#'
#' @import qtl
#' @export
calcQtlMeans<-function(cross, mod, covar, ... ){
  if(!is.null(cross$geno[[1]]$draws)){
    cross<-sim.geno(cross, ...)
  }
  qtlnames<-mod$alt.name
  meanse<-lapply(qtlnames[mod$chr!="X"], function(x) {
    if(!is.null(covar)){
      e<-effectplot(cross, mname1="covar", mark1=covar, geno1=unique(covar),
                    mname2=x,draw=F)$Means
      enames<-unlist(sapply(colnames(e), function(y) {
        y<-strsplit(y,"[.]")[[1]][2]
        paste(rownames(e),y, sep="_")}))
      eout<-data.frame(qtlnames = x, t(as.numeric(e)))
      colnames(eout)[-1]<-enames
    }else{
      eout<-data.frame(qtlnames = x, t(effectplot(cross, pheno.col=phe, mname1=x,draw=F)$Means))
      enames<-as.character(unlist(sapply(colnames(eout)[-1], function(y) {
        strsplit(y,"[.]")[[1]][2]
      })))
      colnames(eout)[-1]<-enames
    }
    return(eout)
  })
  meanse<-do.call(rbind, meanse)
  return(meanse)
}

calcCis<-function(mod, qtlnames, ci.method, drop, prob){
  cis<-data.frame()
  for (j in 1:nqtl(mod)){
    if(ci.method=="drop"){ciout<-lodint(mod,qtl.index=j, expandtomarkers=F, drop=drop)
    }else{
      if(ci.method=="bayes"){ciout<-bayesint(mod,qtl.index=j, expandtomarkers=F, prob=prob)}
    }
    lowmarker<-rownames(ciout)[1]
    highmarker<-rownames(ciout)[3]
    lowposition<-ciout[1,2]
    highposition<-ciout[3,2]
    cis<-rbind(cis,cbind(lowmarker,highmarker,lowposition,highposition))
  }
  return(data.frame(qtlnames, cis))
}


defineCisTrans<-function(out=out, cistrans, cisdist=cisdist){
  out$type[out$type!="epi"]<-"trans"

  if(length(cistrans)==1){
    out$type[cistrans]<-"cis"
  }else{
    is.cis<-which(out$chr == cistrans[1] & (out$pos - cistrans[2]) < cisdist)
    out$type[is.cis][1]<-"cis"
  }
  if("cis" %in% out$type){
    cis.qtl<-out$form.name[out$type=="cis"]
  }else{
    cis.qtl<-NULL
  }

  if(!is.null(cis.qtl) | length(cis.qtl)>0){
    out$type[grepl("[:]", out$form.name) &
               grepl(cis.qtl, out$form.name) &
               grepl("covar", out$form.name)]<-"cis.covar.int"
    out$type[grepl("[:]", out$form.name) &
               !grepl(cis.qtl, out$form.name) &
               grepl("covar", out$form.name)]<-"trans.covar.int"
    out$type[grepl("[:]", out$form.name) &
               grepl(cis.qtl, out$form.name) &
               !grepl("covar", out$form.name)]<-"cis.trans.epi"
  }
  out$type[out$form.name=="covar"]<-"covar"
  return(out)
}

