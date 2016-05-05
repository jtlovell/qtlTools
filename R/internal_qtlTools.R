#' @title A set of functions needed to run the main qtlTools functions
#'
#' @import qtl
#' @import effects
#' @export
calcQtlMeans<-function(cross, mod, covar, dropstats, phe, form, ... ){
  gp<-lapply(mod[[1]], function(x) apply(x,1, function(y) which(y==max(y))))
  geno.names<-colnames(mod[[1]][[1]])
  gp2<-data.frame(do.call(cbind,gp))
  colnames(gp2)<-mod$altname
  for(i in mod$altname) gp2[,i]<-as.factor(gp2[,i])
  phe.num<-pull.pheno(cross, pheno.col=phe)
  phe.num<-phe.num[!is.na(phe.num)]
  rownames(gp2)<-getid(cross)

  if(any(grep(":",rownames(dropstats)))){
    intnames<-rownames(dropstats)[grep(":",rownames(dropstats))]
    hasint<-as.character(sapply(intnames, function(x) strsplit(x, ":")[[1]][1]))
    gpi<-data.frame(gp2[,colnames(gp2) %in% hasint])
    colnames(gpi)<-intnames
    for(i in colnames(gpi)) gpi[,i]<-paste(gpi[,i], covar, sep="_")
    gp2<-cbind(gp2, gpi)
  }
  gp4<-cbind(gp2, y = as.vector(phe.num))

  lm.out<-lm(form, data=data.frame(gp4))
  test<-allEffects(lm.out)
  effects<-lapply(test,function(x) as.numeric(x$fit))
  ses<-lapply(test,function(x) as.numeric(x$se))
  out<-data.frame(qtlnames = mod$name,  do.call(rbind, effects), do.call(rbind, ses))
  colnames(out)[-1]<-c(paste(geno.names,c("effect"),sep="_"), paste(geno.names,c("se"), sep="_"))
  return(out)
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

makeCumPos<-function(chr, pos, gap = NULL, return.centers=FALSE){
  totlength<-sum(tapply(pos, chr, max))
  if(is.null(gap)) gap = totlength*.02
  maxs<-tapply(pos, chr, max)
  starts<-cumsum(c(0,maxs[-length(maxs)]+gap))
  index<-data.frame(chr=names(maxs), starts = starts, stringsAsFactors = F)
  for(i in index$chr){
    pos[chr==i]<-pos[chr==i] + index$starts[index$chr == i]
  }
  if(return.centers){
    ends<-cumsum(c(maxs[1], maxs[-1]+gap))
    out<-(starts+ends)/2
    names(out)<-names(maxs)
    return(out)
  }else{
    return(pos)
  }
}

bootDens<-function(x, binwidth=5, nboot=1000, thresh=.95){
  out<-sapply(1:nboot, function(i){
    n=length(x)
    values=unique(x)
    booted<-sample(values, n, replace=T)
    max(hist(booted, breaks=ceiling(max(values)/binwidth), plot=F)$counts)
  })
  if(!is.null(thresh)){
    quantile(out, thresh)
  }else{
    out
  }
}



