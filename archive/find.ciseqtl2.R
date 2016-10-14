find.ciseqtl2<-function(cross, phe, pos.gene, chr.gene, window = 10, qtlmethod = "hk",output = "model", verbose = T, ...){
  if(verbose) cat("running scanone on", length(phe),"traits\n", sep = " ")
  s1<<-scanone(cross, pheno.col = phe, method = qtlmethod, ...)
  if(verbose) cat("finding local qtl peak within", window, "cM of cis position\n", sep =" ")
  s1.info<-data.frame(s1[,c("chr","pos")], stringsAsFactors=F)
  s1.info$chr <- as.numeric(as.character(s1.info$chr))
  s1.lods<-data.frame(s1[,-c(1:2)])
  if(ncol(s1.lods)==1) colnames(s1.lods)<-phe
  s1.wind = lapply(phe, function(x){
    pos = pos.gene[phe==x]; chr = chr.gene[phe==x]
    cat(pos, chr)
    index<-which(s1.info$chr ==  chr & abs(s1.info$pos-pos)<window)
    s1.out<-cbind(s1.info[index,],lod = s1.lods[index,x])
    return(s1.out[which.max(s1.out$lod)[1],1:2])
  })
  out<-do.call(rbind, s1.wind)
  rownames(out)<-phe
  if(output == "data.frame"){
    return(out)
  }else{
    if(verbose) cat("constructing qtl models\n")
    mtype<-ifelse(qtlmethod == "hk", "prob","draws")
    lout<-lapply(phe, function(x){
      makeqtl(cross, chr = out[x,"chr"], pos =  out[x,"pos"], what = mtype)
    })
    names(lout)<-phe
    return(lout)
  }
}


calcResidPheno<-function(cross, models, formula, covar, phe, verbose = T, ...){
  if(verbose) cat("extracting genotype calls from",length(phe),"models\n", sep = " ")
  phemat<-pull.pheno(cross, pheno.col = phe)
  gcs<-lapply(phe, function(x) {
    gc<-extractGenoCalls(cross, models[[x]], covar)
    gc$y = phemat[,x]
    return(gc)
  })
  if(verbose) cat("calculating residuals from",formula,"\n", sep = " ")
  resids<-lapply(gcs, function(x) lm(as.formula(formula), data = x, ...)$resid)
  if(verbose) cat("compiling results\n", sep = " ")
  resids<-do.call(cbind,resids)
  colnames(resids)<-phe
  return(resids)
}

extractGenoCalls<-function(cross, model, covar = NULL, threshold = 0){
  temp<-model[[1]]
  for(i in 1:length(temp)) temp[[i]][temp[[i]]<=threshold]<-NA

  gp <- lapply(temp, function(x){
    apply(x, 1, function(y) {
      ifelse(sum(is.na(y)) == length(y),NA, which.max(y))
    })
  })
  gp2 <- data.frame(do.call(cbind, gp))
  colnames(gp2) <- model$altname
  if(!is.null(covar)){
    gp2<-cbind(covar, gp2)
  }
  return(gp2)
}
