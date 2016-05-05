#' @title Calculate all confidence intervals from a QTL model
#'
#' @import qtl
#' @export
calcCis<-function(mod, qtlnames, ci.method, drop=1.5, prob=0.95, returnChr=FALSE, returnMaxLod=FALSE, ...){
  cis<-data.frame()
  for (j in 1:nqtl(mod)){
    if(ci.method=="drop"){ciout<-lodint(mod,qtl.index=j, drop=drop, ...)
    }else{
      if(ci.method=="bayes"){ciout<-bayesint(mod,qtl.index=j, prob=prob, ...)}
    }
    lowmarker<-rownames(ciout)[1]
    highmarker<-rownames(ciout)[3]
    lowposition<-as.numeric(ciout[1,2])
    highposition<-as.numeric(ciout[3,2])
    if(returnMaxLod){
      maxLod<-as.numeric(sapply(attr(mod, "lodprofile"), function(x) max(x$lod))[[j]])
    }
    if(returnChr){
      chr=mod$chr[j]
      pos=mod$pos[j]
      if(!returnMaxLod){
        cis<-rbind(cis,data.frame(chr, pos,lowmarker,highmarker,lowposition,highposition))
      }else{
        cis<-rbind(cis,data.frame(chr, pos, maxLod, lowmarker,highmarker,lowposition,highposition))
      }
    }else{
      if(!returnMaxLod){
        cis<-rbind(cis,data.frame(lowmarker,highmarker,lowposition,highposition))
      }else{
        cis<-rbind(cis,data.frame(pos, maxLod, lowmarker,highmarker,lowposition,highposition))
      }

    }
  }
  cis$lowposition<-as.numeric(cis$lowposition)
  cis$highposition<-as.numeric(cis$highposition)
  return(data.frame(qtlnames, cis))
}
