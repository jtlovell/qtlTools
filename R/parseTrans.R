
#' @title Find the location of a cis eQTL
#'
#' @description
#' \code{parseTrans} Find the mapping position of the highest lod score
#' within a window surrounding the known position of a gene.
#'
#' @param cross The qtl cross object with marker names that need to be changed.
#' @param phe Character or numeric vector indicating the phenotype to be tested
#' @param chr.gene The name of the chromosome that contains the focal genes
#' @param pos.gene The suspected mapping position of the gene - may be inferred from bp2cm
#' @param cmWindow The distance from pos.gene that the cis eQTL may reside. Default is 10, but
#' if the confidence of the mapping position is high, a smaller value is suggested.
#' @param outputType Whether to return a data.frame or qtl object (from makeqtl)
#' @param lodThreshold If outputType = data.frame, the presence of a cis.eqtl is defined by whether
#' any marker within cmWindow has a lod score > lodThreshold.
#' @param lodThreshold

#'
#' @return Either a qtl object or a dataframe, indicating the position of a cis eQTL.
#'
#' @import qtl
#' @export
parseTrans<-function(cross, phe,
                     sdrop,
                     formula = NULL,
                     qtl,
                     cisQTL = NULL,
                     covar = NULL,
                     cisdist = 40,
                     verbose = FALSE,
                     retainCis = TRUE){
  if(!is.null(formula)){
    attr(qtl, "formula") <- deparseQTLformula(formula)
  }
  sdrop<-transs$dropstats
  qtl<-transs$model

  if(retainCis & !is.null(cisQTL)){
    cis.chr=cisQTL$chr
    cis.pos=cisQTL$pos
    cis.wind=c(cis.pos-cisdist, cis.pos+cisdist)
    is.cisqtl<-qtl$altname[qtl$chr==cis.chr &
                          qtl$pos>cis.wind[1] &
                          qtl$pos<cis.wind[2]]
    qtl2check<-qtl$altname[qtl$altname != is.cisqtl]
  }else{
    qtl2check<-qtl$altname
  }

  # drop non-significant QTL and interactions...
  for(i in qtl2check){
    qnames<-rownames(sdrop)[grepl(i, rownames(sdrop), fixed=T)]
    pvals<-sdrop$Pvalue.F.[grepl(i, rownames(sdrop), fixed=T)]
    if(all(pvals>0.05) & length(pvals>0)){
      # for(j in qnames){
      #   form.allint<-gsub(paste(j,"+", sep=" "),"",form.allint, fixed=T)
      # }
      # form.allint<-as.formula(form.allint)
      # print(form.allint)
      qnamed<-qnames[!grepl(":",qnames, fixed=T)]
      formula <- as.formula(attr(qtl, "formula"))

      ts<-terms(as.formula(deparseQTLformula(formula)))
      todrop<- grep(qnamed, attr( ts, "term.labels") )
      nt <- drop.terms(ts, dropx = todrop, keep.response = TRUE)
      formula<-deparseQTLformula(reformulate(attr( nt, "term.labels"), response = "y"))

      qtl<-dropfromqtl(qtl, qtl.name=qnamed)

      chr <- qtl$chr
      pos <- qtl$pos

      for(i in 1:nqtl(qtl)){
        an<-qtl$altname
        n<-qtl$name
        formula<-gsub(n,an,formula)
      }
      qtl$name<-qtl$altname
      attr(qtl, "formula") <- deparseQTLformula(formula)

      fall<-fitqtl(cross, pheno.col=phe, formula=formula(qtl), qtl=qtl, covar=covar, method="hk")
      sall<-summary(fall)
      sdrop<-data.frame(sall$result.drop)
    }
  }
  stats<-qtlStats(cross, mod=qtl, phe=phe, cisQTL = cisQTL, covar = covar, form = formula(qtl), calcConfint=F)

  return(list(stats=stats, model=qtl))
}
