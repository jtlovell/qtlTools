#' @title Find the effect of candidate genes based on gene expression.
#'
#' @description
#' \code{covScanQTL} Employ the covariate scan approach (Lovell et al. (2015),
#' Plant Cell), to rank the potential of candidate genes based on their effect
#' on QTL morphology. Can be run for a single phenotype (e.g. Lovell et al. (2015)),
#' or on a set of QTL underlying a trand-band. In the latter case, a boxplot of ranks
#' can be output.
#'
#' @param cross The qtl cross object with marker names that need to be changed.
#' @param pheno.y Character vector specifying the phenotypes modelled. Is passed
#' to the "pheno.col" argument of qtl::scanone.
#' @param pheno.candidates Character vector specifying the names of the expression phenotypes to
#' be extracted from the cross object that serve as the covariate.
#' @param exp.covar data.frame containing the experimental covariate. Passed on to
#' addcovar in qtl::scanone. If QTLxE is TRUE, also passed to intcovar
#' @param mod A qtl model object that can supply additional marker covariates that
#' are held constant during the covariate scan.
#' @param QTLxE Logical, should the focal QTL have an interactive covariate.
#' See exp.covar.
#' @param chr Numeric or character, of length 1, specifying the chromosome on which
#' the QTL peak to scan exists
#' @param pos Numeric, of length 1, specifying the QTL peak position on chromosome chr.
#' @param qtl.method The method passed to scanone
#' @param nperm If permutation tests are desired, specify as >0. These are not currently
#' computationally efficient and take forever when using eQTL data.
#' @param plotit Logical, when more than 1 pheno.y is specified, presents a boxplot of
#' covariate scan ranks.
#'
#' @return Either a list or a dataframe, containing the maximum scanone outputs at
#' chromosome chr and position pos for each phenotype and covariate. If nperm > 0, also
#' gives the permutation results as a second list element.
#'
#' @import qtl
#' @importsfrom reshape2 melt
#' @export
covScanQTL<-function(cross, pheno.y, pheno.candidates, experimental.covar = NULL,
                     mod = NULL, QTLxE = TRUE,
                     chr, pos, qtl.method = "hk", nperm = 0, plotit=TRUE){

  if(!is.null(mod)){
    covBase<-extractGenoCalls(cross, model=mod, covar = experimental.covar, threshold = 0)
  }else{
    covBase<-experimental.covar
  }

  if(QTLxE){
    covBase.int<-experimental.covar
  }else{
    covBase.int<-NULL
  }

  bs<-scanone(cross, pheno.col = pheno.y, method=qtl.method,
              addcovar = covBase, intcovar = covBase.int, chr = chr)

  wh<-which.min(abs(bs$pos - pos))
  bsl<-data.frame(pheno = names(bs[wh,-c(1:2)]), baseLOD = as.numeric(bs[wh,-c(1:2)]))

  maxScans<-data.frame(t(sapply(pheno.candidates, function(x){
    covTemp<-data.frame(covBase, exp = pull.pheno(cross,x))
    cs<-scanone(cross, pheno.col = pheno.y, method=qtl.method,
                addcovar = covTemp, intcovar = covBase.int, chr = chr)
    csl<-as.numeric(cs[wh,-c(1:2)])
    return(csl)
  })))
  colnames(maxScans)<-pheno.y
  bsl.mat<-matrix(rep(bsl$baseLOD,length(pheno.candidates)),byrow=T,nrow=length(pheno.candidates))
  diffScans<- bsl.mat-maxScans
  rankScans<-apply(-diffScans, 2, rank)
  return(list(ranks = rankScans, differences=diffScans, initial.s1max=bsl))
}


