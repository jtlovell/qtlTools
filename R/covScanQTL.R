#' @title Find the effect of candidate genes based on gene expression.
#'
#' @description
#' \code{covScanQTL} Employs the covariate scan approach (Lovell et al. (2015),
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
#' @examples
#' # See tutorial - findingCandidateGenes
#'
#' @import qtl
#' @export
covScanQTL<-function(cross, pheno.y, pheno.candidates,
                     experimental.covar = NULL,
                     mod = NULL, QTLxE = TRUE,
                     chr, pos, qtl.method = "hk",
                     nperm = 0, plotit=TRUE){

  extractGenoCalls<-function(cross, model, covar = NULL, threshold = 0){
    temp<-model[[1]]
    for(i in 1:length(temp)) temp[[i]][temp[[i]]<=threshold]<-NA
    
    gp <- lapply(temp, function(x){
      apply(x, 1, function(y) {
        ifelse(sum(is.na(y)) == length(y),NA, which.max(y))
      })
    })
    gp2 <- data.frame(do.call(cbind, gp))
    colnames(gp2) <- mod$altname
    if(!is.null(covar)){
      gp2<-cbind(covar, gp2)
    }
    return(gp2)
  }
  
  
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
  if(plotit) plot(bs, col = "darkred", main = pheno.y, xlab = paste("Chr",chr,"Map position (cM)"))

  wh<-which.min(abs(bs$pos - pos))
  bsl<-as.numeric(bs[wh,-c(1:2)])

  maxScans<-sapply(pheno.candidates, function(x){
    covTemp<-data.frame(covBase, exp = pull.pheno(cross,x))
    cs<-scanone(cross, pheno.col = pheno.y, method=qtl.method,
                addcovar = covTemp, intcovar = covBase.int, chr = chr)
    if(plotit) plot(cs, add = T, col = rgb(0,0,0,.2), lty=1)
    return(as.numeric(cs[wh,-c(1:2)]))
  })
  diffScans<- bsl-maxScans
  rankScans<-rank(-diffScans)
  if(nperm>0){
    permScans<-sapply(1:nperm, function(x){
      perm.phes<-pull.pheno(cross,pheno.candidates)
      perm.phes<-perm.phes[sample(1:nrow(perm.phes)),]
      sapply(pheno.candidates, function(y){
        covTemp<-data.frame(covBase, exp = perm.phes[,y])
        cs<-scanone(cross, pheno.col = pheno.y, method=qtl.method,
                    addcovar = covTemp, intcovar = covBase.int, chr = chr)
        return(as.numeric(cs[wh,-c(1:2)]))
      })
    })
    ps<-sapply(names(maxScans), 
               function(x) sum(permScans[x,]<=maxScans[x])/nperm)
  }else{
    ps<-NA
  }

  out<-data.frame(phenotype = pheno.y, 
                  candidateID = pheno.candidates, 
                  lodAtPeak = maxScans, 
                  diffAtPeak = diffScans, 
                  rank = rankScans, 
                  perm.p = ps,
                  stringsAsFactors=F)
  rownames(out)<-NULL
  out<-out[order(out$rank),]
  return(out)
}


