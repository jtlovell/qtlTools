#' @title Method to summarize scanone results
#'
#' @description
#' \code{pullSigQTL} Uses qtlpvl to summarize the output of scanone,
#' then culls the output to only significant QTL peaks, based on permutations.
#'
#' @param cross The qtl cross.
#' @param pheno.col Character or numeric vector indicating the phenotype to be tested.
#' @param s1.output The output from scanone
#' @param perm.output The permutation output from scanone
#' @param chr The chromosome to be tested. Defaults to all chromosomes.
#' @param alpha The significance for permutations
#' @param returnQTLModel Logical, should a QTL model be returned (TRUE), or
#' should a culled output from qtlpvl::convert_scan1 be returned (FALSE)?
#' @param ... additional arguments passed on to summary.scanone, such as controlAcrossCol.
#' @return Either QTL models or simplified and converted scanone summary.
#'
#' @examples
#' library(qtlpvl)
#' @export
pullSigQTL<-function(cross, s1.output, perm.output, pheno.col, chr=NULL, alpha = 0.05, returnQTLModel = TRUE,
                     ...){
  if(is.null(chr)) chr <- chrnames(cross)
  maxs<-convert_scan1(s1, phenoname=phes, chr = chr)
  sperms<-summary(perms, alpha = alpha, ...)
  perm.names<-attr(sperms,"dimnames")[[2]]
  sperms<-as.numeric(sperms)
  for(i in pheno.col){
    thresh<-sperms[perm.names==i]
    maxs<-maxs[-which(maxs$pheno == i & maxs$lod1 < thresh),]
  }
  if(returnQTLModel){
    mods<-lapply(pheno.col, function(i){
      temp<-maxs[maxs$pheno == i,]
      if(nrow(temp)==0) {
        out<-"NULL QTL Model"
      }else{
        if("prob" %in% names(cross$geno[[1]])){
          out<-makeqtl(cross, chr = temp$chr, pos = temp$pos, what = "prob")
        }else{
          out<-makeqtl(cross, chr = temp$chr, pos = temp$pos, what = "draws")
        }
      }
    })
    names(mods)<-pheno.col
    return(mods)
  }else{
    return(maxs)
  }
}
