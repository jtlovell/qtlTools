#' @title Assign eQTL as cis/trans
#'
#' @description
#' \code{assignCisTrans} An extension of qtlpvl::get.trans.info, but permitting
#' assignment using confidence intervals.
#' @param s1.output Scanone output to be passed to calcCis
#' @param perms.output Permutation output to be passed to calcCis
#' @param gff.cm The gff with mapping position, as calculated from findGenecM
#' @param useConfInt Logical, should the confidence intervals be used to assign
#' cis/trans to eQTL?
#' @param cmInterval If useConfInt = FALSE, the cM distance the eQTL peak
#' needs to be from the physical position to constitute a cis-eQTL?
#' @param additional parameters passed to calcCis
#'
#'
#' @import qtl
#' @export
assignCisTrans<-function(s1.output, perms.output, gff.cm, useConfInt = TRUE, cmInterval = 10, ...){
  cis<-calcCis(s1.output=s1.output,
               perm.output=perms.output,...)
  gff<-gff.cm[gff.cm$geneID %in% cis$pheno,c("geneID","chr","pos")]
  names(gff)<-c("pheno","chr0","pos0")
  out<-merge(gff,cis, by="pheno")
  if(!useConfInt){
    out$is.cis <- (out$chr0 == out$chr) & abs(out$pos0 - out$pos) < cmInterval
  }else{
    out$is.cis<- (out$chr0 == out$chr) & (out$pos0 >= out$lowposition) &
      (out$pos0 <= out$highposition)
  }
  out$is.trans <- !out$is.cis
  return(out)
}
