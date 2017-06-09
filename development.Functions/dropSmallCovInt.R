#' @title Reduce covariate interactions following stepwiseqtl.int
#'
#' @description
#' under development. Use with care.
#'
#' \code{dropSmallCovInt}  Drops non-significant QTL*covariate interactions
#' @import qtl
#' @export
dropSmallCovInt<-function(cross, qtl, formula, covar, pheno.col, method = "hk",
                          thresh.diff){
  f<-data.frame(summary(fitqtl(cross, qtl = qtl, formula = formula,
                               covar = covar, method = method,
                               pheno.col = pheno.col))$result.drop)

  f$id = rownames(f)
  for(i in 1:length(qtl$name)) f$id = gsub(qtl$name[i],qtl$altname[i],f$id, fixed = T)
  tocheck<-c(":","covar")
  checked<-sapply(rownames(f), function(x) all(grepl(tocheck[1],x,fixed=T),grepl(tocheck[2],x,fixed=T)))

  fout<-f[!checked | f[,"LOD"]>thresh.diff,]
  newform = paste("y ~ ", paste(fout$id, collapse = " + "))
  attr(qtl,"formula")<-newform
  newqtl = refineqtl(cross, pheno.col = pheno.col, qtl = qtl,
                     formula = formula(qtl), method = "hk", covar = covar, verbose=F)
  return(newqtl)
}
