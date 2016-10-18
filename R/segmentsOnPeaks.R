#' @title Single QTL confidence intervals on LOD curves
#'
#' @description
#' \code{segmentsOnPeaks} Overlay confidence intervals onto a LOD plot,
#' from as either plotLodProfile or plot.scanone.
#' @param cross R/qtl cross object
#' @param s1.output Scanone output to be passed to calcCis
#' @param perms.output Permutation output to be passed to calcCis
#' @param chr If passing a scanone analysis, what chromosomes should plot be subset to?
#' Defaults to all chromosomes
#' @param mod If plotting qtl intervals from a QTL model, pass model here.
#' Ensure that refineqtl, or stepwise w/ keep.lodprofile = T has been run
#' so that plotLodProfile can be called successfully.
#' @param showallchr If passing a model, should all chromosomes be plotted,
#' even if they do not contain a QTL?
#' @param qtlnames passed to calcCis
#' @param lodint passed to calcCis, should droplod intervals be calculated?
#' @param drop passed to calcCis, if lodint=TRUE, what drop level?
#' @param prob passed to calcCis, if lodint=FALSE, what baye probability level?
#' @param expandtomarkers passed to calcCis, should interval be expanded to nearest
#' marker?
#' @param pt.pch pch parameter to points
#' @param pt.cex cex parameter to points
#' @param pt.col color parameter to points
#' @param int.y the y intercept of the segments
#' @param add should segements be added to an existing plot?
#' @param ... additional parameters passed to segments, e.g. col, lty, etc.
#'
#' @examples
#' library(qtlTools)
#' data(fake.bc)
#' cross<-fake.bc
#' cross <- calc.genoprob(cross, step=2.5)
#'
#' # Make a QTL object and formula
#' mod <- makeqtl(cross, chr = c(2,5), pos = c(40,25), what = "prob")
#' sex <- data.frame(sex = cross$phe$sex)
#' nform <- "y ~ Q1 + Q2 + Q1*sex + sex"
#' #Calculate lodprofiles for confidence interval estimation
#' mod <- refineqtl(cross, mod, pheno.col = "pheno1",
#'                  qtl = mod, formula = nform, covar = sex, method="hk")
#' segmentsOnPeaks(cross, mod = mod)
#' segmentsOnPeaks(cross, mod = mod, showallchr=F)
#' segmentsOnPeaks(cross, mod = mod, showallchr=F, add = T, drop = 3, col = "purple", int.y=.5)
#'
#' s1<-scanone(cross, method="hk", pheno.col="pheno1")
#' perm<-scanone(cross, n.perm=100, method="hk",pheno.col="pheno1", verbose=F)
#'
#' segmentsOnPeaks(cross, s1.output=s1, perm.output=perm)
#' segmentsOnPeaks(cross, s1.output=s1, perm.output=perm, chr = c(2,5))
#'
#' @import qtl
#' @export
#'

segmentsOnPeaks<-function(cross, s1.output, perm.output, chr = NULL,
                          mod = NULL, showallchr=TRUE,
                          qtlnames = NULL, lodint = TRUE, drop = 1.5,
                          prob = 0.95, expandtomarkers = FALSE,
                          pt.pch = 8, pt.cex = .8, pt.col = "red", int.y = 0, add=FALSE, ...){
  if(!is.null(mod)){
    if(!add){
      plotLodProfile(mod, showallchr=showallchr)
    }
    if(showallchr){
      chrs<-chrnames(cross)
    }else{
      chrs<-mod$chr
    }

    cis<-calcCis(mod=mod, qtlnames=qtlname, lodint=lodint,
                 drop=drop, prob=prob, expandtomarkers=expandtomarkers)
  }else{
    if(ncol(s1.output)>3)
      stop("only one phenotype can be tested at a time")
    if(is.null(chr)){
      chrs<-chrnames(cross)
    }else{
      chrs<-chr
      s1.output<-s1.output[s1.output$chr %in% chrs, ]
      cis<-cis[cis$chr %in% chrs, ]
    }
    cis<-calcCis(s1.output=s1.output, perm.output=perm.output,
                 lodint=lodint, drop=drop, prob=prob,
                 expandtomarkers=expandtomarkers)

    if(!add){
      plot(s1.output)
    }
  }

  m<-pull.map(cross, chr=chrs, as.table=T)
  class(m)<-c("scanone", "data.frame")

  for(i in 1:nrow(cis)){
    points(xaxisloc.scanone(m, cis$chr[i], cis$pos[i]), int.y,
           pch = pt.pch, cex=pt.cex, col=pt.col)
    segments(xaxisloc.scanone(m, cis$chr[i], cis$lowposition[i]), int.y,
             xaxisloc.scanone(m, cis$chr[i], cis$highposition[i]), int.y,
             ...)
  }
}
