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
#'
#' cis<-calcCis(mod)
#' segmentsOnPeaks(cross, mod = mod)
#' plotLodProfile(mod, showallchr=F)
#' segmentsOnPeaks(cross, qtlModel = mod, calcCisOutput = cis,  showallchr=F)
#'
#' plotLodProfile(mod, showallchr=T)
#' segmentsOnPeaks(cross, qtlModel = mod, calcCisOutput = cis,  showallchr=T)
#' segmentsOnPeaks(cross, mod = mod, showallchr=F, add = T, drop = 3, col = "purple", int.y=.5)
#'
#' s1<-scanone(cross, method="hk", pheno.col="pheno1")
#' perm<-scanone(cross, n.perm=100, method="hk",pheno.col="pheno1", verbose=F)
#' cis<-calcCis(s1.output=s1, perm.output=perm)
#' plot(s1)
#' segmentsOnPeaks(cross, s1.output=s1, calcCisOutput = cis)
#' plot(s1, chr = c(2,5))
#' segmentsOnPeaks(cross, s1.output=s1, calcCisOutput = cis, chr = c(2,5))
#'
#' @import qtl
#' @export
#'

segmentsOnPeaks<-function(cross, calcCisOutput = NULL,
                          ci.chr = NULL, peak = NULL, l = NULL, h = NULL,
                          s1.output = NULL, qtlModel = NULL,
                          showallchr = FALSE, chr = NULL,
                          pt.pch = 8, pt.cex = .8, pt.col = "red", int.y = 0, ...){

  if(is.null(s1.output) & is.null(qtlModel) |
     !is.null(s1.output) & !is.null(qtlModel))
    stop("either s1.output or qtlModel must be provided")

  if(!is.null(qtlModel)){
    if(showallchr){
      chrs<-chrnames(cross)
    }else{
      chrs<-mod$chr
    }
  }else{
    if(is.null(chr)){
      chrs<-chrnames(cross)
    }else{
      chrs<-chr
      s1.output<-s1.output[s1.output$chr %in% chrs, ]
      cis<-cis[cis$chr %in% chrs, ]
    }
  }

  m<-pull.map(cross, chr=chrs, as.table=T)
  class(m)<-c("scanone", "data.frame")
  if(is.null(calcCisOutput) & is.null(ci.chr)
     & is.null(peak) & is.null(l) & is.null(h))
    stop("if calcCisOutput is not provided, ci.chr, peak, l and h must be specified")
  if(!is.null(calcCisOutput)){
    cis<-calcCisOutput
  }else{
    if(length(unique(c(length(ci.chr),
                       length(peak),
                       length(l),
                       length(h)))) != 1)
      stop("ci.chr, peak, l and h must be the same length")
    cis<-data.frame(chr = ci.chr, pos = peak,
                    lowposition = l, highposition = h)
  }
  for(i in 1:nrow(cis)){
    points(xaxisloc.scanone(m, cis$chr[i], cis$pos[i]), int.y,
           pch = pt.pch, cex=pt.cex, col=pt.col)
    segments(xaxisloc.scanone(m, cis$chr[i], cis$lowposition[i]), int.y,
             xaxisloc.scanone(m, cis$chr[i], cis$highposition[i]), int.y,
             ...)
  }
}
