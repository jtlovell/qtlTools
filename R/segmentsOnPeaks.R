#' @title Single QTL confidence intervals on LOD curves
#'
#' @description
#' \code{segmentsOnPeaks} Overlay confidence intervals onto a LOD plot,
#' from as either plotLodProfile or plot.scanone.
#' @param cross R/qtl cross object
#' @param calcCisOutput The output from qtlTools::calcCis. If this is not
#' supplied, ci.chr, peak, l, h must be. If calcCisOutput, these arguments are
#' ignored
#' @param ci.chr vector of chromosomes for segments
#' @param peak vector of QTL peak positions (cM) for segments
#' @param l vector of QTL lower confidence interval bounds (cM) for segments
#' @param h vector of QTL upper confidence interval bounds (cM) for segments
#' @param s1.output If plotting on a scanone LOD plot, provide the scanone object
#' used previously.
#' @param chr If plotting on a specified set of chromosomes via scanone, replicate
#' the chr call here.
#' @param mod If plotting on lodProfile, pass the qtl model object that was input into
#' plotLodProfile
#' @param showallchr Must match call to plotLodProfile if qtl is specified
#' @param pt.pch pch parameter to points
#' @param pt.cex cex parameter to points
#' @param pt.col color parameter to points
#' @param int.y the y intercept of the segments
#' @param add should segements be added to an existing plot?
#' @param ... additional parameters passed to segments, e.g. col, lty, etc.
#'
#' @examples
#' \dontrun{
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
#' plotLodProfile(mod, showallchr=FALSE)
#' segmentsOnPeaks(cross, qtl = mod, calcCisOutput = cis,  showallchr=FALSE)
#'
#' plotLodProfile(mod, showallchr=TRUE)
#' segmentsOnPeaks(cross, qtl = mod, calcCisOutput = cis, showallchr=TRUE)
#' cis<-calcCis(mod, drop = 3)
#' segmentsOnPeaks(cross, qtl = mod, calcCisOutput = cis, showallchr=TRUE,
#'    col = "purple", int.y=.5)
#'
#' s1<-scanone(cross, method="hk", pheno.col="pheno1")
#' perm<-scanone(cross, n.perm=100, method="hk",pheno.col="pheno1", verbose=FALSE)
#' cis<-calcCis(s1.output=s1, perm.output=perm)
#' plot(s1)
#' segmentsOnPeaks(cross, s1.output=s1, calcCisOutput = cis)
#' plot(s1, chr = c(2,5))
#' segmentsOnPeaks(cross, s1.output=s1, calcCisOutput = cis, chr = c(2,5))
#' }
#' @import qtl
#' @export
#'

segmentsOnPeaks<-function(cross, calcCisOutput = NULL,
                          ci.chr = NULL, peak = NULL, l = NULL, h = NULL,
                          s1.output = NULL, mod = NULL,
                          showallchr = FALSE, chr = NULL,
                          pt.pch = 8, pt.cex = .8, pt.col = "red", int.y = 0, ...){

  if(is.null(s1.output) & is.null(mod) |
     !is.null(s1.output) & !is.null(mod))
    stop("either s1.output or qtl model (mod) must be provided")

  if(!is.null(mod)){
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
