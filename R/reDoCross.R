#' @title Run a set of functions that improve the genetic map
#'
#' @description
#' \code{reDoCross} Drop close markers, re-order markers via ripple, then re-check
#' for close markers. Returns a cross object with improved marker order and overly
#' dense marker regions removed.
#'
#' @param cross The QTL cross object.
#' @param window Passed to ripple. Smaller values (e.g. 3) run more quickly but explore
#' fewer possible marker orders. Larger values (e.g. 6) should only be used on maps with
#' few markers. Large numbers of markers with large windows will take a long time.
#' @param min.distance The minimum distance between markers in the final map. If two markers
#' are closer than this value, drop the one with fewer NAs and better segregation distortion.
#' @param map.function The map function to pass on to est.map
#' @param sex.sp Should sex-specific maps be estimated in a 4-way cross? Passed on to est.map.
#' @param verbose Should updates be printed ot the terminal?
#' @param initialEstMap Should the map be estimated before any functions are run?
#'
#' @details If any chromosomes are sex-specific, re-class such chromosomes after running
#' newLG: using class(cross$geno[["Xchrom"]]) <- "S" or whatever. The new map has
#' arbitrary cM positions where each marker is separated by 10cM. Run est.map to
#' get true cM positions.
#'
#' @return  A cross object with markers reordered.
#'
#' @examples
#' library(qtlTools)
#' data(fake.f2)
#' cross<-fake.f2
#' \dontrun{
#' map=est.map(cross, map.function = "haldane")
#' cross=replace.map(cross, map)
#' cross2<-reDoCross(cross, map.function = "haldane")
#' plot.map(cross, cross2)
#' }
#'
#' @import qtl
#' @export
reDoCross<-function(cross,  window = 5, min.distance = 1,
                    map.function = "kosambi", sex.sp=F,
                    verbose=T, initialEstMap = TRUE,...){
  loadNamespace("qtl")
  if(initialEstMap){
    if(verbose) cat("initial estimation of map \n")
    map<-est.map(cross, map.function = map.function, sex.sp=sex.sp, ...)
    cross<-replace.map(cross, map)
  }

  orig.len<-sum(chrlen(cross))
  orig.nmar<-sum(nmar(cross))

  if(verbose) cat("drop close markers")
  cross<-repPickMarkerSubset(cross = cross, verbose=F, min.distance = min.distance)
  if(verbose) cat(" --- tossed", orig.nmar-sum(nmar(cross)),"markers\n")

  if(verbose) cat("ripple marker order")
  cross<-repRipple(cross, window = window, verbose = T, ripVerb=T)
  if(verbose) cat(" --- reduced map size by", orig.len-sum(chrlen(cross)),"\n")

  orig.len<-sum(chrlen(cross))
  orig.nmar<-sum(nmar(cross))

  if(verbose) cat("drop close markers")
  cross<-repPickMarkerSubset(cross = cross, verbose=F, min.distance = min.distance)
  if(verbose) cat(" --- tossed", orig.nmar-sum(nmar(cross)),"markers\n")

  if(verbose) cat("final map estimation \n")
  map<-est.map(cross, map.function = map.function, sex.sp=sex.sp, ...)
  cross<-replace.map(cross, map)
  return(cross)
}
