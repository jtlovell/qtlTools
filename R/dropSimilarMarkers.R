#' @title Method to improve a genetic map.
#'
#' @description
#' \code{dropSimilarMarkers} finds markers that have a small recombination fraction and
#' drops the one with combined greater segregation distortion and/or missing data. ***Note:
#' if the cross object has many markers (>1000), avoid running this function on more than one
#' chromosome at a time. Also make sure to run est.rf first and use re.est.map = FALSE***
#'
#' @param cross The qtl cross object to search
#' @param chr The chromosome to scan. Can be a vector of chromosome names or a single name.
#' If NULL, run on all chromosomes.
#' @param rf.threshold The recombination fraction threshold to drop a marker. If est.rf has
#' not been run on cross, it will be done so automatically. See qtl::est.rf for details
#' @param sd.weight The weighting of segregation distortion rank in dropping a marker.
#' Higher values relative to na.weight increase the weight of the sd rank. Setting a value
#' of 0 removes sd as a factor in choosing the best marker.
#' @param na.weight Same as sd.weight, but for the number of NAs.
#' @param keepEnds Logical, should markers on the ends of the chromosomes always be retained?
#' @param doNotDrop Character vector of markers to retain no matter their rfs.
#' @param verbose Logical, should updates be printed?
#' @param ... if recombination fractions are not included in the cross object,
#' pass on additional arguments to est.rf.
#'
#' @return A new cross object with the similar markers dropped.
#'
#' @examples
#' set.seed(42)
#' map<-sim.map(len = c(50,20), n.mar = c(20,30), include.x=F)
#' cross0<-sim.cross(map, n.ind=50, type="f2", map.function="kosambi", error.prob=.01, missing.prob = .05)
#' cross0<-est.rf(cross0)
#' cross1<-dropSimilarMarkers(cross0)
#' cross2<-dropSimilarMarkers(cross0, keepEnds=T)
#' par(mfrow=c(2,1))
#' plot.map(cross0, cross1, main = "comparison of full and culled maps")
#' plot.map(cross0, cross2, main = "comparison of full and culled maps")
#'
#' @import qtl
#' @export

dropSimilarMarkers<-function(cross,
                             chr = NULL,
                             rf.threshold=0.02,
                             sd.weight=1,
                             na.weight=1,
                             keepEnds = FALSE,
                             doNotDrop = NULL,
                             verbose=TRUE,
                             ...){
  # 1. Get the rfs, geno table and chromosomes in order
  if(is.null(chr)) chr<-chrnames(cross)
  gt<-geno.table(cross, chr = chr)
  if(!"rf" %in% names(cross)){
    if(verbose) cat("running est.rf\n")
    cross<-est.rf(cross, chr = chr, ...)
  }
  rf<-pull.rf(cross, chr = chr, what = "rf")
  rf[!upper.tri(rf)]<-1

  # 2. drop the markers to retain from the matrix
  if(!is.null(doNotDrop)){
    if(verbose) cat("retaining markers: ",paste(doNotDrop, collapse=", "),"\n")
    dnd.index<-which(colnames(rf) %in% doNotDrop)
    rf<-rf[-dnd.index, -dnd.index]
  }
  if(keepEnds){
    if(verbose) cat("retaining first and last markers on each chromosome\n")
    tokeep<-as.character(
      unlist(
        lapply(pull.map(cross),function(x)
          c(names(x)[1], names(x)[length(x)]))))
    ends.index<-which(colnames(rf) %in% tokeep)
    rf<-rf[-ends.index, -ends.index]
  }

  # 3. Loop through the rf matrix, dropping one of the two markers with the lowest
  # recombination fraction.
  nmarstart<-sum(nmar(cross))
  if(verbose) cat("initial n markers: ",nmarstart,"\n")
  while(min(rf)<rf.threshold){
    worst<-colnames(rf)[which(rf == min(rf, na.rm=TRUE), arr.ind=T)[1,]]
    gtm<-gt[worst,]
    ord.missing<-rank(gtm$missing)
    ord.sd<-rank(-gtm$P.value)
    prod<-(ord.missing*sd.weight)+(ord.sd*sd.weight)
    badmars<-rownames(gtm)[-which.min(prod)]

    which.todrop<-which(colnames(rf) == badmars)
    rf<-rf[-which.todrop,-which.todrop]

    cross<-drop.markers(cross, markers = badmars)
  }
  nmarend<-sum(nmar(cross))
  if(verbose) cat("final n markers: ",nmarend,"\n")
  return(cross)
}
