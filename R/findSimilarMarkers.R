#' @title Method to improve a genetic map.
#'
#' @description
#' \code{findSimilarMarkers} Find markers that have a small recombination fraction and
#' drop the one with combined greater segregation distortion and/or missing data. ***Note:
#' if the cross object has many markers (>1000), avoid running this function on more than one
#' chromosome at a time. Also make sure to run est.rf first and use re.est.map = FALSE***
#'
#' @param cross The qtl cross object to search
#' @param chr The chromosome to scan. Can be a vector of chromosome names or a single name.
#' @param rf.threshold The recombination fraction threshold to drop a marker. If est.rf has
#' not been run on cross, it will be done so automatically. See qtl::est.rf for details
#' @param sd.weight The weighting of segregation distortion rank in dropping a marker.
#' Higher values relative to na.weight increase the weight of the sd rank. Setting a value
#' of 0 removes sd as a factor in choosing the best marker.
#' @param na.weight Same as sd.weight, but for the number of nas.
#' @param drop.similar.markers Should the similar markers be dropped and a subsetted cross
#' be returned?
#' @param re.est.map Should the map be re-estimated after drop.similar.markers?
#' @param ... Additional arguments passed to est.map if re.est.map = TRUE
#' @return Either a character vector of the markers that should be dropped (if drop.similar.markers = FALSE),
#' or a new cross object with the similar markers dropped.
#'
#' @examples
#' set.seed(42)
#' map<-sim.map(len = c(50,20), n.mar = c(10,20), include.x=F)
#' cross0<-sim.cross(map, n.ind=50, type="f2", map.function="kosambi", error.prob=.01, missing.prob = .05)
#' cross0<-est.rf(cross0)
#' cross1<-findSimilarMarkers(cross0, error.prob=0.001, map.function="kosambi", rf.threshold = 0.005)
#' cross2<-findSimilarMarkers(cross0, error.prob=0.001, map.function="kosambi", rf.threshold = 0.005, keepEnds=TRUE)
#' par(mfrow=c(2,1))
#' plot.map(cross0, cross1, main = "comparison of full and culled maps")
#' plot.map(cross0, cross2, main = "comparison of full and culled maps")
#'
#' @import qtl
#' @export

findSimilarMarkers<-function(cross,
                             chr = NULL,
                             rf.threshold=0.02,
                             sd.weight=1,
                             na.weight=1,
                             drop.similar.markers = TRUE,
                             re.est.map = TRUE,
                             verbose = TRUE,
                             keepEnds = FALSE,
                             ...){
  if(is.null(chr)) chr <- chrnames(cross)
  if(verbose) cat("scanning chr: ", paste(chr, collapse = ","),"\n")
  todrop<-lapply(chr, function(i){
    gt<-geno.table(cross, chr=i)
    ms<-markernames(cross, chr =i)
    rf<-pull.rf(cross, chr = i, what = "rf")
    if(keepEnds){
      rf<-rf[-1,-1]
      rf<-rf[-nrow(rf),-ncol(rf)]
    }
    rf[!upper.tri(rf)]<-1
    bads<-vector()
    if(min(rf)>=rf.threshold){
      return(NULL)
    }else{
      while(min(rf)<rf.threshold){
        worst<-colnames(rf)[which(rf == min(rf, na.rm=TRUE), arr.ind=T)[1,]]
        gtm<-gt[worst,]
        ord.missing<-rank(gtm$missing)
        ord.sd<-rank(-gtm$P.value)
        prod<-(ord.missing*sd.weight)+(ord.sd*sd.weight)
        todrop<-rownames(gtm)[-which.min(prod)]
        bads<-c(bads, todrop)
        which.todrop<-which(colnames(rf) == todrop)
        rf<-rf[-which.todrop,-which.todrop]
      }
      return(bads)
    }
  })
  todrop<-unlist(todrop)
  if(drop.similar.markers){
    cross<-drop.markers(cross, markers = todrop)
    if(re.est.map){
      if(verbose) cat("re-estimating genetic map\n")
      em<-est.map(cross, ...)
      cross<-replace.map(cross, map = em)
    }
    return(cross)
  }else{
    return(todrop)
  }
}
