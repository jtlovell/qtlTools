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
#' @import qtl
#' @export

findSimilarMarkers<-function(cross,
                             chr = NULL,
                             rf.threshold=0.02,
                             sd.weight=1,
                             na.weight=1,
                             drop.similar.markers = TRUE,
                             re.est.map = TRUE,
                             ...){
  if(is.null(chr)) chr <- chrnames(cross)
  gt<-geno.table(cross, chr=chr)
  ms<-markernames(cross, chr =chr)
  rf<-pull.rf(cross, chr = chr, what = "rf")
  rf[lower.tri(rf)]<-NA
  nondup<-sapply(ms, function(x){
    wh<-names(which(rf[x,]<rf.threshold))
    if(length(wh>0)){
      gtm<-gt[wh,]
      ord.missing<-rank(gtm$missing)
      ord.sd<-rank(-gtm$P.value)
      prod<-(ord.missing*sd.weight)+(ord.sd*sd.weight)
      rownames(gtm)[which.min(prod)]
    }else{
      x
    }
  })
  nondup<-unique(nondup)
  nondup<-nondup[order(nondup)]
  todrop<-markernames(cross, chr = chr)[!markernames(cross, chr = chr) %in% nondup]
  if(drop.similar.markers){
    cross<-drop.markers(cross, markers = todrop)
    if(re.est.map){
      em<-est.map(cross, ...)
      cross<-replace.map(cross, map = em)
    }
    return(cross)
  }else{
    return(todrop)
  }
}
