#' @title Find markers that have inferred positions within gaps in a map
#'
#' @description
#' \code{fillGapsInMap} Using known cM positions and chromosomes,
#' determine which markers could be used to fill garps in a map.
#'
#' @param cross The QTL cross object.
#' @param minGapSize The smallest gap size to check for markers. Gaps in the map
#' smaller than this will be ignored.
#' @param marker.names The names of markers to check
#' @param marker.chr The inferred chromsomes of the markers - must match chrnames of
#' the cross.
#' @param marker.pos The inferred cM positions of the markers.
#' @param return4newLG Logical, Should the gaps be combined with the genetic map to pipe
#' straight into the newLG function?
#' @param ... Not currently in use.
#' @details A simple function to parse output from inferMarker postion
#'
#' @return If return4newLG = TRUE, a named and ordered list of all markers, including
#' those in the map. This output can be input directly into the newLG argument 'markerList'.
#' If return4newLG = FALSE, returns a simple named list of the markers that could be added.
#' List element names match the chromosomes that they belong to.
#'
#' @examples
#' library(qtlTools)
#' data(fake.bc)
#' cross<-fake.bc
#' \dontrun{
#' ... more here ...
#' }
#' @import qtl
#' @export
fillGapsInMap<-function(cross,
                        minGapSize = 1,
                        marker.names,
                        marker.chr,
                        marker.pos,
                        return4newLG = F, ...){
  out<-lapply(chrnames(cross), function(x){
    map<-pullMap(cross, chr = x)
    gaps<-diff(map$pos)

    if(any(gaps>minGapSize)){
      wh<- which(gaps>minGapSize)
      gap.start<-map$pos[wh]
      gap.end<-map$pos[wh+1]
      out.chr<-sapply(1:length(wh), function(y){
        marker.names[marker.chr == x &
                       marker.pos>gap.start[y] &
                       marker.pos<gap.end[y] ]
      })
      return(as.character(out.chr))
    }
  })

  names(out)<-chrnames(cross)
  if(return4newLG){
    newmars<-unlist(out)
    outn<-data.frame(marker.name = marker.names,
                     chr = marker.chr,
                     pos = marker.pos, stringsAsFactors=F)
    outn<-outn[outn$marker.names %in% newmars,]
    outd<-rbind(pullMap(cross), outn)
    outd<-outd[with(outd, order(chr, pos)),]
    outd$chr<-as.character(outd$chr)
    out<-lapply(chrnames(cross), function(x)
      as.characterout(d$marker.name[d$chr == x]))
    names(out)<-chrnames(cross)
  }
  return(out)
}
