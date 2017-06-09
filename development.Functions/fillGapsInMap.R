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
                        newMarGap = 1,
                        marker.names,
                        marker.chr,
                        marker.pos,
                        marker.lod,
                        return4newLG = F,
                        which.position = 3,
                        ...){
  out<-lapply(chrnames(cross), function(x){
    map<-pullMap(cross, chr = x)
    colnames(map)[which.position]<-"pos"
    map<-map[,c("marker.name", "chr","pos")]
    gaps<-diff(map$pos)

    if(any(gaps>minGapSize)){
      wh<- which(gaps>minGapSize)
      gap.start<-map$pos[wh]
      gap.end<-map$pos[wh+1]
      out.chr<-lapply(1:length(wh), function(y){
        sq<-seq(from = ceiling(gap.start[y]),
            to = floor(gap.end[y]),
            by = newMarGap)
        sapply(1:(length(sq)-1), function(z){
          gsz<-sq[z]
          gez<-sq[z+1]
          
          mout<-marker.names[marker.chr == x &
                         marker.pos>gsz &
                         marker.pos<gez]
          mlods<-marker.names[marker.chr == x &
                                marker.pos>gsz &
                                marker.pos<gez]
          return(mout[which.max(mlods)])
        })
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
