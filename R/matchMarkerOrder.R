#' @title Orient chromosomes to match some previously found marker order.
#'
#' @description
#' \code{matchMarkerOrder} Use previously known marker orders (e.g. from a
#' annotation) to inform the marker orientation along a chromosome. Does not re-order
#' markers within a linkage group, but only flips orders of chromosomes that seem to be
#' inverted relative to the annotation.
#'
#' @param cross The QTL cross object. If marker names
#' are formated as physicalChromosome_physicalBasePairPosition ("_" separated), this
#' is the only required argument. *Note, the physical Chromosome must match exactly
#' the chromosome names of the cross.*
#' If there are additional characters (e.g. "Chr01" = "1"), specify the character
#' string to sub stripped in marker.char.rm. In this example, marker.char.rm = "Chr0".
#' @param physicalMarkerOrder If marker names are not as specified above, this is a vector of
#' marker names in the order they should exist in the physical genome.
#'
#' @details For each chromosome, fits a linear model and asks if the slope < 0. If so,
#' flips the maker order on that chromosome.
#'
#' @return A cross object with flipped marker order on chromosomes that seem to be in the
#' reverse orientation of the annotation.
#'
#' @examples
#' \dontrun{
#' library(qtlTools)
#' set.seed(42)
#' map<-sim.map(len = c(100,200,50), n.mar = c(50,100,25),
#' include.x = F, eq.spacing=T)
#' cross<-sim.cross(map, type = "riself",
#' map.function = "kosambi", error.prob = 0)
#' cross<-est.rf(cross)
#' markerlist<-lapply(chrnames(cross), function(x)
#' sample(markernames(cross, chr =x)))
#' names(markerlist)<-as.character(chrnames(cross))
#' cross.rand<-newLG(cross, markerList = markerlist)
#' cross.ord<-tspOrder(cross = cross.rand,
#'                     max.rf = .5,
#'                     concorde_path = "/Users/John/Documents/concorde/TSP")
#' cross.ord<-replace.map(cross.ord, est.map(cross.ord, map.function = "kosambi"))
#' orig.marker.order<-lapply(chrnames(cross), function(x)
#' markernames(cross, chr = x))
#' names(orig.marker.order)<-chrnames(cross)
#'
#' cross.match<-matchMarkerOrder(cross = cross.ord,
#'   physicalMarkerOrder = orig.marker.order)
#' }
#' @import qtl
#' @export

matchMarkerOrder<-function(cross, physicalMarkerOrder = NULL, plotit = T){

  if(is.null(physicalMarkerOrder)){
    chr<-splitText(markernames(cross))
    pos<-as.numeric(splitText(markernames(cross),num = 2))
    physicalMarkerOrder<-markernames(cross)[order(chr,pos)]
  }
  if(plotit){
    plot(match(markernames(cross), physicalMarkerOrder),1:length(physicalMarkerOrder),
         xlab = "order of markers", ylab = "physical order of markers",
         pch = 16, col = "grey", cex = .5)
  }
  flipIt<-sapply(chrnames(cross), function(x){
    marker.id = markernames(cross, chr = x)
    phys.order = match(markernames(cross, chr = x),physicalMarkerOrder)
    orig.order = 1:length(markernames(cross, chr = x))
    out<-lm(phys.order~orig.order)$coefficients["orig.order"]
    return(out<0)
  })
  toflip<-chrnames(cross)[flipIt]
  for(i in toflip) cross<-flip.order(cross, chr = i)

  if(plotit){
    points(match(markernames(cross), physicalMarkerOrder),1:length(physicalMarkerOrder),
         xlab = "order of markers", ylab = "physical order of markers",
         pch = 16, col = "black", cex = .2)
    legend("bottomright",c("initial order","corrected order"),
           pch = c(16,16), col = c("grey","black"), pt.cex = c(.5,.2))
  }
  return(cross)
}
