#' @title Order markers within linkage groups
#'
#' @description
#' \code{tspOrder} Apply a traveling salesperson problem solver to find the
#' shortest path through the recombination fraction matrix. Requires that
#' the TSP R package and the TSP solver "Concorde" program are installed. See details.
#'
#' @param cross The QTL cross object.
#' @param concorde_path Required. The directory containing concorde executables.
#' @param return If "cross" is specified, pass the marker order through newLG, which quickly
#' reorders the markers of the original cross to follow the TSP order. Otherwise,
#' return a named list of markers in order for each chromosome.
#' @param ... Additional arguments passed on to solve_TSP that are then paseed to TSP::contol.
#' @details This function relies on the TSP R packages to perform the TSP solvers. See documentation
#' therein, esspecially the function TSP::solve_TSP. We permit inference of marker order within a
#' Hamiltonian circuit by adding a dummy node that has 0 distance to all other nodes. This
#' allows for discrete start and end points and is more appropriate for genetic map construction
#' than forcing a complete route through all markers.
#'
#' To install the concorde program:
#' http://www.math.uwaterloo.ca/tsp/concorde/downloads/codes/src/co031219.tgz.
#' If on a mac, this is not super easy.
#' See https://qmha.wordpress.com/2015/08/20/installing-concorde-on-mac-os-x/
#' for details on the best way to do this.
#'
#' @return Either a cross object with reordered markers, or a named list
#' of markers in a new order.
#'
#' @examples
#' \dontrun{
#' library(qtlTools)
#' set.seed(42)
#' map<-sim.map(len = 200, n.mar = 200, include.x = F, eq.spacing=T)
#' cross<-sim.cross(map, type = "riself", map.function = "kosambi", error.prob = 0)
#' cross<-est.rf(cross)
#' markerlist<-lapply(chrnames(cross), function(x) sample(markernames(cross, chr =x)))
#' names(markerlist)<-as.character(chrnames(cross))
#' cross.rand<-newLG(cross, markerList = markerlist)
#' set.seed(42)
#' cross.ord<-tspOrder(cross = cross.rand,
#'   max.rf = .5,
#'   concorde_path = "/Users/John/Documents/concorde/TSP")
#'
#' plot(match(markernames(cross), markernames(cross.ord)),
#' xlab = "position in similated map", ylab = "position in ordered map")
#' }
#' @import qtl
#' @import TSP
#' @export
tspOrder<-function(cross,
                   concorde_path,
                   return = "cross"){
  if(!requireNamespace("TSP", quietly = TRUE)){
    stop("install the TSP package to use tspOrder\n")
  }else{
    requireNamespace("TSP", quietly = TRUE)
  }
  if(is.null(concorde_path)){
    stop("if method = concorde, concorde_path must specify directory of program\n")
  }

  concorde_path(concorde_path)

  rf<-data.matrix(pull.rf(cross, what = "rf"))
  class(rf)<-"matrix"
  diag(rf)<-0

  markerList<-lapply(chrnames(cross), function(x) markernames(cross, chr = x))
  names(markerList)<-chrnames(cross)

  chr.ord<-lapply(names(markerList),function(x){
    mnames<-markerList[[x]]
    totsp<-rf[mnames,mnames]
    chr.tsp<-TSP(totsp)

    tsp <- insert_dummy(chr.tsp, label = "cut")
    tour <- solve_TSP(tsp, method="concorde")
    path <- cut_tour(tour, "cut")

    ord<-labels(path)
    return(ord)
  })
  names(chr.ord)<-names(markerList)
  if(return != "cross"){
    return(chr.ord)
  }else{
    out<-newLG(cross = cross, markerList = chr.ord, keep.rf = T)
    return(out)
  }
}
