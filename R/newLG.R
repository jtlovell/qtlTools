#' @title Re-order markers and re-form linkage groups
#'
#' @description
#' \code{newLG} Using a named list of markers (list element names = LG names),
#' re-assign markers to linkage groups and re-orient markers. Most of the code is
#' taken from qtl::formLinkageGroups.#'
#' Note: all intermediate calculations are dropped. All chromosomes are set as
#' autosomes.
#'
#' @param cross The QTL cross object.
#' @param markerList A named list where element names are the chromosome IDs. Each
#' element must contain a character vector specifying marker names. Markers in
#' cross, but not in markerList will be dropped. Any markers in markerList that are not
#' in markernames(cross) will be dropped.
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
#' data(fake.bc)
#' cross<-fake.bc
#' \dontrun{
#' ... more here ...
#' }
#' @import qtl
#' @export
newLG<-function(cross, markerList){
  markerList<-lapply(markerList, function(x) x[x %in% markernames(cross)])
  
  if ("rf" %in% names(cross)){
    has.rf<-TRUE
    mars<-unlist(markerList)
    o <- match(mars, colnames(cross$rf))
    rf<-cross$rf[o,o]
  }else{
    has.rf<-FALSE
  }
  cross<-clean(cross)
  n.mar<-nmar(cross)
  crosstype <- class(cross)[1]
  g <- pull.geno(cross)
  cross$geno <- vector("list", length(markerList))
  names(cross$geno) <- names(markerList)
  for (i in names(markerList)) {
    mars<-markerList[[i]]
    cross$geno[[i]]$data <- g[, mars, drop = FALSE]
    cross$geno[[i]]$map <- seq(0, by = 10, length = length(mars))
    if (crosstype == "4way") {
      cross$geno[[i]]$map <- rbind(cross$geno[[i]]$map,
                                   cross$geno[[i]]$map)
      colnames(cross$geno[[i]]$map) <- colnames(cross$geno[[i]]$data)
    }else{
      names(cross$geno[[i]]$map) <- colnames(cross$geno[[i]]$data)
    }
    class(cross$geno[[i]]) <- "A"
  }
  
  if(has.rf) {
    cross$rf <- rf
  }
  return(cross)
}
