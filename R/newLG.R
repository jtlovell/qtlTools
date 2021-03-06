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
#' fake.f2<-est.rf(fake.f2)
#' cross<-fake.f2
#' #Perturb the marker order and chromosome names
#' markerlist<-lapply(chrnames(cross), function(x) markernames(cross, chr =x))
#' names(markerlist) = sample(chrnames(cross),replace = F)
#' markerList<-lapply(markerlist, function(x) sample(x, replace = F))
#' #Make new cross with new marker order and chromosome names
#' cross2<-newLG(cross, markerList = markerList)
#'
#' # completely perturb the genotype matrix, then pull back together
#' markerlist = list("1" = sample(markernames(cross), replace = F))
#' cross.noOrder<-newLG(cross, markerList = markerlist)
#' lgmar <- formLinkageGroups(cross.noOrder)
#' print(tab<-table(lgmar))
#' marlist<-lapply(1:length(unique(lgmar$LG)),
#'   function(x) rownames(lgmar)[lgmar$LG==x])
#' cross2<-newLG(cross = cross.noOrder,
#'               markerList = marlist)
#' }
#' @import qtl
#' @export
newLG<-function(cross, markerList, keep.rf = FALSE){
  if(any(is.null(names(markerList)))){
    names(markerList)<-as.character(1:length(markerList))
  }
  # drop markers not in markerList
  newmars<-unlist(markerList)
  if(any(duplicated(newmars))){
    stop("duplicated markers found in markerList, all markers must be unique\n")
  }

  if (!("rf" %in% names(cross)) & keep.rf) {
    warning("Running est.rf.")
    cross <- est.rf(cross)
  }

  if(any(!newmars %in% markernames(cross))){
    stop("some markers in list are not in the cross, dropthem\n")
  }

  if(any(!markernames(cross) %in% newmars)){
    todrop<-markernames(cross)[!markernames(cross) %in% newmars]
    cross<-drop.markers(cross, markers = todrop)
  }

  n.mar <- nmar(cross)
  tot.mar <- totmar(cross)
  if(keep.rf) {
    rf <- cross$rf
    diagrf <- diag(rf)
    if (ncol(rf) != tot.mar)
      stop("dimension of recombination fractions inconsistent with no. markers in cross.")
    onlylod <- attr(cross$rf, "onlylod")

    lod <- rf
    lod[lower.tri(rf)] <- t(rf)[lower.tri(rf)]
    rf[upper.tri(rf)] <- t(rf)[upper.tri(rf)]
    diag(rf) <- 1
    diag(lod) <- 0
  }

  marnam <- markernames(cross)
  chrstart <- rep(names(cross$geno), n.mar)
  ingrp <- 1:tot.mar
  chrnum<-1:length(markerList)
  revgrp <- rep(chrnum,sapply(markerList, length))

  cross <- clean(cross)
  chrtype <- rep(sapply(cross$geno, class), n.mar)
  crosstype <- class(cross)[1]
  g <- pull.geno(cross)
  cross$geno <- vector("list", max(revgrp))
  names(cross$geno) <- 1:max(revgrp)

  for (i in 1:max(revgrp)) {
    cross$geno[[i]]$data <- g[, markerList[[i]], drop = FALSE]
    cross$geno[[i]]$map <- seq(0, by = 10, length = length(markerList[[i]]))
    if (crosstype == "4way") {
      cross$geno[[i]]$map <- rbind(cross$geno[[i]]$map,
                                   cross$geno[[i]]$map)
      colnames(cross$geno[[i]]$map) <- colnames(cross$geno[[i]]$data)
    }else{
      names(cross$geno[[i]]$map) <- colnames(cross$geno[[i]]$data)
    }
    thechrtype <- unique(chrtype[revgrp == i])
    if (length(thechrtype) > 1){
      warning("Problem with linkage group ", i, ": A or X?\\n",
              paste(thechrtype, collapse = " "))
    }else{
      class(cross$geno[[i]]) <- thechrtype
    }
  }
  mname <- markernames(cross)
  m <- match(mname, marnam)
  if(keep.rf) {
    rf <- rf[m, m]
    lod <- lod[m, m]
    rf[upper.tri(rf)] <- lod[upper.tri(lod)]
    diag(rf) <- diagrf[m]
    cross$rf <- rf
  }
  return(cross)
}
