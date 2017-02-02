#' @title Find the best mapping position for a marker
#'
#' @description
#' \code{inferMarkerPos} Used to orient scaffolds or find markers that can fill
#' gaps in a map.
#'
#' @param cross The QTL cross object.
#' @param marker.matrix A numeric genotype matrix where the column names are markers and
#' the row names are identical to the line IDs in the cross object
#' @param initial.geno.density The initial distance between markers. Higher values speeds up
#' the initial scan for the best chromosome. Corresponds to the min.distance argument
#' of repPickMarkerSubset.
#' @param final.step The density of pseudomarkers on the final scan. Corresponds to the
#' step argument of calc.genoprob. Smaller values get more precision on estimates of
#' marker position, but slow down the algorithm.
#' @param est.ci Logical, should confidence intervals around the inferred position be
#' estimated?
#' @param ... Additional arguments passed to scanone. Cannot be "chr" or "method" as
#' these are hard coded in the function.
#' @details There are many methods to find the most likely position of a new marker.
#' This approach attempts to speed up detection of the most likely position of a marker
#' by using the haley-knott qtl scan to find the strongest correlation between markers
#' in the genetic map and the numeric-coded genotypes of a new marker.
#' This function can be used independently to order scaffolds during genome assembly,
#' or to improve a genetic map.
#' In the latter case, the results can be piped into fillGapsInMap to find the markers
#' that will do the best job of increasing marker density. Once these markers are known
#' it is best to go back and reconstruct the genetic map. However, if only a few markers
#' are to be added, one could use qtl::addmarkers.
#'
#' @return A data frame with the most likely position of the marker.
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
inferMarkerPos<-function(cross, marker.matrix,
                         initial.geno.density = 20,
                         final.step = .1,
                         est.ci = FALSE,
                         ...){

  if(!"prob" %in% names(cross$geno[[1]]))
    stop("must run calc.genoprob first\n")

  error.prob<-attr(cross$geno[[1]]$prob,"error.prob")
  off.end<-attr(cross$geno[[1]]$prob,"off.end")
  map.function<-attr(cross$geno[[1]]$prob,"map.function")

  mcross<-cross

  id<-as.character(getid(cross))
  if(nind(cross) != nrow(marker.matrix))
    stop("the number of rows in marker.matrix != number of indviduals in the cross\n")
  if(!identical(rownames(marker.matrix), id))
    warning("marker rownames are not the same as cross IDs ...
            proceeding assuming marker rows correspond to cross IDs\n")
  if(is.null(colnames(marker.matrix))){
    warning(paste0("no marker names (colnames(marker.matrix)) provided ...
                   assigning names as m1:m",ncol(marker.matrix)))
    colnames(marker.matris)<-paste0("m",1:ncol(marker.matrix))
  }

  if(verbose) cat("preparing cross objects for initial and final scans\n")
  mnames<-colnames(marker.matrix)
  if(any(mcross$pheno %in% mnames)){
    mcross$pheno<-data.frame(mcross$pheno[,-which(mcross$pheno %in% mnames)], marker.matrix)
  }else{
    mcross$pheno<-data.frame(mcross$pheno, marker.matrix)
  }

  icross<-repPickMarkerSubset(cross, na.weight = 2, sd.weight = 1,
                              min.distance = initial.geno.density, verbose = FALSE)
  fcross<-clean(mcross)
  fcross<-calc.genoprob(icross, stepwidth = "fixed", step = final.step,
                        error.prob = error.prob, off.end = off.end, map.function = map.function)

  if(verbose) cat("running initial scanone to find best chromosome and position\n")
  s1s.i<-scanone(icross, pheno.col = mnames, method = "hk", ...)
  chrs<-s1s.i[,"chr"]
  pos<-s1s.i[,"pos"]
  wh<-apply(s1s.i[,mnames], 2, which.max)
  chr.max<-chrs[wh]
  pos.max<-pos[wh]

  if(verbose) cat("running final scanone to refine best position\n")
  out.f<-lapply(chrnames(fcross), function(x){
    if(verbose) cat("analyzing chr",x,"\n")
    phes<-names(chr.max)[chr.max==x]
    s1s.f<-scanone(fcross, pheno.col = phes, method = "hk", chr = x, ...)
    wh.f<-apply(s1s.f[,phes], 2, which.max)
    pos.max<-pos[wh.f]
    lod.max<-apply(s1s.f[,phes], 2, function(y) y[which.max(y)])
    out<-data.frame(marker.name = phes,
                    chr = x,
                    pos = pos.max,
                    lod = lod.max,
                    stringsAsFactors = F)

    if(est.ci){
      if(verbose) cat("calculating confidence intervals of best position\n")
      out.ci<-sapply(1:length(phes), function(y){
        ciout<-lodint(results = s1s.f, chr=x,lodcolumn=y,
                      drop=1.5, expandtomarkers=F)
        return(as.numeric(ciout[c(1,3),2]))
      })
      out<-data.frame(out,
                      lowCI = out.ci[,1],
                      highCI = out.ci[,2],
                      stringsAsFactors = F)
    }

    return(out)
  })

  inferred.pos<-do.call(rbind, out.f)
  inferred.pos<-inferred.pos[with(inferred.pos, order(chr, pos)),]
  return(inferred.pos)
}
