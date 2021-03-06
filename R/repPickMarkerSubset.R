#' @title Method to improve an estimated genetic map.
#'
#' @description
#' \code{repPickMarkerSubset} finds the best marker within a cM window.
#'
#' @param cross The qtl cross object to search
#' @param chr The chromosome to scan. Can be a vector of chromosome names or a single name.
#' If NULL, run on all chromosomes.
#' @param sd.weight The weighting of segregation distortion rank in dropping a marker.
#' Higher values relative to na.weight increase the weight of the sd rank. Setting a value
#' of 0 removes sd as a factor in choosing the best marker. Only used if cross is not a 4way.
#' @param na.weight Same as sd.weight, but for the number of NAs.
#' Only used if cross is not a 4way.
#' @param min.distance Logical, should markers on the ends of the chromosomes always be retained?
#' @param verbose Logical, should updates be printed?
#'
#' @return A new cross object with the similar markers dropped.
#'
#' @examples
#' set.seed(42)
#' library(qtlTools)
#' map<-sim.map(len = c(50,20), n.mar = c(20,30), include.x=FALSE)
#' cross0<-sim.cross(map, n.ind=50, type="f2", map.function="kosambi",
#'    error.prob=.01, missing.prob = .05)
#' cross0<-est.rf(cross0)
#' cross1<-repPickMarkerSubset(cross0)
#' par(mfrow=c(2,1))
#' plot.map(cross0, cross1, main = "comparison of full and culled maps")
#'
#' @import qtl
#' @export
repPickMarkerSubset<-function(cross, chr = NULL,
                               na.weight = 2, sd.weight=1,
                               min.distance = 1, verbose=TRUE){
  loadNamespace("qtl")
  if(is.null(chr)) chr<-chrnames(cross)
  smallmarkers<-unlist(
    sapply(chr, function(i){
      gt<-geno.table(cross, chr = i)

      if(class(cross)[1]=="4way"){
        gt.names<-colnames(gt)[3:6]
        gt$P.value<-apply(gt[,gt.names],1,function(x)
          min(x, na.rm = T)/nind(cross))
      }
      gt$rank.p<-with(gt, rank(rank(-P.value, ties.method = "min")*sd.weight))
      gt$rank.sd<-with(gt, rank(rank(missing, ties.method = "min")*na.weight))
      wts<-with(gt, rank(rank.p+rank.sd))
      if(class(cross)[1]=="4way"){
        pickMarkerSubset(pull.map(cross)[[i]][1,], min.distance = min.distance, wts)
      }else{
        pickMarkerSubset(pull.map(cross)[[i]], min.distance = min.distance, wts)
      }
  }))
  dupmar<-markernames(cross)[!markernames(cross) %in% smallmarkers]
  cross<-drop.markers(cross, dupmar)
  if(verbose) cat("dropped ", length(dupmar), " markers")
  return(cross)
}
