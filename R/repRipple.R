#' @title Method to improve a genetic map.
#'
#' @description
#' \code{repRipple} iteratively runs ripple,
#' returning a cross object with the best order.
#'
#' @param cross The qtl cross object to search
#' @param chr The chromosome to scan. Can be a vector of chromosome names or a single name.
#' @param window The window size passed on to ripple
#' @param method The method passed on to ripple
#' @param re.est.map Should the map be re-estimated after drop.similar.markers?
#' @param map.function If re.est.map = TRUE, what map function should be used?
#' @param error.prob If re.est.map = TRUE, what error probability should be used?
#' @param verbose logical, should updates be reported?
#' @param ... Additional arguments passed to qtl::ripple
#' @return A new cross object with the optimal rippled marker order
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' map<-sim.map(len = c(50,20), n.mar = c(20,40), include.x=FALSE)
#' cross0<-sim.cross(map, n.ind=50, type="f2", map.function="kosambi")
#' ## Mess up the order on chr1
#' mars<-1:nmar(cross0)[1]
#' set.seed(42)
#' badorder<-order(jitter(mars, factor = 10))
#' cross0<-switch.order(cross0, chr = 1, order = badorder,
#'                      error.prob=0.001, map.function="kosambi")
#' newmap<-est.map(cross0, error.prob=0.001, map.function="kosambi")
#' cross0 <- replace.map(cross0, newmap)
#' cross0<-est.rf(cross0)
#' cross1<-repRipple(cross0, chr = 1, error.prob=0.001, map.function="kosambi",window = 6)
#' plot.rf(cross0, chr = 1, main = "recombination fractions before ripple")
#' plot.rf(cross1, chr = 1, main = "recombination fractions after ripple")
#' }
#' @import qtl
#' @export
repRipple<-function(cross, chr = NULL, window = 5,
                    method = "countxo", verbose = T,
                    map.function = "kosambi", sex.sp=F, clean1st = FALSE, ripVerb = TRUE, ...){
  loadNamespace("qtl")
  if(clean1st) cross<-clean(cross)
  if(is.null(chr)){
    chr<-chrnames(cross)
  }
  for(j in chr){
    if(length(markernames(cross, chr = j))>=3){
      if(verbose) cat(j,"...")
      new.xo<-0
      orig.xo<-1
      while(new.xo<orig.xo){
        mars<-lapply(chrnames(cross), function(x) markernames(cross,chr = x))
        names(mars)<-chrnames(cross)
        verb<-ifelse(new.xo == 0 & ripVerb , TRUE, FALSE)

        s<-summary(ripple(cross, chr = j, window = window,
                          method = method, verbose = verb))
        best<-as.numeric(which.min(s[,ncol(s)])[1])
        ord<-s[best,-ncol(s)]
        orig.xo<-s[1,ncol(s)]
        new.xo<- s[best,ncol(s)]
        if(verbose){
          if(orig.xo == new.xo){
            cat("no reduction in XOs found\n")
          }else{
            cat("orig n XO = ", orig.xo, "new n XO = ",new.xo,"\n")
          }
        }
        mars[[j]]<-mars[[j]][ord]
        cross<-newLG(cross = cross, markerList = mars)
      }
    }
  }

  if(verbose) cat("final map estimation")
  map<-est.map(cross, map.function = map.function, sex.sp=sex.sp,...)
  cross<-replace.map(cross, map)
  return(cross)
}
