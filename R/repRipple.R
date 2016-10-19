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
#' @param makePlots Should the recombination fraction be iteratively plotted?
#' @param repeatloop Should ripple be run once, or iteratively?
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
repRipple<-function(cross, chr = NULL, window=6,
                    repeatloop = TRUE, makePlots = FALSE,
                    method = "countxo", re.est.map = TRUE,
                    map.function = "kosambi", error.prob = 0.001,
                    verbose = TRUE,
                    ...){
  if(is.null(chr)) chr <- chrnames(cross)
  for(i in chr){
    if(verbose) cat("running ripple for chromosome: ", i,"\n")
    diff<-1
    reps<-1

    while(diff>0){
      if(verbose & repeatloop) cat("ripple run #", reps, "... ")
      if(makePlots) plot.rf(cross,
                            chr = i,
                            main = paste("chr",i,"before ripple"))
      rip<-summary(ripple(cross,
                          chr = i,
                          window = window,
                          method = method,
                          verbose = F,
                          ...))
      diff<-rip[1,ncol(rip)] - rip[2,ncol(rip)]
      if(diff>0){
        cross <- switch.order(cross, i, rip[2,-ncol(rip)])
        if(verbose) cat("n crossovers reduced by", diff,"\n")
        if(makePlots) plot.rf(cross,
                              chr = i,
                              main = paste("chr",i,"after ripple"))
      }else{
        if(verbose) cat("n crossovers not reduced\n")
      }
      reps<-reps+1
      if(!repeatloop) diff <- 0
    }
  }
  if(re.est.map){
    if(verbose) cat("re-estimating genetic map\n")
    em<-est.map(cross, map.function = map.function, error.prob = error.prob)
    cross<-replace.map(cross, map = em)
  }
  return(cross)
}
