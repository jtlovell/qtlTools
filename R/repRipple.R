#' @title Method to improve a genetic map.
#'
#' @description
#' \code{repRipple} Iteratively run ripple.
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
#' @param ... Additional arguments passed to qtl::ripple
#' @return A new cross object with the optimal rippled marker order
#'
#' @import qtl
#' @export
repRipple<-function(cross, chr = NULL, window=6, repeatloop = TRUE, makePlots = FALSE, method = "countxo",
                       re.est.map = TRUE, map.function = "kosambi", error.prob = 0.001, verbose = T,
                       ...){
  if(is.null(chr)) chr <- chrnames(cross)
  for(i in chr){
    if(verbose) cat("running ripple for chromosome: ", i,"\n")
    diff<-1
    reps<-1
    if(repeatloop){
      while(diff>0){
        if(verbose) cat("ripple run #", reps, "... ")
        reps<-reps+1
        if(makePlots) plot.rf(cross,
                              chr = i,
                              main = paste("chr",i,"before ripple"))

        rip<-ripple(cross,
                    chr = i,
                    window = window,
                    method = method,
                    verbose = F,
                    ...)
        index<-nmar(cross)[which(chrnames(cross)==i)]+1
        diff<-rip[1,index] - rip[2,index]
        if(diff>0){
          cross <- switch.order(cross, i, rip[2,])
          if(verbose) cat("n crossovers reduced by", diff,"\n")
        }else{
          if(verbose) cat("n crossovers not reduced\n")
        }

        if(makePlots) plot.rf(cross,
                              chr = i,
                              main = paste("chr",i,"after ripple"))
      }
    }else{
      if(makePlots) plot.rf(cross,
                            chr = i,
                            main = paste("chr",i,"before ripple"))

      rip<-ripple(cross,
                  chr = i,
                  window = window,
                  method = method,
                  verbose = F,
                  ...)

      cross <- switch.order(cross, i, rip[2,])

      if(makePlots) plot.rf(cross,
                            chr = i,
                            main = paste("chr",i,"after ripple"))
    }
  }
  if(re.est.map){
    if(verbose) cat("re-estimating genetic map\n")
    em<-est.map(cross, map.function = map.function, error.prob = error.prob)
    cross<-replace.map(cross, map = em)
  }
  return(cross)
}
