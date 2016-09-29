#' @title Method to improve a genetic map.
#'
#' @description
#' \code{repeatRipple} Iteratively run ripple.
#'
#' @param cross The qtl cross object to search
#' @param chr The chromosome to scan. Can be a vector of chromosome names or a single name.
#' @param window The window size passed on to ripple
#' @param method The method passed on to ripple
#' @param makePlots Should the recombination fraction be iteratively plotted?
#' @param repeatloop Should ripple be run once, or iteratively?
#' @param ... Additional arguments passed to qtl::ripple
#' @return A new cross object with the optimal rippled marker order
#'
#' @import qtl
#' @export
repeatRipple<-function(cross, chr = NULL, window=6, repeatloop = TRUE, makePlots = FALSE, method = "countxo", ...){
  if(is.null(chr)) chr <- chrnames(cross)
  for(i in chr){
    diff<-1
    if(repeatloop){
      while(diff>0){
        if(makePlots) plot.rf(cross,
                              chr = i,
                              main = paste("chr",i,"before ripple"))

        rip<-ripple(cross,
                    chr = i,
                    window = window,
                    method = method,
                    ...)

        cross <- switch.order(cross, i, rip[2,])
        index<-nmar(cross)[which(chrnames(cross)==i)]+1
        diff<-rip[1,index] - rip[2,index]

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
                  ...)

      cross <- switch.order(cross, i, rip[2,])

      if(makePlots) plot.rf(cross,
                            chr = i,
                            main = paste("chr",i,"after ripple"))
    }
  }
  return(cross)
}
