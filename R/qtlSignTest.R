#' @title Conduct the sign test for expression bias
#'
#' @description
#' \code{qtlSignTest} Uses the Fraser etc. sign test to assess neutral evolution
#' for a set of expression phenotypes
#'
#' @param effect Signed cis-eQTL or allele-specific expression log2-fold change data
#' @param trans.effect If effect is signed cis-eQTL effect, one can test the independence
#' of cis and trans eQTL evolution. Here provide effect effect of trans QTL.
#' @param id label for output
#' @param verbose should results be printed to console?
#' @param ... Additional arguments to pass on to fisher.test / binom.test, such as
#' null hypotheses.
#'
#' @details See ...
#' @return The results of a binomial test (or Fisher's exact test if trans.dir is specified).
#' @export

qtlSignTest<-function(effect = NULL, trans.effect = NULL, verbose=T, ...){
  if(is.null(trans.effect)){
    if(verbose)
      cat("running sign test as: bias of single effect / gene\n")
    effect<-effect[!is.na(effect) & is.finite(effect) & !is.nan(effect)]
    nup = sum(effect>0)
    ntot = sum(effect!=0)
    out<-binom.test(nup,ntot, p = 0.5, alternative = "two.sided", ...)
    if(verbose)
      cat("n. up regulated = ", nup, " out of ", ntot,
          " genes (",round(nup/ntot,3),"%)\npvalue = ",out$p.value, sep = "")
    ret<-with(out, data.frame(n.up = nup, ntotal = ntot,
                              prob.success = estimate,
                              ci.low = conf.int[1], ci.high = conf.int[2],
                              p.value = p.value))
  }else{
    if(!is.null(trans.effect)){
      if(verbose)
        cat("running sign test as independence of cis and trans effect directions\n")
      effect<-effect[!is.na(effect) & is.finite(effect) & !is.nan(effect)]
      trans.effect<-trans.effect[!is.na(trans.effect) & is.finite(trans.effect) & !is.nan(trans.effect)]
      upup = sum(effect>0 & trans.effect>0)
      updown = sum(effect>0 & trans.effect<0)
      downup = sum(effect<0 & trans.effect>0)
      downdown = sum(effect<0 & trans.effect<0)
      ctsame = sum(upup, downdown)
      ctdiff = sum(updown, downup)
      out<-binom.test(ctsame,sum(ctdiff,ctsame), p = 0.5, alternative = "two.sided", ...)
      if(verbose)
        cat("cis&trans same = ", ctsame,
            "; cis&trans different = ", ctdiff,
            "; ratio = ", round(out$estimate,2),
            "; p.value = ",out$p.value, sep = "")
      ret<-with(out, data.frame(ctsame = ctsame, ctdiff = ctdiff,
                                prob.success = estimate,
                                ci.low = conf.int[1], ci.high = conf.int[2],
                                p.value = p.value))
    }else{
      stop("a vector of effects is required\n")
    }
  }
  return(ret)
}
