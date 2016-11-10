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
    out<-binom.test(nup,ntot, ...)
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
      nup = sum(effect>0)
      ndown = sum(effect<0)
      nup.trans = sum(trans.effect>0)
      ndown.trans = sum(trans.effect<0)
      out<-fisher.test(rbind(c(nup,ndown), c(nup.trans,ndown.trans)), ...)
      if(verbose)
        cat("cis = ", round(nup/ndown,2),
            "; trans = ", round(nup.trans/ndown.trans,2),
            "; odds ratio = ", round(out$estimate,2),
            "; p.value = ",out$p.value, sep = "")
      ret<-with(out, data.frame(n.up_cis = nup, n.down_cis = ndown,
                                n.up_trans = nup.trans, n.down_trans = ndown.trans,
                                odds.ratio = estimate,
                                ci.low = conf.int[1], ci.high = conf.int[2],
                                p.value = p.value))
    }else{
      stop("a vector of effects is required\n")
    }
  }
  return(ret)
}
