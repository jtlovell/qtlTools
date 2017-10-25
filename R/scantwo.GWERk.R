#' @title Scanone qtl analysis permitting GWERk analysis.
#'
#' @description
#' \code{scantwo.GWERk}  see qtl::scanone for details. GWERk = 0 runs normal scanone
#' permutations. Otherwise GWERk perms are used.
#' @import qtl
#' @export
scantwo.GWERk <- function(cross, n.perm = 100, ..., gwerk = 1, use.scantwopermhk = T){
  if(use.scantwopermhk){
    pl<-lapply(chrnames(cross), function(x){
      set.seed(42)
      return(scantwopermhk(cross, chr = x, n.perm = n.perm, ...))
    })
    fulls<-sapply(pl, function(x) x[[1]])
    wh.2keep <- apply(fulls, 1, function(x) order(-x)[gwerk+1])
    pl.out<-lapply(1:length(wh.2keep), function(x) {
      tmp<-pl[[wh.2keep[x]]]
      return(lapply(tmp, function(y) y[x]))
    })
    out<-lapply(1:length(pl[[1]]), function(x){
      tmp<-sapply(pl.out, function(y) {
        unlist(y[[x]])
      })
      tmp<-data.frame(tmp)
      colnames(tmp)<-colnames(pl[[1]][[1]])
      return(as.matrix(tmp))
    })
    names(out)<-names(pl[[1]])
    class(out)<-"scantwoperm"
    return(out)
  }
}
