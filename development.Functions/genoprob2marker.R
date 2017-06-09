#' @title Get allelic/genotype probabilities from QTL models
#'
#' @description
#' \code{genoprob2marker} Not intended for external use.
#'
#' @param cross The qtl cross object with marker names that need to be changed.
#' @param useGrid Should the genotype probability matrix be reduced to a grid?
#' @param chr Should the analysis be limited to specific chromosomes.
#' @param qtl If extracting data from a single qtl model, supply here.
#' @param prob.thresh What should the probability of an allele need to be to be called?
#' @param returnNumeric Should genotype calls or probabilities be returned?
#' @param genoNum The numeric encoding for each genotype.
#' @param ... additional arguments passed on to rrBLUP::mixed.solve
#'
#' @import qtl
#' @export

genoprob2marker<-function(cross, useGrid = TRUE, prob.thresh = 0, chr = NULL, qtl = NULL,
                          returnNumeric=TRUE, genoNum=NULL){
  to.gp<-function(m, prob.thresh = 0, genoNum = NULL, returnNumeric=TRUE, is.qtl = TRUE){
    if(is.qtl) m<-m[[1]]
    gp<-lapply(m, function(x) apply(x,1, function(y) {
      if(max(y) < prob.thresh){
        return(NA)
      }else{
        if(length(colnames(x)[which(y==max(y))])>1){
          return(NA)
        }else{
          if(returnNumeric){
            if(is.null(genoNum)){
              genoNum<-seq(from = 0, to = 1, length.out = ncol(x))
            }
            out<-sum(genoNum*y)
            return(out)
          }else{
            return(colnames(x)[which(y==max(y))])
          }
        }
      }
    }))
    gp<-data.frame(do.call(cbind,gp))
    if(!returnNumeric){
      for(i in colnames(gp)) gp[,i]<-as.character(gp[,i])
    }
    return(gp)
  }

  if ("prob" %in% names(cross$geno[[1]]) & useGrid) {
    stepwidth <- attr(cross$geno[[1]]$prob, "stepwidth")
    if (stepwidth != "fixed") {
      stop("You need to have run calc.genoprob with stepwidth=\"fixed\".")
    }
  }
  if(!is.null(qtl)){
    returnNumeric<-FALSE
    gp<-to.gp(m = qtl,
              prob.thresh = prob.thresh,
              genoNum = genoNum,
              returnNumeric=returnNumeric)
    colnames(gp)<-qtl$altname
  }else{
    if(is.null(chr)){
      chr = chrnames(cross)[chrnames(cross) != "X"]
    }
    if(class(cross)[1] == "4way" & useGrid){
      s1<-scanone(cross, method = "hk", chr=chr)
      s1<-s1[grep(".loc",rownames(s1), fixed=T),]
    }else{
      if(useGrid){
        grid<-reduce2grid(cross)
        s1<-scanone(grid, method = "hk", chr=chr)
      }else{
        s1<-scanone(cross, method = "hk", chr=chr)
      }
    }
    m<-makeqtl(cross, chr = s1$chr,
               pos = s1$pos,
               what="prob")
    if(class(cross)[1] != "4way"){
      returnNumeric<-TRUE
      gp<-to.gp(m = m,
                prob.thresh = prob.thresh,
                genoNum = genoNum,
                returnNumeric=returnNumeric)
    }else{
      m1<-lapply(m[[1]], function(x) x[,1:2])
      m2<-lapply(m[[1]], function(x) x[,3:4])
      gp1<-to.gp(m = m1, is.qtl = FALSE,
                 prob.thresh = prob.thresh,
                 genoNum = genoNum,
                 returnNumeric=returnNumeric)
      gp2<-to.gp(m = m2, is.qtl = FALSE,
                 prob.thresh = prob.thresh,
                 genoNum = genoNum,
                 returnNumeric=returnNumeric)
      gp<-cbind(gp1, gp2)
    }
  }
  return(gp)
}
