#' @title Quantile normalize vector of numeric data
#'
#' @description
#' \code{quantnorm} Rank-normalize (quantile normalize) data to enforce strict normality
#' of data. Useful if assuming identical distributions across many response variables
#'
#' @param x The numeric vector to quantile normalize
#' @param drop0s Logical, should 0's be transformed to NAs
#'
#' @return a vector of quantile normalized data
#'
#' @examples
#' library(qtlTools)
#' data(fake.bc)
#' cross<-fake.bc
#' par(mfrow = c(2,1))
#' hist(pull.pheno(cross, pheno.col = "pheno2"), breaks=20,
#'    main = "before quantile normalization",
#'    xlab = "raw pheno1")
#' cross<-transformPheno(cross, pheno.col = "pheno1", quantnorm)
#' hist(pull.pheno(cross, pheno.col = "pheno1"), breaks=20,
#'    main = "after quantile normalization",
#'    xlab = "qn pheno1")
#'
#' @export
quantnorm<-function(x, drop0s=FALSE) {
  x<-as.numeric(x)
  if(drop0s){
    x[x==0]<-NA
  }
  n=sum(!is.na(x),na.rm=T)
  x=rank(x)/(n+1)
  x=qnorm(x)
  x[is.infinite(x)]=NA
  x[x=="NaN"]<-NA
  x
}
