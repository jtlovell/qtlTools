#' @title Quantile normalize vector of numeric data
#'
#' @description
#' \code{quantnorm} Rank-normalize (quantile normalize) data to enforce strict normality
#' of data. Useful if assuming identical distributions across many response variables
#'
#' @param x The numeric vector to quantile normalize
#'
#' @return a vector of quantile normalized data
#' @examples
#' par(mfrow = c(2,1))
#' hist(x <- sample(1:100, 50), breaks = 12)
#' hist(quantnorm(x), breaks = 12)
#'
#' library(qtl)
#' data(fake.f2)
#' hist(pull.pheno(fake.f2, "phenotype"))
#' hist(quantnorm(pull.pheno(fake.f2, "phenotype")))
#' fake.f2<-transformPheno(fake.f2, pheno.col = "phenotype", transf = quantnorm)
#' hist(pull.pheno(fake.f2, "phenotype"))
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
