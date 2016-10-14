#' @title Find the location of a cis eQTL
#'
#' @description
#' \code{scan4trans} Find the mapping position of the highest lod score
#' within a window surrounding the known position of a gene.
#'
#' @param cross The qtl cross object with marker names that need to be changed.
#' @param phe Character or numeric vector indicating the phenotype to be tested
#' @param chr.gene The name of the chromosome that contains the focal genes
#' @param pos.gene The suspected mapping position of the gene - may be inferred from bp2cm
#' @param cmWindow The distance from pos.gene that the cis eQTL may reside. Default is 10, but
#' if the confidence of the mapping position is high, a smaller value is suggested.
#' @param outputType Whether to return a data.frame or qtl object (from makeqtl)
#' @param lodThreshold If outputType = data.frame, the presence of a cis.eqtl is defined by whether
#' any marker within cmWindow has a lod score > lodThreshold.
#' @param lodThreshold

#'
#' @return Either a qtl object or a dataframe, indicating the position of a cis eQTL.
#'
#' @import qtl
#' @export

scan4trans<-function(cross,
                     phe,
                     perms = NULL,
                     nperms = 1000,
                     k = 0,
                     cisQTL = NULL,
                     cisGenotypeCovariate = TRUE,
                     covar = NULL,
                     cisWindow = 40,
                     verbose = FALSE,
                     cisname = "Q1",
                     refine = TRUE,
                     ...){
  if(verbose) cat("scanning for trans eQTL\n")

  if(cisGenotypeCovariate & !is.null(cisQTL)){
    gp<-apply(data.frame(cisQTL[[1]]),1,function(x) which(x==max(x)))
    gps<-data.frame(covar=covar$covar, gp)
  }else{
    gps<-covar
  }

  s1.trans<-scanone(cross, model="normal",
                    method="hk", addcovar=gps, intcovar=covar, pheno.col=phe)

  if(is.null(perms)){
    if(verbose) cat("running permutations\n")
    perms <- GWERk(cross.in=cross, phe=phe, nperms=100, addcovar=gps, intcovar=covar, k=k, method="hk")
  }
  if(class(perms)!="scanoneperm"){
    perms=matrix(perms)
    colnames(perms)[1]<-"lod"
    class(perms)<-"scanoneperm"
  }
  if(!is.null(cisQTL)){
    cisChr=cisQTL$chr
    cisPos=cisQTL$pos
    cisWind<-c(cisPos-cisWindow, cisPos+cisWindow)
    s1.trans$lod[s1.trans$chr == cisChr & s1.trans$pos>cisWind[1] & s1.trans$pos<cisWind[2]]<-0
  }else{
    cisChr=NULL
    cisPos=NULL
  }

  sumtrans<-data.frame(summary(s1.trans, perms=perms, pvalues=T))
  sigtrans<-sumtrans[sumtrans$pval<=0.05,]

  chr <- c(cisChr, sigtrans$chr)
  pos <- c(cisPos, sigtrans$pos)

  if(length(chr)==0){
    qtl=NULL; sdrop=NULL;s1.trans=NULL
  }else{
    o <- order(factor(chr, levels=names(cross$geno)), chr)
    qtl <- makeqtl(cross, chr[o], pos[o], what="prob")

    formula<-paste("y ~ ",
                   paste(qtl$altname, collapse=" + "),
                   " + ",
                   paste(paste(qtl$altname, "covar", sep = ":"), collapse = " + "),
                   " + covar",
                   sep = "")
    qtl$name<-qtl$altname

    attr(qtl, "formula") <- deparseQTLformula(formula)

    fall<-fitqtl(cross, pheno.col=phe, formula=formula(qtl), qtl=qtl, covar=covar, method="hk")
    sall<-summary(fall)
    sdrop<-data.frame(sall$result.drop)
  }
  return(list(model=qtl, dropstats=sdrop, s1 = s1.trans))
}
