#' @title Find the location of a cis eQTL
#'
#' @description
#' \code{find.ciseqtl} Find the mapping position of the highest lod score
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
#' @examples
#' library(qtl)
#' data(fake.f2)
#' write.cross(fake.f2, format = "csv", filestem = "test") ### Must be "csv" or "csvs" formats.
#' old.x <- markernames(fake.f2)[grep("X", markernames(fake.f2))]
#' new.x <- c("X1", "X2", "X3")
#' renameMarkers(cross = fake.f2, crossfile = "test.csv",
#'    oldnames = old.x, newnames = new.x, outputName = "new.fake.f2.csv")
#' new.fake.f2 <- read.cross(file = "new.fake.f2.csv", format = "csv")
#' cbind(markernames(fake.f2), markernames(new.fake.f2))
#' @import qtl
#' @export



scan4trans<-function(cross,
                     phe,
                     perms = NULL,
                     nperms = 1000
                     k = 0,
                     cisQTL = NULL,
                     cisGenotypeCovariate = NULL,
                     covar = NULL,
                     cisWindow = 40,
                     verbose = FALSE,
                     ...){
  if(verbose) cat("scanning for trans eQTL\n")
  gp<-data.frame(m[[1]])
  gps<-data.frame(covar=covar$covar, gp[,1:ncol(gp)-1])
  s1.trans<-scanone(cross, model="normal",
                    method="hk", addcovar=gps, intcovar=covar, pheno.col=phe)

  if(is.null(perms)){
    if(verbose) cat("running permutations\n")
    perms <- GWERk(cross.in=cross2, phe="Unique.ID", nperms=1000, addcovar=gps, intcovar=covar, k=k)
  }

  if(!is.null(cisQTL)){
    cisChr=cisQTL$chr
    cisPos=cisQTL$pos
    cisWind<-c(cisPos-cisWindow, cisPos+cisWindow)
    s1.trans$lod[s1.trans$chr == cisChr & s1.trans$pos>cisWind[1] & s1.trans$pos<cisWind[2]]<-0
  }
  summary(s1.trans)


  gp<-data.frame(cis.scans[[phe]][[1]])
  gps<-data.frame(covar=covar$covar, gp[,1:ncol(gp)-1])

  s1<-scanone(cross, model="normal",
                method="hk", addcovar=gps, intcovar=covar, pheno.col=phe)
  perms<-scanone(cross, model="normal",
              method="hk", addcovar=gps, intcovar=gps, pheno.col=phe, n.perm=100)
  perms<-GWERk(cross.in = cross, phe = phe, k = 0, nperms = 1000,
                           printUpdate = TRUE, sameSeed = FALSE,
                           method="hk", addcovar = covar, intcovar = covar)
  colnames(perms)<-"lod"
  class(perms)<-"scanoneperm"
  summary(s1, perms=perms, pvalues=TRUE)
  plot(s1)
  add


  ind.ciseqtl(cross=cross, phe=x, pos.gene=annot$predictedCM, chr.gene=annot$lg,
              method="hk", addcovar=covar, intcovar=covar, outputType="qtl")

  covar2<-cbind()
  if(is.null(perms)){
    sig1<-data.frame(summary(s1, ci.function = "bayesint", format = "tabByCol")
    sig1<-sig1[sig1$lod >=qtlThreshold,]
    sig1$p.value = NA
  }else{
    sig1<-data.frame(summary(s1, perms=perms, pvalues=TRUE))
    sig1<-sig1[sig1$lod >=qtlThreshold,]
  }
  summary(s1)
  s1<-scanone(cross,
              pheno.col = phe,
              chr = chr.gene, ...)
  best.pos <- s1[s1$pos <= pos.gene + cmWindow &
                   s1$pos >= pos.gene - cmWindow, ]
  best.lod <- max(best.pos$lod)
  best.pos <- best.pos$pos[best.pos$lod == max(best.pos$lod)]
  if(outputType=="data.frame"){
    has.cis.qtl <- ifelse(best.lod > lodThreshold, "yes", "no")
    out<-data.frame(phe, chr.gene, pos.gene, best.pos, best.lod, has.cis.qtl)
    colnames(out)<-c("phenotype","chr","pos.gene","pos.cis.qtl","lod","has.cis.qtl")
    return(out)
  }else{
    if(!is.null(cross$geno[[1]]$draws)){
      m<-makeqtl(cross, chr=chr.gene, pos=best.pos, what = "draws")
    }else{
      if(!is.null(cross$geno[[1]]$prob)){
        m<-makeqtl(cross, chr=chr.gene, pos=best.pos, what = "prob")
      }else{
        cross<-calc.genoprob(cross)
        m<-makeqtl(cross, chr=chr.gene, pos=best.pos, what = "prob")
      }
    }
    return(m)
  }
}
