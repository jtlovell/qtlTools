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
find.ciseqtl<-function(cross, phe, pos.gene, chr.gene, cmWindow = 10,
                       lodThreshold = 3, outputType = "data.frame", ...){
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
