#' @title Find the location of a cis eQTL
#'
#' @description
#' \code{findCisEQTL} find the mapping position of the highest QTL peak
#' within a window surrounding the known position of a gene.
#'
#' @param cross The qtl cross object.
#' @param pullSigQTL.output The results of pullSigQTL with returnQTLModel = FALSE.
#' @param findGenePos.output The findGenePos output, containing the mapping position of all
#' genes.
#' @param gff The .gff file containing information about each gene. This object must be of
#' standard format containing field described in details.
#' @param marker.info The qtlTools standard data.frame containing map and physiical position of
#' markers. See details. The base-pair positions of the markers must be known.
#' @param cmWindow The distance from pos.gene that the cis eQTL may reside. Default is 10, but
#' if the confidence of the mapping position is high, a smaller value is suggested.
#' @param useConfInt Logical, should the confidence iterval be used to determine cis QTL? If so,
#' cmWindow is ignored
#' @param outputType Whether to return a data.frame or qtl object (from makeqtl)
#' @param ... Additional arguments passed to scanone.
#'
#' @details
#' gff fields - the names of the fields can differ, but column numbers must follow:
#' 1. seqname: name of the chromosome or scaffold
#' 2. source: name of the program that generated this feature,
#'           or the data source (database or project name)
#' 3. feature: feature type name, e.g. Gene, Variation, Similarity
#'          **Note** the term "Gene" must be present in this column
#' 4. start: Start position of the feature, with sequence numbering starting at 1.
#' 5. end: End position of the feature, with sequence numbering starting at 1.
#' 6. score: A floating point value.
#' 7. strand: defined as + (forward) or - (reverse).
#' 8. frame: One of '0', '1' or '2'.
#'          '0' indicates that the first base of the feature is the first base of a codon,
#'          '1' that the second base is the first base of a codon, and so on..
#' 9. attribute: A semicolon-separated list of tag-value pairs,
#'          providing additional information about each feature.
#'
#' marker.info fields - names must match exactly. The first three fields can be
#' generated using qtlTools::pullMap(cross)
#' 1. marker.name: Marker names (rownames from pull.map with as.table=T)
#' 2. chr: the chromosome of the marker
#' 3. pos: the mapping position of the marker
#' 4. bp: the base-pair position of the marker
#'
#'
#' @return Either a qtl object or a dataframe, indicating the position of a cis eQTL.
#'
#' @import qtl
#' @export
find.ciseqtl<-function(cross, pheno.col, pos.gene, chr.gene, cmWindow = 10,
                       lodThreshold = 3, outputType = "data.frame", ...){
  s1<-scanone(cross,
              pheno.col = pheno.col,
              chr = chr.gene, ...)
  best.pos <- s1[s1$pos <= pos.gene + cmWindow &
                   s1$pos >= pos.gene - cmWindow, ]
  best.lod <- max(best.pos$lod)
  best.pos <- best.pos$pos[best.pos$lod == max(best.pos$lod)]

  has.cis.qtl <- ifelse(best.lod > lodThreshold, "yes", "no")
  out<-data.frame(pheno.col, chr.gene, pos.gene, best.pos, best.lod, has.cis.qtl)
  colnames(out)<-c("phenotype","chr","pos.gene","pos.cis.qtl","lod","has.cis.qtl")
  if(!is.null(cross$geno[[1]]$draws)){
    m<-makeqtl(cross, chr=chr.gene, pos=best.pos, what = "draws")
  }else{
    if(!is.null(cross$geno[[1]]$prob)){
      m<-makeqtl(cross, chr=chr.gene, pos=best.pos, what = "prob")
    }else{
      cross<-calc.genoprob(cross)
      m<-makeqtl(cross, chr=chr.gene, pos=best.pos, what = "prob")
    }
    return(list(cisQTL=m, cisDF=out, s1=s1))
  }
}
