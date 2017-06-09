#' @title Find the cM location of genes.
#'
#' @description
#' \code{findGenecM} Using the physical position of genetic markers,
#' infer the mapping position of every gene.
#'
#' @param cross The qtl cross object.
#' @param gff The .gff file containing information about each gene. This object must be of
#' standard format containing field described in details.
#' @param gffCols If the gff file does not follow the standard format, this vector
#' specifies the chr,feature, start, end, strand and attribute columns of the supplied gff-like file.
#' @param attributeParse Character vector of strings to drop from the first element of
#' the attribute column. Defaults to "ID=".
#' @param seqnameParse Character vector of strings in the seqname gff column to remove
#' to make the cross chromosome names match the gff seqname. Defaults to c("Chr","scaffold")
#' @param marker.info The qtlTools standard data.frame containing map and physiical position of
#' markers. See details. The base-pair positions of the markers must be known.
#' @param dropNonColinearMarkers logical, should markers that are not in the right bp order
#' be dropped?
#' @param verbose Logical, should updates be reported?
#' @param ... Not currently in use.
#'
#' @details
#' standard gff fields are as follows:
#' \enumerate{
#'    \item{seqname: name of the chromosome or scaffold}
#'    \item{source: name of the program that generated this feature,
#'           or the data source (database or project name)}
#'    \item{feature: feature type name, e.g. Gene, Variation, Similarity
#'          **Note** the term "Gene" must be present in this column}
#'    \item{start: Start position of the feature,
#'          with sequence numbering starting at 1.}
#'    \item{end: End position of the feature,
#'          with sequence numbering starting at 1.}
#'    \item{score: A floating point value.}
#'    \item{strand: defined as + (forward) or - (reverse)}
#'    \item{frame: One of '0', '1' or '2'.
#'          '0' indicates that the first base of the feature is the first base of a codon,
#'          '1' that the second base is the first base of a codon, and so on..}
#'    \item{attribute: A semicolon-separated list of tag-value pairs,
#'          providing additional information about each feature.}
#' }
#' marker.info fields - names must match exactly. The first three fields can be
#' generated using qtlTools::pullMap(cross)
#' \enumerate{
#'      \item{marker.name: Marker names (rownames from pull.map with as.table=T)}
#'      \item{chr: the chromosome of the marker}
#'      \item{pos: the mapping position of the marker}
#'      \item{bp: the base-pair position of the marker}
#' }
#'
#'
#' @return The gff file, with three added columns:
#' \itemize{
#'   \item{"geneID"}{the parsed name of the attribute}
#'   \item{"bp"}{the average base pair position}
#'   \item{"pos"}{the inferred cM position}
#' }
#'
#' @examples
#' \dontrun{
#' library(qtl)
#' library(qtlTools)
#' data(multitrait)
#' map<-pullMap(multitrait)
#' #simulate the bp positions of the markers with
#' #low recombination at the center of the chromosome
#' map$bp<-0
#' for(i in unique(map$chr)){
#'   n<-sum(map$chr==i)
#'   p<-sin((1:n/n)*pi)
#'   map$bp[map$chr==i]<-cumsum(p*1000000)
#' }
#'
#'
#' #simulate a gff w/ 1000 genes
#' gff<-data.frame(chr = rep(paste0("scaffold_",1:5),each = 200),
#'    feature = rep("gene",1000),
#'    start = rep(seq(from = 0, to = max(map$bp), length = 200), 5),
#'    end = rep(seq(from = 0, to = max(map$bp), length = 200))+1000,
#'    strand = rep("+",1000),
#'    attribute = paste0("gene",1:1000,";","gene",1:1000,".1"), stringsAsFactors=F)
#' test<-findGenecM(cross = multitrait, marker.info = map, gff = gff,
#'    gffCols = c("chr","feature","start","end","strand","attribute"))
#'
#' par(mfrow=c(3,2))
#' for(i in unique(map$chr)){
#' with(test[test$chr==i,], plot(pos,bp, col="grey"))
#' with(map[map$chr==i,], points(pos,bp, col=i, pch = 19, cex=.8))
#' }
#' }
#' @import qtl
#' @export

findGenecM<-function(cross, marker.info, gff, gffCols = NULL,
                     attributeParse = c("ID="),seqnameParse = c("Chr","scaffold_"),
                     dropNonColinearMarkers=TRUE, verbose = TRUE,...){

  dropNonColMar<-function(map, plotit = F, maxBpPercM = 5000000){
    tdiff<-function(y){
      d1<-diff(c(y,y[length(y)]))
      d2<-diff(c(0,y))
      return(d1<0|d2<0)
    }
    tmp<-map
    tmp$ord<-1:nrow(tmp)
    good<-unlist(lapply(unique(map$chr), function(x){
      tc<-tmp[tmp$chr == x,]
      tc<-tc[order(tc$pos),]
      tc$bpPcm<-
      d<-tdiff(tc$bp)
      bads<-numeric()
      while(any(d)){
        bads<-c(bads, tc$ord[d])
        tc<-tc[!d,]
        d<-tdiff(tc$bp)
      }
      return(tc$ord)
    }))
    return(map[tmp$ord %in% good,])
  }
  if(dropNonColinearMarkers){
    marker.info<-dropNonColMar(marker.info)
  }

  if(is.null(gffCols) & ncol(gff) != 9)
    stop("gff file must be of standard format, see details")
  if(!all(c("marker.name","chr","pos","bp") %in% colnames(marker.info)))
    stop("marker.info must be of standard format, see details")

  if(!is.null(gffCols)){
    if(length(gffCols) != 6)
      stop("if supplied, gffCols must be a vector length 6")
    if(!all(gffCols %in% colnames(gff)))
      stop("gffCols must be a vector that matches column names in the gff file")
    gff<-data.frame(chr = gff[,gffCols[1]],
                    source = NA,
                    feature = gff[,gffCols[2]],
                    start = gff[,gffCols[3]],
                    end = gff[,gffCols[4]],
                    score = NA,
                    strand = gff[,gffCols[5]],
                    frame = NA,
                    attribute = gff[,gffCols[6]],
                    stringsAsFactors=F)
  }

  if(verbose) cat("parsing attributes column of gff file\n")
  gff$geneID<-sapply(gff[,9], function(x) strsplit(x,";")[[1]][1])
  for(i in attributeParse){
    gff$geneID<-gsub(i,"",gff$geneID, fixed=T)
  }

  if(verbose) cat("culling chromosomes to those in the cross\n")
  colnames(gff)<-c("chr","source","feature","start","end",
                   "score","strand","frame","attribute","geneID")
  for(i in seqnameParse){
    gff$chr<-gsub(i,"",gff$chr, fixed=T)
  }
  
  gff$chr<-as.character(gff$chr)
  gff<-gff[gff$chr %in% as.character(chrnames(cross)),]

  gff$bp<-(gff[,4]+gff[,5])/2
  if(verbose) cat("inferring mapping position for:\n")
  out<-lapply(unique(gff$chr), function(i){
    tgff<-gff[gff$chr==i,]
    if(verbose) cat("chr ",i, " (n. features = ", nrow(tgff),")\n", sep="")
    tmap<-marker.info[marker.info$chr == i,]

    outint<-lapply(2:nrow(tmap), function(j) {
      bpcm<-tmap[(j-1):j,c("pos","bp")]
      gffbp<-tgff[tgff$bp >= min(bpcm$bp) &
                    tgff$bp < max(bpcm$bp),]
      if(nrow(gffbp)>=1){
        mod<-lm(pos~bp, data = bpcm)
        gffbp$pos = predict(mod, newdata = gffbp)
      }
      return(gffbp)
    })
    outchr <- do.call(rbind, outint)
    if(any(tgff$bp < min(tmap$bp))){
      outlow <- tgff[tgff$bp < min(tmap$bp), ]
      outlow$pos <- 0
      outchr<-rbind(outlow, outchr)
    }
    if(any(tgff$bp > max(tmap$bp))){
      outhi <- tgff[tgff$bp >= max(tmap$bp), ]
      outhi$pos <- max(tmap$pos)
      outchr<-rbind(outchr, outhi)
    }
    return(outchr)
    })

  return(do.call(rbind,out))
}
