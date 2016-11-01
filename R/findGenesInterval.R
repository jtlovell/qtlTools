#' @title Find candidate genes for QTL.
#'
#' @description
#' \code{findGenesInterval} Use the mapping position (inferred) of each gene to
#' find potential candidates under a QTL confidence interval
#'
#' @param findGenecM.output Output from findGenecM
#' If supplied, geneID, geneChr, geneBp and genecM are ignored
#' @param calcCis.output Output from calcCis.
#' If supplied, chr, lowposition and highposition are ignored
#' @param qtlname Identifier of each qtl
#' @param chr QTL Chromosome
#' @param lowposition Lower confidence interval bound (cM)
#' @param highposition Upper confidence interval bound (cM)
#' @param geneID The name of each gene
#' @param geneChr The chromosome of each gene
#' @param geneBp The basepair position of each gene (not currently in use)
#' @param genecM The mapping position of each gene
#' @param ... Not currently in use.
#'
#' @return A named list where each element corresponds to a provided QTL
#' confidence interval. Character vectors in each element display potential
#' candidate genes for each QTL.
#'
#' @export

findGenesInterval<-function(findGenecM.output = NULL, calcCis.output = NULL,
                            qtlname = NULL, chr, lowposition, highposition,
                            geneID, geneChr, geneBp, genecM){
  if(is.null(calcCis.output)){
    if(any(is.null(c(chr,l,h))))
      stop("if calcCis.output is not specified, chr, l and h must be\n")
    if(is.null(qtlname)){
      qtlname<-1:length(chr)
    }
    ci<-data.frame(qtlname = qtlname,
                   chr = chr,
                   lowposition = lowposition,
                   highposition = highposition,
                   stringsAsFactors = F)
  }else{
    if(is.null(qtlname)){
      calcCis.output$qtlname<-with(calcCis.output, paste(chr, round(pos,0), sep = "@"))
    }else{
      calcCis.output$qtlname<-qtlname
    }
    ci<-calcCis.output[,c("qtlname","chr","lowposition","highposition")]
  }

  if(is.null(findGenecM.output)){
    if(any(is.null(c(geneChr,genecM))))
      stop("if findGenecM.output is not specified, geneChr and genecM must be\n")

    gff<-data.frame(geneID = geneID,
                    chr = geneChr,
                    bp = geneBp,
                    pos = genecM,
                   stringsAsFactors = F)
  }else{
    gff<-findGenecM.output[,c("geneID","chr","bp","pos")]
  }
  out<-lapply(1:nrow(ci), function(x){
    tqtl<-ci$qtlname[x]
    tchr<-ci$chr[x]
    tl<-ci$lowposition[x]
    th<-ci$highposition[x]
    tgff<-gff[with(gff,
                   chr == tchr &
                     pos>=tl &
                     pos<=th),]
    return(tgff$geneID)
  })
  if(length(out)==1){
    return(out[[1]])
  }else{
    names(out)<-ci$qtlname
  }
  return(out)
}
