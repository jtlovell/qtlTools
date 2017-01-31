#' @title Format a genotype matrix to a R/qtl readable file
#'
#' @description
#' \code{geno2cross} Combine marker information and genotype data into
#' a R/qtl "csv" file.
#'
#' @param genomat A matrix of individual (row) * marker (col) genotype calls.
#' @param chr A vector of chromosome IDs. If null, taken from the column names
#' @param pos A vector of marker positions. Must be able to be coerced to
#' numeric (does not contain letters / special characters).
#' If null, taken from the column names .
#' @param id Individual identifiers. If NULL, taken from row names of genomat.
#' @param chrpos.sep If chr and pos are NULL, use this character to parse
#' marker names into chromosome and position information. For example, if
#' chrpos.sep = "_", a marker named 01_100000 would be parsed into chr = 01
#' pos = 100000.
#' @param crossfile The name of the cross file to be written. Defaults to store
#' in current working directory. Specify full path if you want to write to another
#' location.
#'
#' @details More here
#'
#' @return  Nothing - writes a file called cross.file
#'
#' @examples
#' \dontrun{
#' ... more here ...
#' }
#' @import qtl
#' @export
geno2cross<-function(genomat,
                     chr=NULL, pos=NULL, id=NULL,
                     chrpos.sep = "_",
                     crossfile = "cross.csv"){
  if(is.null(chr) & is.null(pos)){
    if(!all(grepl(chrpos.sep,colnames(genomat), fixed=T)))
      stop("if chr and pos are not specified, they must be provided in the column names
           of genomat, separated by chrpos.sep character\n")
    chr<-sapply(colnames(genomat), function(x) strsplit(x, chrpos.sep)[[1]][1])
    pos<-sapply(colnames(genomat), function(x) strsplit(x, chrpos.sep)[[1]][2])
  }

  if(is.null(id)){
    if(is.null(rownames(genomat)))
      stop("if id is not specified, genomat rownames must be supplied\n")
  }else{
    rownames(genomat)<-id
  }

  towrite<-genomat[,order(chr, as.numeric(pos))]
  chr<-chr[order(chr, as.numeric(pos))]
  pos<-pos[order(chr, as.numeric(pos))]
  info<-t(cbind(chr,pos))
  towrite2<-rbind(info,towrite)

  id[1:2]<-""
  towrite2<-cbind(id, towrite2)
  write.csv(towrite2,
            file = crossfile,
            quote = F, row.names=F)
  cat("cross file written to", crossfile,"\n")
}
