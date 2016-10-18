#' @title Method to summarize scanone results
#'
#' @description
#' \code{pullSigQTL} Uses qtlpvl to summarize the output of scanone,
#' then culls the output to only significant QTL peaks, based on permutations.
#'
#' @param cross The qtl cross.
#' @param pheno.col Character or numeric vector indicating the phenotype to be tested.
#' @param s1.output The output from scanone
#' @param perm.output The permutation output from scanone
#' @param chr The chromosome to be tested. Defaults to all chromosomes.
#' @param alpha The significance for permutations
#' @param returnQTLModel Logical, should a QTL model be returned (TRUE), or
#' should a culled output from qtlpvl::convert_scan1 be returned (FALSE)?
#' @param ... additional arguments passed on to summary.scanone,
#' such as controlAcrossCol.
#' @return Either QTL models or simplified and converted scanone summary.
#'
#' @examples
#' library(qtlTools)
#' data(fake.bc)
#' cross<-fake.bc
#' cross <- calc.genoprob(cross, step=2.5)
#' s1<-scanone(cross, method="hk", pheno.col=c("pheno1", "pheno2"))
#' perm<-scanone(cross, n.perm=100, method="hk",
#'    pheno.col=c("pheno1", "pheno2"), verbose=F)
#' pullSigQTL(cross, s1.output=s1, perm.output=perm)
#' pullSigQTL(cross, s1.output=s1, perm.output=perm, returnQTLModel=F)
#'
#' @export
pullSigQTL<-function(cross, s1.output, perm.output,
                     pheno.col=NULL, chr=NULL,
                     alpha = 0.05, returnQTLModel = TRUE,
                     ...){

  ########
  ########
  # convert_scan1, from qtlpvl...
  convert_scan1<- function(out, phenoname, chr = NULL)  {
    if(length(phenoname)>1){
      CHR <- out[, "chr"]
      POS <- out[, "pos"]
      out <- out[, phenoname]
      maxLOD <- matrix(NA, length(phenoname), length(chr))
      colnames(maxLOD) <- chr
      rownames(maxLOD) <- phenoname
      maxPOS <- maxLOD
      for (i in 1:length(chr)) {
        index <- CHR == chr[i]
        LOD <- out[index, ]
        maxLOD[, i] <- apply(LOD, 2, max)
        maxPOS[, i] <- POS[index][apply(LOD, 2, which.max)]
      }
      out1 <- data.frame(pheno = rep(phenoname, ncol(maxLOD)),
                         chr = rep(colnames(maxLOD), each = nrow(maxLOD)), lod1 = c(maxLOD),
                         stringsAsFactors = FALSE)
      out1$pos <- c(maxPOS)
    }else{
      out<-data.frame(summary(out))
      colnames(out)[3]<-"lod1"
      out$pheno<-phenoname

      out1 <- out[,c("pheno","chr","lod1","pos")]
    }

    return(out1)
  }
  ########
  ########

  if(is.null(chr)) chr <- chrnames(cross)
  if(is.null(pheno.col)) pheno.col<- names(s1.output)[-c(1:2)]
  maxs<-convert_scan1(s1.output, phenoname=pheno.col, chr = chr)
  sperms<-summary(perm.output, alpha = alpha, ...)
  perm.names<-attr(sperms,"dimnames")[[2]]
  sperms<-as.numeric(sperms)
  for(i in pheno.col){
    thresh<-sperms[perm.names==i]
    maxs<-maxs[-which(maxs$pheno == i & maxs$lod1 < thresh),]
  }
  if(returnQTLModel){
    mods<-lapply(pheno.col, function(i){
      temp<-maxs[maxs$pheno == i,]
      if(nrow(temp)==0) {
        out<-"NULL QTL Model"
      }else{
        if("prob" %in% names(cross$geno[[1]])){
          out<-makeqtl(cross, chr = temp$chr, pos = temp$pos, what = "prob")
        }else{
          out<-makeqtl(cross, chr = temp$chr, pos = temp$pos, what = "draws")
        }
      }
    })
    names(mods)<-pheno.col
    return(mods)
  }else{
    return(maxs)
  }
}
