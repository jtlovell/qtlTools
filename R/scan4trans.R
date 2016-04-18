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
                     ...){
  if(verbose) cat("scanning for trans eQTL\n")

  if(cisGenotypeCovariate){
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

  chrs <- c(cisChr, sigtrans$chr)
  poss <- c(cisPos, sigtrans$pos)
  qnames<-1:(1+nrow(sigtrans))
  qnames<-paste("Q", qnames, sep="")
  mall<-makeqtl(cross, chr=chrs, pos=poss, qtl.name=qnames, what = "prob")

  form.allint<-paste("y ~ ",
                     paste(qnames, collapse=" + "),
                     " + ",
                     paste(paste(qnames, "covar", sep = ":"), collapse = " + "),
                     " + covar",
                     sep = "")

  fall<-fitqtl(cross, pheno.col=phe, formula=form.allint, qtl=mall, covar=covar, method="hk")

  if(refine){
    mall<-refineqtl(cross, qtl = mall, formula = form.allint, pheno.col=phe, covar=covar, method = "hk",
                    verbose=F)
  }

  #summary(fitqtl(cross, pheno.col = phe, qtl = mall, formula = form.allint, method="hk", covar=covar))

  sall<-summary(fall)
  sdrop<-data.frame(sall$result.drop)
  nonsigQTL<-rownames(sdrop)[sdrop$Pvalue.F. == max(sdrop$Pvalue.F.) & sdrop$Pvalue.F. > 0.05 & rownames(sdrop)!="covar"][1]

  while(length(nonsigQTL)>0 & gsub(" ","", form.allint) != "y~covar"){
    print(nonsigQTL)

    if(nonsigQTL %in% mall$name &
       !paste(nonsigQTL, "covar", sep = ":") %in% nonsigQTL &
       paste(nonsigQTL, "covar", sep = ":") %in% rownames(sdrop)){
      nonsigQTL<-gsub(nonsigQTL,paste(i, "covar", sep = ":"), nonsigQTL)
    }

    form.allint<-gsub(paste(nonsigQTL,"+",sep=" "),"", form.allint, fixed=TRUE)

    if(gsub(" ","", form.allint) == "y~covar"){
      fall<-NULL
    }else{
      if(refine){
        small<-refineqtl(cross, qtl = mall, formula = form.allint,
                         pheno.col=phe, covar=covar, method = "hk",verbose=F)
      }
      fall<-fitqtl(cross, pheno.col=phe, formula=form.allint, qtl=mall, covar=covar, method="hk")
      sall<-summary(fall)
      sdrop<-data.frame(sall$result.drop)
      nonsigQTL<-rownames(sdrop)[sdrop$Pvalue.F. == max(sdrop$Pvalue.F.) & sdrop$Pvalue.F. > 0.05 & rownames(sdrop)!="covar"]
    }
  }
  if(verbose) cat("compiling stats")
  if(is.null(fall)){
    return(list(stats=NULL, model = NULL, formula = NULL))
  }else{
    stats<-qtlStats(cross, mod=mall, phe=phe, cisname = "Q1", covar = covar, form = form.allint)
    return(list(stats=stats, model = mall, formula = form.allint))
  }
}
