#' @title Find the effect of candidate genes based on gene expression.
#'
#' @description
#' \code{covScanQTL} Employs the covariate scan approach (Lovell et al. (2015),
#' Plant Cell), to rank the potential of candidate genes based on their effect
#' on QTL morphology. Can be run for a single phenotype (e.g. Lovell et al. (2015)),
#' or on a set of QTL underlying a trand-band. In the latter case, a boxplot of ranks
#' can be output.
#'
#' @param cross The qtl cross object with marker names that need to be changed.
#' @param pheno.col Character vector specifying the phenotypes modelled.
#' @param expression.covariates Numeric matrix with normalized expression values
#' of candidate genes. Rows must exactly match the individuals in the cross
#' @param addcovar data.frame of experimental additive covariates (see scanone)
#' @param intcovar data.frame of experimental interactive covariates (see scanone)
#' @param qtl A qtl model that supplies the position of the focal QTL
#' (see focalqtl.index). If qtl only contains a single locus, this is assumed
#' to be the qtl of interest. If a multiple qtl model is supplied, the genotypes
#' of those QTL are inferred using the Viterbi algorithm and added to the addcovar
#' dataframe of covariates.
#' @param focalqtl.index Numeric, which qtl in the model is the the QTL to test?
#' @param which.epiqtl Numeric, which qtl (if any) in the model should be forced
#' to interact with the focal qtl. The inferred genotypes of this marker are added to
#' the intcovar dataframe.
#' @param qtl.method The method passed to scanone
#' @param nperm If permutation tests are desired, specify as >0. These are not currently
#' computationally efficient and take forever when using eQTL data.
#' @param plotit Logical, when more than 1 pheno.y is specified, presents a boxplot of
#' covariate scan ranks.
#' @param verbose Logical, should updates be printed?
#'
#' @return A dataframe, containing the maximum scanone outputs at
#' chromosome chr and position pos for each phenotype and covariate.
#'
#' @examples
#' \dontrun{
#' data(multitrait)
#' cross<-subset(multitrait, ind = !is.na(pull.pheno(multitrait, 1)))
#' phe<-pull.pheno(cross, 1)
#' mult.fact<-exp(seq(from = 0, to = 50, length.out = 50))
#' #ssimulate some gene expression data, some of which are correlated with the phenotype
#' facs<-sapply(1:50, function(x){
#'   scale(sapply(scale(phe), function(y) rnorm(n = 1, mean = y, sd = mult.fact[x])))
#'   })
#' plot(sapply(1:50, function(x) cor(phe, facs[,x])),
#'      ylab = "cor. coef. (expression ~ chr5 QTL genotype)",
#'      xlab = "gene id")
#'
#' expression.covariates = facs
#' colnames(expression.covariates)<-paste0("gene",1:ncol(expression.covariates))
#' qtl = makeqtl(cross, chr = max(s1)$chr, pos = max(s1)$pos, what = "prob")
#' test<-covScanQTL(cross = cross,
#'                  pheno.col = 1,
#'                  qtl = qtl,
#'                  expression.covariates = expression.covariates,
#'                  qtl.method = "hk",
#'                  nperm = 100)
#' qtl2 = makeqtl(cross, chr = summary(s1)$chr[4:5],
#'                pos = summary(s1)$pos[4:5], what = "prob")
#' test2<-covScanQTL(cross = cross, pheno.col = 1,qtl = qtl2,
#'                  which.epiqtl = 1,
#'                  focalqtl.index = 2,
#'                  expression.covariates = expression.covariates,
#'                  qtl.method = "hk",
#'                  nperm = 100)
#'
#' @import qtl
#' @export
covScanQTL<-function(cross,
                     pheno.col = 1,
                     qtl,
                     addcovar = NULL,
                     intcovar = NULL,
                     which.epiqtl = NULL,
                     focalqtl.index = 1,
                     expression.covariates,
                     qtl.method = "hk",
                     nperm = 0,
                     which.perm = NA,
                     plotit=TRUE,
                     verbose = TRUE){

  if("prob" %in% names(cross$geno[[1]])){
    atr<-attributes(cross$geno[[1]]$prob)
    genotypes<-attr(cross$geno[[1]]$prob,"dimnames")[[3]]
  }else{
    if("draws" %in% names(cross$geno[[1]])){
      atr<-attributes(cross$geno[[1]]$draws)
      tmp<-calc.genoprob(cross)
      genotypes<-attr(tmp$geno[[1]]$prob,"dimnames")[[3]]
    }else{
      stop("run either calc.genoprob or sim.geno first.\n")
    }
  }

  if(verbose) cat("creating the marker covariate dataset\n")

  if(nqtl(qtl)>1){
    ag<-argmax.geno(cross, step = atr$step, error.prob = atr$error.prob,
                    off.end = atr$off.end, map.function = atr$map.function,
                    stepwidth = "fixed")
    gp<-pull.argmaxgeno(ag, include.pos.info=F)
    gp.info<-pull.argmaxgeno(ag, include.pos.info=T,rotate=TRUE)[,1:3]

    marker.names<-sapply(1:nqtl(qtl), function(x)
      gp.info$marker[gp.info$chr == qtl$chr[x] & gp.info$pos == qtl$pos[x]])

    if(!is.null(which.epiqtl)){
      epi.marker.names<-sapply(which.epiqtl, function(x)
        gp.info$marker[gp.info$chr == qtl$chr[x] & gp.info$pos == qtl$pos[x]])
      epicov<-data.frame(gp[,epi.marker.names])
      colnames(epicov)<-epi.marker.names
      intcov<-data.frame(intcovar, epicov)
    }else{
      intcov <- intcovar
    }
    addmar<-marker.names[-focalqtl.index]
    addgencov<-data.frame(gp[,addmar])
    colnames(addgencov)<-addmar
    addcov<-cbind(addcovar, addgencov)
  }else{
    addcov <- addcovar
    intcov <- intcovar
  }

  focal.chr<-qtl$chr[focalqtl.index]
  focal.pos<-qtl$pos[focalqtl.index]

  if(verbose) cat("running covariate scan\n")

  bs<-scanone(cross, pheno.col = pheno.col, method=qtl.method,
              addcovar = addcov, intcovar = intcov, chr = focal.chr)

  if(plotit) plot(bs, col = "darkred", main = pheno.col, type = "n",
                  xlab = paste("Chr",chr,"Map position (cM)"))

  wh<-which.min(abs(bs$pos - focal.pos))
  bsl<-as.numeric(bs[wh,-c(1:2)])

  maxScans<-apply(expression.covariates,2, function(x){
    covTemp<-cbind(addcov, x)
    cs<-scanone(cross, pheno.col = pheno.col,
                method=qtl.method,
                addcovar = covTemp,
                intcovar = intcov,
                chr = focal.chr)
    if(plotit) plot(cs, add = T, col = rgb(0,0,0,.2), lty=1)
    return(as.numeric(cs[wh,-c(1:2)]))
  })

  if(plotit) plot(bs, col = "red", main = pheno.col, lty = 2,add = T)

  diffScans<- bsl-maxScans
  rankScans<-rank(-diffScans)

  if(nperm>0){
    if(verbose) cat("running permutation: ")
    permScans<-sapply(1:nperm, function(x){
      if(verbose) if(x %% 10 == 0) cat(x,"")
      apply(expression.covariates,2, function(y){
        covTemp<-cbind(addcov,sample(y))
        cs<-scanone(cross, pheno.col = pheno.col, method=qtl.method,
                    addcovar = covTemp, intcovar = intcov, chr = focal.chr)
        return(as.numeric(cs[wh,-c(1:2)]))
      })
    })
    ps<-sapply(names(maxScans),
               function(x) sum(permScans[x,]<=maxScans[x])/nperm)
    if(verbose) cat("\ndone")
  }else{
    ps<-NA
  }

  out<-data.frame(phenotype = pheno.col,
                  candidateID = names(maxScans),
                  lodAtPeak = maxScans,
                  diffAtPeak = diffScans,
                  rank = rankScans,
                  perm.p = ps,
                  stringsAsFactors=F)
  out<-out[order(out$rank),]
  rownames(out)<-NULL
  return(out)
}
