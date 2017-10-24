#' @title Make a scanone plot of genotype means.
#'
#' @description
#' \code{meanScan} Like qtl:effectscan, but plots the genotypic means across the
#' marker / pseudomarker grid.
#'
#' @param cross The qtl cross. Must contain genotype probabilities or imputed genotypes
#' @param pheno.col Character or numeric vector indicating the phenotype to be tested.
#' Only 1 phenotype can be tested at a time.
#' @param covar dataframe of covariates like addcovar in qtl::scanone.
#' @param chr The chromosomes to analyze. If not supplied, all are tested.
#' @param leg.pos The character of numeric (length 2) position of the legend. Default is "topright".
#' @param leg.inset The numeric inset of the legend from the plot border.
#' @param leg.bty The border setting for the legend
#' @param cols optional vector of line colors.
#' @param ltys optional vector of linetypes
#' @param ylim optional limits for the y axis
#' @param ylab optional label for the y axis
#' @param plotit Logical, should a plot be made?
#' @param draw.legend Logical, should a legend be drawn?
#' @param ... additional arguments passed on to plot.scanone
#' @details Calculates genotypic means at each marker/psuedomarker, plots them.
#' @return The plot (if plotit) and an object of class scanone containing the means
#' for each combination of genotype and covariate.
#'
#' @examples
#' \dontrun{
#' library(qtlTools)
#' library(qtl)
#' data(fake.bc)
#' cross<-fake.bc
#' cross<-calc.genoprob(cross)
#' ms<-meanScan(cross, covar = NULL, ylim = c(5,7), leg.inset = .05)
#'
#' # Some simple customization
#' ms<-meanScan(cross, covar = NULL, cols = c("black","orange"), ltys = c(1,2),
#'    leg.pos = "top", leg.inset = 0.02)
#'
#' sex<-data.frame(sex = ifelse(pull.pheno(cross, "sex") == 0,"F","M"))
#' ms<-meanScan(cross, covar = sex, chr = c(2,7), leg.pos = "right", leg.inset = .1)
#' ms<-meanScan(cross, covar = sex, chr = c(2,7), leg.pos = "right", leg.inset = .1,
#'    cols = c("black","grey","darkblue","lightblue"))
#'
#' cross.fem<-subset(cross, ind = sex[,1] == "F")
#' cross.male<-subset(cross, ind = sex[,1] == "M")
#' par(mfrow = c(2,1))
#' ms <-meanScan(cross.fem, leg.pos = "topright", leg.inset = .05, main = "females only")
#' ms <-meanScan(cross.male, leg.pos = "topright", leg.inset = .05, main = "males only")
#'
#' data(fake.f2)
#' cross<-fake.f2
#' cross<-calc.genoprob(cross)
#' ms <- meanScan(cross, covar = NULL, leg.pos = "bottomright", leg.inset = .1)
#' sex<-data.frame(sex = ifelse(pull.pheno(cross, "sex") == 0,"F","M"))
#' ms<-meanScan(cross, covar = sex, chr = c(1,13))
#'
#' data(fake.4way)
#' cross<-fake.4way
#' cross<-calc.genoprob(cross)
#' ms<-meanScan(cross, covar = NULL)
#' sex<-data.frame(sex = ifelse(pull.pheno(cross, "sex") == 0,"F","M"))
#' ms<-meanScan(cross, covar = sex, chr = c(2,7))
#' }
#'
#' @import qtl
#' @export
meanScan<-function(cross, pheno.col = 1,
                   covar = NULL, chr = NULL, plotit = T, plotOverallMean = T){
  if(!is.null(chr)) cross<-subset(cross, chr = chr)
  if(!pheno.col %in% phenames(cross)) pheno.col<-phenames(cross)[pheno.col]


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
  ag<-argmax.geno(cross, step = atr$step, error.prob = atr$error.prob,
                  off.end = atr$off.end, map.function = atr$map.function,
                  stepwidth = atr$stepwidth)

  gen<-pull.argmaxgeno(ag, include.pos.info=F)
  map<-pull.argmaxgeno(ag, include.pos.info=T, rotate = T)[,1:3]
  
  phe<-pull.pheno(cross, pheno.col)
  genids<-attr(cross$geno[[1]]$prob,"dimnames")[[3]]
  if(!is.null(covar)){
    if(length(covar)!=nind(cross))
      stop("covar must be a vector of length = nind (cross)\n")
    
    means<-data.frame(t(apply(gen,2, function(x){
      tapply(phe, list(x,covar), mean)
    })))
    colnames(means)<-unlist(lapply(unique(covar), function(x) 
      paste0(genids,"_",x)))
  }else{
    means<-data.frame(t(apply(gen,2, function(x){
      tapply(phe, x,mean)
    })))
    colnames(means)<-genids
  }

  
  
  
  out<-cbind(map, means)
  out$marker<-NULL
  class(out)<-c("scanone","data.frame")
  if(plotit){
    cols = c("darkred","darkgrey","cornflowerblue","darkorange")
    colnam<-c("red","grey","blue","orange")
    par(mfrow = c(length(unique(covar)),1))
    for(i in unique(covar)){
      j<-paste0("_",i)
      wh<-grep(j,colnames(out))
      lims = c(min(as.matrix(out[,wh])),
               max(as.matrix(out[,wh])))
      plot(out, lod = wh-2, ylim = lims, ylab = "genotype_means",
           col = cols[1:length(wh)],
           main = paste(paste(colnam[1:length(wh)],colnames(out)[wh], sep = ":"), 
                        collapse = ", "))
      if(plotOverallMean){
        abline(h = colMeans(out[,wh]), lwd = .5, lty = 2, col = cols[1:length(wh)])
      }
    }
  }
  return(out)
}
