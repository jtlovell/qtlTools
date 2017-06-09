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
                   covar = NULL, chr = NULL,
                   leg.pos = "topright", leg.inset = 0.001, leg.bty = "n",
                   cols = NULL, ltys = NULL, ylim = NULL, ylab =NULL,
                   plotit = TRUE, draw.legend = TRUE, ...){
  if(!is.null(chr)) cross<-subset(cross, chr = chr)
  if(!pheno.col %in% phenames(cross)) pheno.col<-phenames(cross)[pheno.col]
  if(is.null(covar)){
    formula = paste0(pheno.col, " ~ QTL")
  }else{
    formula = paste0(pheno.col, " ~ ", "QTL + ", paste(colnames(covar), collapse = " + "))
  }

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

  s1<-scanone(cross,method = ifelse("prob" %in% names(cross$geno[[1]]),"hk","imp"))[,-3]
  mars<-gsub("[[:punct:]]", "", colnames(gen))
  rownames(s1)<-mars
  colnames(gen)<-mars

  dat<-cbind(pull.pheno(cross, pheno.col), gen, covar)
  colnames(dat)[1]<-pheno.col

  out<-lapply(mars, function(x){
    form <- as.formula(gsub("QTL",x,formula, fixed=T))
    a<-aggregate(form, data = dat, mean)
    res<-all.vars(form)[1]
    of<-all.vars(form)[-1]
    of<-of[of!=x]
    a$mar.id = paste("geno",genotypes[a[,x]], sep = ":")
    if(length(of)>=1){
      a$covar.id<-sapply(1:nrow(a), function(y) paste(paste(of, collapse = "."),
                                                      paste(a[,of][y], collapse = "."),sep = ":"))
      a$ids<-paste(a$mar.id, a$covar.id, sep = " ")
    }else{
      a$ids<-a$mar.id
    }
    num<-a[,res]
    names(num)<-a$ids
    return(num)
  })
  n <- max(sapply(out, length))
  out1 <- do.call(rbind, lapply(out, `[`, seq_len(n)))
  for(i in colnames(out1)) s1[,i]<-out1[,i]
  if(is.null(ylim)) ylim = c(min(out1, na.rm=T),max(out1, na.rm=T))
  if(is.null(cols)) cols = highContrastColors(ncol(s1))
  if(is.null(ltys)) ltys = rep(1, ncol(s1))
  if(is.null(ylab)) ylab = paste0(pheno.col," mean")
  if(plotit){
    plot(s1, type = "n", ylim = ylim, ylab = ylab, ...)
    for(i in 1:(ncol(s1)-2)) plot(s1, lodcolumn = i, col = cols[i],lty=ltys[i], add = T, ...)
    if(draw.legend){
      legend(leg.pos, inset = leg.inset, colnames(out1), col = cols, lty = ltys, bty = leg.bty)
    }
  }
  return(s1)
}
