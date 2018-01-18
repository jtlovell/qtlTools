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
#' ms<-meanScan(cross, covar = NULL)
#'
#' # Some simple customization
#' ms<-meanScan(cross, covar = NULL, col = c("black","orange"),
#'    leg.pos = "top")
#'
#' sex<-data.frame(sex = ifelse(pull.pheno(cross, "sex") == 0,"F","M"))
#' ms<-meanScan(cross, covar = sex, leg.pos = "right")
#' ms<-meanScan(cross, covar = sex, chr = c(2,7), leg.pos = "right",
#'    col = c("black","grey","darkblue","lightblue"))
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
#' ms<-meanScan(cross, covar = sex)
#'
#' data(fake.4way)
#' cross<-fake.4way
#' cross<-calc.genoprob(cross)
#' ms<-meanScan(cross, covar = NULL)
#' sex<-data.frame(sex = ifelse(pull.pheno(cross, "sex") == 0,"F","M"))
#' ms<-meanScan(cross, covar = sex, chr = c(2,7))
#'
#' data("multitrait")
#' cross =multitrait
#' cross<-calc.genoprob(cross)
#' ms<-meanScan(cross, covar = NULL)
#' ms<-meanScan(cross, covar = NULL, sw.width = 5)
#' }
#'
#' @import qtl
#' @export
meanScan<-function(cross, pheno.col = 1, chr = NULL, covar = NULL,
                   mean.FUN = function(x) mean(x, na.rm = T),
                   se.FUN = function(x) sd(x, na.rm = T)/sqrt(sum(!is.na(x))),
                   sw.width = 1, plot.se = T, se.alpha = .2,
                   col = NULL, leg.pos = "topright", ...){

  # -- Drop non -autozomes
  clss= sapply(cross$geno, function(x) attr(x, "class"))
  if(any(clss!="A")){
    cat("found non-autozomes, dropping these\n")
    cross = subset(cross, chr = clss=="A")
  }

  # -- Subset cross if chr specified
  if(!is.null(chr)) cross<-subset(cross, chr = chr)

  # -- Get phenotype names in order
  if(is.numeric(pheno.col)){
    pheno.col<-phenames(cross)[pheno.col]
  }
  if(!pheno.col %in% phenames(cross)) stop("phenotype ", pheno.col, " not found\n")

  # -- Find the names of the genotypes
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

  # -- pull the best-guess genotypes
  ag<-argmax.geno(cross, step = atr$step, error.prob = atr$error.prob,
                  off.end = atr$off.end, map.function = atr$map.function,
                  stepwidth = atr$stepwidth)
  gen<-pull.argmaxgeno(ag, include.pos.info=F)
  map<-pull.argmaxgeno(ag, include.pos.info=T, rotate = T)[,1:3]

  # - get the phenotype data
  phe<-pull.pheno(cross, pheno.col)

  if(is.null(covar)){
    covar = data.frame(fake = rep("fake",length(phe)))
  }
  if(nrow(covar)!=nind(cross))
    stop("covar must be a vector of length = nind (cross)\n")
  if(!is.factor(covar[,1])){
    covar[,1] = as.factor(covar[,1])
  }
  cov.levs = levels(covar[,1])

  tmp = tapply(phe,data.frame(gen = as.factor(gen[,1]),covar), mean)

  gen2 = gen
  for(i in 1:length(genotypes)){
    for(j in colnames(gen2)){
      gen2[,j][gen2[,j]==i]<-genotypes[i]
    }
  }

  #for(i in colnames(gen2)) gen2[,i]<-factor(gen2[,i], levels = genotypes)

  ou = tapply(phe, data.frame(gen2[,1],covar), mean.FUN)
  tap.names = apply(expand.grid(colnames(ou),rownames(ou)),1,
                    function(x) paste(x, collapse = "_"))
  tap.names<-gsub("fake_","",tap.names)
  tap.names<-gsub("_fake","",tap.names)
  tap.names<-gsub("fake","",tap.names)

  means<-data.frame(t(apply(gen2,2, function(x){
    tapply(phe, data.frame(x,covar), mean.FUN)
  })))
  colnames(means) <- tap.names

  ses<-data.frame(t(apply(gen2,2, function(x){
    tapply(phe, data.frame(x,covar), se.FUN)
  })))
  colnames(ses) <- tap.names

  hasZoo = requireNamespace("zoo", quietly = TRUE)
  if(sw.width>1 & !hasZoo){
    warning("install the zoo package to conduct sliding window averaging\n")
    sw.width = 1
  }

  if(sw.width>1){
    requireNamespace("zoo")
    means2 = zoo::rollapply(means,width = sw.width, by.column = T,
                            partial = T, fill = NA, FUN = mean.FUN)
    rownames(means2) = rownames(means)
    ses2 = zoo::rollapply(ses,width = sw.width, by.column = T,
                          partial = T, fill = NA, FUN = mean.FUN)
    rownames(ses2) = rownames(ses)
  }else{
    means2 = means
    ses2 = ses
  }

  out.means<-cbind(map, means2)
  out.ses = cbind(map, ses2)
  out.means$marker<-NULL
  out.ses$marker<-NULL
  class(out.means)<-c("scanone","data.frame")
  class(out.ses)<-c("scanone","data.frame")

  uplim = means2+ses2
  lowlim = means2-ses2

  out.uplim<-cbind(map, uplim)
  out.lowlim = cbind(map, lowlim)
  out.uplim$marker<-NULL
  out.lowlim$marker<-NULL
  class(out.uplim)<-c("scanone","data.frame")
  class(out.lowlim)<-c("scanone","data.frame")


  if(is.null(col)){
    col.base = list(a = c("pink","red4"),
                    b = c("cornflowerblue","darkblue"),
                    c = c("lightgrey","black"),
                    d = c("gold","darkorange"),
                    e = c("chartreuse","darkgreen"))
    if(length(genotypes)>5) stop("if n genotypes > 5, manually specify colors")
    pals = lapply(1:length(genotypes), function(x) colorRampPalette(col.base[[x]]))
    cov.levs= length(levels(covar[,1]))
    cols = lapply(pals, function(x) x(cov.levs))
    col = as.character(unlist(sapply(cols, function(x)  sapply(1:cov.levs, function(y) x[y]))))
  }
  if(length(col)!=ncol(means2)) stop("length of color vector must match the number of covariate / genotype combinations\n")


  lims = c(min(as.matrix(lowlim), na.rm = T), max(as.matrix(uplim), na.rm = T))
  me = mean(as.matrix(means2), na.rm = T)
  plot(out.means, type = "n", ylim = lims, ylab = "mean phenotypic value", ...)

  poss = xaxisloc.scanone(out.uplim, thechr = out.uplim[,1], thepos = out.uplim[,2])

  ul = cbind(map$chr, poss, uplim)
  dl = cbind(map$chr,poss, lowlim)


  for(i in 1:ncol(means2)){
    if(plot.se){
      for(j in chrnames(cross)){
        u = ul[ul[,1]==j,]
        d = dl[dl[,1]==j,]
        polygon(x = c(u[,2],rev(d[,2])),
                y = c(u[,i+2], rev(d[,i+2])),
                col = adjustcolor(col[i],alpha.f = se.alpha),
                border = NA)
      }
    }
    plot(out.means, lod = i, col = col[i], add = T)
  }


  legend(leg.pos, colnames(means2), col = col, lwd = 1, cex = .8, bty = "n")

  return(list(mean.scanone = out.means,
              se.scanone = out.ses))
}
