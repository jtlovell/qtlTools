#' @title Make a scanone plot of genotype means.
#'
#' @description
#' \code{meanScan} Like qtl:effectscan, but plots the genotypic means across the
#' marker / pseudomarker grid.
#'
#' @param cross The qtl cross. Must contain genotype probabilities, calculated via calc.genoprob.
#' If the cross has an X chromosome, this is dropped, as the genotypes are different.
#' @param pheno.col Character or numeric vector indicating the phenotype to be tested.
#' Only 1 phenotype can be tested at a time.
#' @param covar dataframe of covariates with names that match terms in the formula.
#' Each column must be either a character or factor. If a numeric vector was used to fit the
#' model, convert it to a factor by as.factor.
#' @param prob.thresh The genotype probability threshold required to call a genotype. If set at
#' .5 (default) all individuals are assigned genotype calls, otherwise, those with a probability
#' < prob.threshold are called as NA. Values closer to 1 are more stringent.
#' @param leg.pos The character of numeric (length 2) position of the legend. Default is "topright".
#' @param leg.inset The numeric inset of the legend from the plot border.
#' @param leg.bty The border setting for the legend
#' @param cols The vector of colors for the lines in the plot (by genotype)
#' @param ltys The vector of linetypes
#' @param set.mfrow Logical, Should the plotting window be set to the number of covariates?
#' @param ... additional arguments passed on to plot.scanone
#' @details Calculates genotypic means at each marker/psuedomarker, plots them.
#' @return The plot
#'
#' @examples
#' \dontrun{
#' library(qtlTools)
#' library(qtl)
#' data(fake.bc)
#' cross<-fake.bc
#' cross<-calc.genoprob(cross)
#' meanScan(cross, covar = NULL, ylim = c(5,7))
#'
#' # Some simple customization
#' meanScan(cross, covar = NULL, cols = c("black","orange"), ltys = c(1,2),
#'    leg.pos = "top", leg.inset = 0.02)
#'
#' sex<-data.frame(sex = ifelse(pull.pheno(cross, "sex") == 0,"F","M"))
#' meanScan(cross, covar = sex, chr = c(2,7))
#'
#' data(fake.f2)
#' cross<-fake.f2
#' cross<-calc.genoprob(cross)
#' meanScan(cross, covar = NULL)
#' sex<-data.frame(sex = ifelse(pull.pheno(cross, "sex") == 0,"F","M"))
#' meanScan(cross, covar = sex, chr = c(1,13))
#'
#' data(fake.4way)
#' cross<-fake.4way
#' cross<-calc.genoprob(cross)
#' meanScan(cross, covar = NULL)
#' sex<-data.frame(sex = ifelse(pull.pheno(cross, "sex") == 0,"F","M"))
#' meanScan(cross, covar = sex, chr = c(2,7))
#' }
#'
#' @import qtl
#' @export
meanScan<-function(cross, pheno.col = 1,  covar = NULL,
                   prob.thresh = 0, set.mfrow=FALSE,
                   leg.pos = "topright", leg.inset = 0.001, leg.bty = "n",
                   cols = NULL, ltys = NULL, ylim = NULL, ...){

  if("X" %in% chrnames(cross)){
    cat("dropping X chromosome")
    cross<-subset(cross, chr = chrnames(cross)[chrnames(cross)!="X"])
  }
  if (!("prob" %in% names(cross$geno[[1]])))
    stop("You must first run calc.genoprob.")

  if(length(pheno.col) != 1)
    stop("pheno.col (pheno.col) must be a numeric or character vector of length 1.")
  if(is.numeric(pheno.col)){
    pheno.col<-phenames(cross)[pheno.col]
  }
  # 1. get the covariate data in order
  if(is.null(covar)){
    covar<-factor(rep(1, nind(cross)))
    covar.name<-levels(covar)
    covar.id<-NULL
  }else{
    if(class(covar) == "data.frame"){
      covar.id<-names(covar)
      if(ncol(covar)!=1) stop("only single covariates permitted")
      covar<-as.factor(covar[,1])
    }else{
      covar<-as.factor(covar)
      covar.id<-"covar"
    }
    covar.name<-levels(covar)
  }

  # 2. Extract the genotype probabilities
  s1<-scanone(cross, method = "hk")
  m<-makeqtl(cross, chr = s1$chr,
             pos = s1$pos,
             what="prob")

  # 3. infer the genotypes for each individuals * (pseudo)marker
  gp<-lapply(m[[1]], function(x) apply(x,1, function(y) {
    if(max(y) < prob.thresh){
      return(NA)
    }else{
      if(length(colnames(x)[which(y==max(y))])>1){
        return(NA)
      }else{
        return(colnames(x)[which(y==max(y))])
      }
    }
  }))
  gp<-data.frame(do.call(cbind,gp))
  for(i in colnames(gp)) gp[,i]<-as.character(gp[,i])

  # 4. add in phenotype and covariate data
  y<-pull.pheno(cross, pheno.col=pheno.col)
  if(set.mfrow) par(mfrow = c(length(unique(covar)),1))

  # 5. Loop through the covariates, plotting the effects for each.
  for(i in covar.name){
    tmp<-gp[covar == i,]
    tmp.y<-y[covar == i]
    out<-data.frame(t(apply(tmp,2,function(x){
      tapply(tmp.y, x, mean, na.rm=T)
    })))
    if(is.null(covar.id)){
      covar.title = NULL
    }else{
      covar.title = paste0(covar.id,": ", i)
    }

    nnames<-ncol(out)
    mnames<-colnames(out)
    cnames<-1:(nnames)

    for(j in mnames)  s1[,j]<-out[,j]

    if(is.null(ylim)){
      ylim = c(min(as.matrix(out[,cnames])),
               max(as.matrix(out[,cnames])))
    }
    plot(s1, type = "n", ylim = ylim,
         main = covar.title, ylab = pheno.col, ...)
    s1$lod<-NULL

    if(is.null(cols)){
      cols<-rainbow(nnames)
    }
    if(length(cols)!=nnames) stop("vector of colors must match the number of genotypes\n")
    if(is.null(ltys)){
      ltys<-rep(1, nnames)
    }
    if(length(ltys)!=nnames) stop("vector of linetypes must match the number of genotypes\n")
    for(i in mnames) s1[,i]<-out[,i]
    for(i in 1:nnames){
      plot(s1, lodcolumn = cnames[i] , col = cols[i], add=T, lty=ltys[i],...)
    }

    legend(leg.pos, mnames, col = cols, lty=ltys, inset = leg.inset, bty=leg.bty)
  }
}
