## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.width = 7, fig.height = 7)
library(knitr)
library(qtl)
library(qtlTools)
library(plyr)

## ----env.set.up, eval = FALSE--------------------------------------------
#  library(devtools)
#  install_github("jtlovell/qtlTools", build_vignettes = TRUE)
#  library(qtlTools)
#  library(qtl)

## ------------------------------------------------------------------------
data(multitrait)
cross<-multitrait
cross<-calc.genoprob(cross)

## ------------------------------------------------------------------------
phes <- c("X3.Hydroxypropyl", "X3.Butenyl", "X3.Methylsulfinylpropyl", "X5.Methylsulfinylpentyl", "X6.Methylthiohexyl", "X4.Benzoyloxybutyl", "Quercetin.deoxyhexosyl.dihexoside")

## ------------------------------------------------------------------------
modlist<-list(X3.Hydroxypropyl=makeqtl(cross, chr = c(4,5), pos=c(7,35), what="prob"),
              X3.Butenyl=makeqtl(cross, chr = c(4,5), pos=c(7,35), what="prob"),
              X3.Methylsulfinylpropyl=makeqtl(cross, chr = c(4,5), pos=c(0,35), what="prob"),
              X5.Methylsulfinylpentyl=makeqtl(cross, chr = c(1,5), pos=c(105,40), what="prob"),
              X6.Methylthiohexyl=makeqtl(cross, chr = c(4,5), pos=c(9,35), what="prob"),
              X4.Benzoyloxybutyl=makeqtl(cross, chr = c(5), pos=c(35), what="prob"),
              Quercetin.deoxyhexosyl.dihexoside=makeqtl(cross, chr = c(4,5), pos=c(13,35), what="prob"))

## ------------------------------------------------------------------------
refout<-lapply(phes, function(x){
  refineqtl(cross, qtl=modlist[[x]], pheno.col=x, method="hk", verbose=F)
})
names(refout)<-phes
models.out<-refout

## ------------------------------------------------------------------------
stepout<-lapply(phes, function(i) {
  stepwiseqtl(cross, pheno.col=i, method="hk", max.qtl=4, 
              additive.only=TRUE, penalties = c(2.5,3.5,1.5), verbose=F)
})
names(stepout)<-phes
models.out<-stepout

## ------------------------------------------------------------------------
cis<-lapply(models.out, function(x) {
  calcCis(mod=x, lodint=FALSE,  prob=0.99, expandtomarkers=T)
})

## ------------------------------------------------------------------------
cis.df<-ldply(cis, data.frame)
colnames(cis.df)[1]<-"phe"

