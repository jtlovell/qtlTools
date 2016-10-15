## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.width = 7, fig.height = 7, warning = FALSE)
library(knitr)
library(qtl)
library(qtlTools)
library(plyr)

## ----env.set.up, eval = FALSE--------------------------------------------
#  library(devtools)
#  install_github("jtlovell/qtlTools")

## ----env.set.up2---------------------------------------------------------
library(qtlTools)
library(qtl)
library(ggplot2)
library(lsmeans)
library(qtlpvl)

## ------------------------------------------------------------------------
data(multitrait)
cross<-multitrait
cross<-calc.genoprob(cross)
phes<-phenames(cross)

## ------------------------------------------------------------------------
s1<-scanone(cross, pheno.col = phes, method = "hk")
cols<-rainbow(length(phes))
plot(s1, type="n", ylim = c(0,max(as.matrix(s1[,-c(1:2)]))), 
     ylab = "LOD Score", main = "all phenotypes plotted together")
for(i in 1:length(phes)) plot(s1, add=T, lodcolumn = i, col = cols[i])

## ------------------------------------------------------------------------
perms<-scanone(cross, pheno.col = phes, method = "hk", n.perm = 100, verbose=F)
print(pullSigQTL(cross, pheno.col=phes,
                 s1.output = s1,  perm.output = perms, 
                 returnQTLModel=FALSE, alpha = 0.05,
                 controlAcrossCol = TRUE))

## ------------------------------------------------------------------------
mods<-pullSigQTL(cross, pheno.col=phes,
                 s1.output = s1,  perm.output = perms, 
                 returnQTLModel=TRUE, alpha = 0.05,
                 controlAcrossCol = TRUE)

## ------------------------------------------------------------------------
print(mods[1:5])

## ------------------------------------------------------------------------
cis<-calcCis(s1.output=s1, perm.output=perms)
head(cis)

## ------------------------------------------------------------------------
segmentsOnMap(cross, calcCisResults = cis, legendCex = .35)

## ---- warning = FALSE----------------------------------------------------
mod<-mods[["X3.Methylsulfinylpropyl"]]
form<-paste("y ~",paste(mod$altname, collapse = " + "))
mod<-refineqtl(cross, qtl = mod, formula = form, verbose=F, method="hk")
qtlStats(cross, 
         pheno.col = "X3.Methylsulfinylpropyl", 
         form = form, 
         mod = mod)

## ------------------------------------------------------------------------
mods<-mods[sapply(mods,length)!=1] # a NULL QTL model has length = 1, toss these

## ---- warning = FALSE----------------------------------------------------
stats<-lapply(names(mods), function(x){
  mod<-mods[[x]]
  form<-paste("y ~",paste(mod$altname, collapse = " + "))
  mod<-refineqtl(cross, qtl = mod, formula = form, 
                verbose=F, method="hk", pheno.col = x)
  qtlStats(cross, 
         pheno.col = x, 
         form = form, 
         mod = mod)
})
stats<-do.call(rbind, stats)
print(head(stats))

## ------------------------------------------------------------------------
segmentsOnMap(cross=cross, phe=stats$phenotype, 
              chr=stats$chr, 
              l = stats$lowposition, lwd=2,
              h =stats$highposition,legendCex=.35, leg.inset=.01)

