## ----setup, include=FALSE, message=FALSE---------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.width = 7, fig.height = 7, warning = FALSE)
library(knitr)

## ----env.set.up, warning = F, message=FALSE------------------------------
library(devtools)
install_github("jtlovell/qtlTools")
library(qtlTools)

## ----loadmultitrait------------------------------------------------------
data(multitrait)

## ----fakemap-------------------------------------------------------------
map<-pullMap(multitrait)
map$bp<-0
for(i in unique(map$chr)){
  n<-sum(map$chr==i)
  p<-sin((1:n/n)*pi)
  map$bp[map$chr==i]<-cumsum(p*1000000)
}

## ----fakegff-------------------------------------------------------------
gff<-data.frame(chr = rep(paste0("scaffold_",1:5),each = 200),
   feature = rep("gene",1000),
   start = rep(seq(from = 0, to = max(map$bp), length = 200), 5),
   end = rep(seq(from = 0, to = max(map$bp), length = 200))+1000,
   strand = rep("+",1000),
   attribute = paste0("gene",1:1000,";","gene",1:1000,".1"), stringsAsFactors=F)

## ------------------------------------------------------------------------
geneCM<-findGenecM(cross = multitrait, marker.info = map, gff = gff,
   gffCols = c("chr","feature","start","end","strand","attribute"))

## ----plotrecom-----------------------------------------------------------
par(mfrow=c(3,2))
for(i in unique(map$chr)){
  with(geneCM[geneCM$chr==i,], plot(pos,bp, col="grey", 
                                main = "cM and bp positions of genes and markers",
                                ylab = "physical position (bp)",
                                xlab = "mapping position (cM)"))
  with(map[map$chr==i,], points(pos,bp, col=i, pch = 19, cex=.8))
}

## ------------------------------------------------------------------------
s1<-scanone(multitrait, method="hk", pheno.col=1)
perm<-scanone(multitrait, n.perm=100, method="hk",pheno.col=1, verbose=FALSE)
cis<-calcCis(cross = multitrait, s1.output=s1, perm.output=perm, drop=5)

par(mfrow = c(1,1))
plot(s1)
segmentsOnPeaks(multitrait, s1.output=s1, calcCisOutput = cis, int.y = 13.1)

## ------------------------------------------------------------------------
candGenes<-findGenesInterval(findGenecM.output = geneCM, calcCis.output = cis)
print(candGenes)

