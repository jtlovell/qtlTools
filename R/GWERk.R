Permute the phenotypic trait values among the individuals and generate a set of null statistics; 

repeat B times to get B sets of null statistics. 

To obtain stable estimates, Churchill and Doerge (1994) suggested choosing 
B ≥ 1000 for α = 0.05 and as many as B = 10,000 for more extreme critical values, 
such as α = 0.01.

For each set of null statistics, obtain the peaks statistics and order them. 
Take the (k + 1)th highest peak statistic for each of the B null sets and store them together. 
For example, if we take k = 0, we pick the maximum test statistic for each null set. 
If we choose k = 1, we take the second highest peak statistic for each of the B permutations.
Sort the B (k + 1)th peak statistics in descending order. 
The 100(1 − α) percentile is the estimated GWERk threshold value. 
For example, if we choose B = 1000, k = 1, and α = 0.05, we would take the 50th highest 
value among all 1000 second-highest peak statistics.
setwd("~/Desktop/eqtl_v8")
load("eqtl2015_v8input.RData")
library(qtl)
phes=c("Pahal.A00006", "Pahal.A00007", "Pahal.A00008",
       "Pahal.A00010", "Pahal.A00012", "Pahal.A00015")
cross<-calc.genoprob(cross)
covar=data.frame(covar = as.numeric(ifelse(pull.pheno(cross, pheno.col="Treatment") == "Wet",1,2)))
gks<-GWERk(cross.in=cross, phe=phes, method="hk",addcovar=covar, intcovar=covar, nperms=100)
perms<-scanone(cross, pheno.col=phes, method="hk",addcovar=covar, intcovar=covar, n.perm=100)
perms2<-data.frame(perms)
par(mfrow=c(3,2))
for(i in phes){
  plot(gks[,i][order(gks[,i])], perms2[,i][order(perms2[,i])], 
       ylim=c(0, max(perms[,i])), xlim=c(0, max(gks[,i])), ylab="trad.perms", xlab="GWERk", main=i)
  abline(a=0,b=1)
}
library(plyr)
GWERk<-function(cross.in, phe, k=1, nperms=1000, printUpdate=TRUE, sameSeed = TRUE, ...){
  
  out<-lapply(phe, function(i){
    if(sameSeed) set.seed(42)
    if(printUpdate) cat(i, "\n")
    
    phe.in<-pull.pheno(cross.in, pheno.col=i)
    
    make.perm<-lapply(1:nperms, function(x)  sample(phe.in, size=length(phe.in)))
    
    permd.phe<-data.frame(do.call(cbind,make.perm))
    colnames(permd.phe)<-paste("perm",1:nperms,sep="_")
    
    cross2<-cross.in
    cross2$pheno <- cbind(cross2$pheno, permd.phe)
    
    perm<-scanone(cross2, pheno.col=colnames(permd.phe), ...)
    
    sperm<-summary(perm)[,colnames(permd.phe)]
    sum.perm<-as.numeric(apply(sperm, 2, function(x) {
      for (i in 1:k) x<-x[-which(x==max(x))]
      return(max(x))
    }))
    return(sum.perm)
  })
  names(out)<-phe
  return(do.call(cbind, out))
}