library(devtools)
install_github("jtlovell/qtltools")
library(qtlTools)

library(qtl)
source("https://raw.githubusercontent.com/kbroman/qtl/master/R/util.R")
source("https://raw.githubusercontent.com/kbroman/qtl/master/R/calc.pairprob.R")
source("https://raw.githubusercontent.com/kbroman/qtl/master/R/xchr.R")

data(fake.f2)
covar<-data.frame(covar = fake.f2$phe$sex)
fake.f2<-calc.genoprob(fake.f2)

#######################################
# scanone permutation comparison
#######################################
# normal
set.seed(42)
perms0<-scanone(fake.f2, pheno.col="phenotype", addcovar=covar,
                intcovar=covar, perm.strata=covar[,1],
                n.perm=100, verbose = F)
summary(perms0)

#GWERk
set.seed(42)
perms1<-scanone.GWERk(fake.f2, pheno.col="phenotype", addcovar=covar,
                intcovar=covar, perm.strata=covar[,1],
                n.perm=100, GWERk=1, verbose = F)
summary(perms1)

plot(as.numeric(perms0), as.numeric(perms1), xlab="standard perms", ylab = "GWER perms")
abline(a=0,b=1, lty=3)

#######################################
# scantwo permutation comparison
#######################################
# normal
set.seed(42)
perms2.0<-scantwo(fake.f2, chr=c(1,2), pheno.col="phenotype", addcovar=covar,
                intcovar=covar, perm.strata=covar[,1],
                n.perm=50, verbose = T)
summary(perms2.0)

#GWERk
set.seed(42)
perms2.1<-scantwo.GWERk1(fake.f2, chr=c(1,2), pheno.col="phenotype", addcovar=covar,
                       intcovar=covar, perm.strata=covar[,1],
                       n.perm=50, verbose = T)
summary(perms2.1)

plot(as.numeric(perms2.0$full), as.numeric(perms2.1$full), xlab="standard perms", ylab = "GWER perms")
abline(a=0,b=1, lty=3)

#######################################
# stepwiseqtl comparison
#######################################
set.seed(42)
perms2.1<-scantwo.GWERk1(fake.f2, pheno.col="phenotype", addcovar=covar,
                         intcovar=covar, perm.strata=covar[,1],
                         n.perm=10, verbose = T)
summary(perms2.1)
pens = calc.penalties(perms2.1, alpha = 0.05)

outsw <- stepwiseqtl(fake.f2, pheno.col="phenotype", max.qtl=3, method="hk", keeptrace=TRUE,
                     penalties = pens)
summary(outsw)
summary(fitqtl(fake.f2, pheno.col="phenotype", formula=formula(outsw), qtl = outsw), method="hk")


outsw.int <- stepwiseqtl.int(fake.f2, pheno.col="phenotype", max.qtl=3, method="hk", keeptrace=TRUE, covar=covar)
summary(outsw.int)
summary(fitqtl(fake.f2, pheno.col="phenotype", formula=formula(outsw.int), qtl = outsw.int, covar = covar, method="hk"))

