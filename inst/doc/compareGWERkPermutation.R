## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.width = 7, fig.height = 7, warning = FALSE)
library(knitr)
library(qtl)
library(qtlTools)

## ----env.set.up----------------------------------------------------------
library(devtools)
install_github("jtlovell/qtlTools")

## ----env.set.up2---------------------------------------------------------
library(qtlTools)
library(qtl)

## ----get source----------------------------------------------------------
source("https://raw.githubusercontent.com/kbroman/qtl/master/R/util.R")
source("https://raw.githubusercontent.com/kbroman/qtl/master/R/calc.pairprob.R")
source("https://raw.githubusercontent.com/kbroman/qtl/master/R/xchr.R")

## ----load data-----------------------------------------------------------
data(fake.f2)
covar<-data.frame(covar = fake.f2$phe$sex)
fake.f2<-calc.genoprob(fake.f2)

## ----scanone normal------------------------------------------------------
set.seed(42)
perms0<-scanone(fake.f2, pheno.col="phenotype", addcovar=covar,
                intcovar=covar, perm.strata=covar[,1],
                n.perm=100, verbose = F)
summary(perms0)

## ----scanone gwerk-------------------------------------------------------
set.seed(42)
perms1<-scanone.GWERk(fake.f2, pheno.col="phenotype", 
                      addcovar=covar, intcovar=covar, 
                      perm.strata=covar[,1],
                      n.perm=100, GWERk=1, verbose = F)
summary(perms1)

plot(as.numeric(perms0), as.numeric(perms1), 
     xlab="standard perms", ylab = "GWER perms", 
     main = "scanone permutation comparison")
abline(a=0,b=1, lty=3)

