## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.width = 5, fig.height = 5)
library(knitr)
library(qtl)
library(qtlTools)

## ----get qtlTools, eval = FALSE------------------------------------------
#  library(qtl)
#  library(devtools)
#  install_github("jtlovell/qtlTools")
#  library(qtlTools)

## ----sim map-------------------------------------------------------------
# Set the seed so that all our maps are identical
set.seed(42) 

# Simulate a 2-chromosome map with 25 markers on each 100cM linkage group
map<-sim.map(len = c(100,100), n.mar = c(25,25), include.x=F)

# Simulate 1 QTL on each chromosome
# Here, we set the bigger one on Chr 1 to be purely additive
#       the second, smaller one is dominant
cross<-sim.cross(map, n.ind=200, type="f2", map.function="kosambi",
                 model = rbind(c(1,10,10,1),c(2,90,5,5)))

## ------------------------------------------------------------------------
plot.map(cross)

## ------------------------------------------------------------------------
plot.rf(cross)

## ------------------------------------------------------------------------
cross<-calc.genoprob(cross, 
                     step = 2, # This means that we want 
                               # pseudomarkers every 2cM
                     error.prob = 0.001)

## ------------------------------------------------------------------------
meanScan(cross, pheno.col = "phenotype")

## ------------------------------------------------------------------------
geno.image(cross, chr = 1, reorder = FALSE)

## ------------------------------------------------------------------------
geno.image(cross, chr = 1, reorder = TRUE)

## ------------------------------------------------------------------------
set.seed(42) 
map<-sim.map(len = c(100,100), n.mar = c(25,25), include.x=F)
cross<-sim.cross(map, n.ind=200, type="f2", map.function="kosambi",
                 model = rbind(c(1,10,1,1),c(1,90,.5,0),c(2,90,.5,1)))
cross<-calc.genoprob(cross, step=2, error.prob = 0.001)

## ------------------------------------------------------------------------
s1<-scanone(cross, method = "hk")
summary(s1)
plot(s1, main = "one-way QTL scan")

## ------------------------------------------------------------------------
perms<-scanone(cross, n.perm=100, verbose=F, method = "hk")
plot(perms)
abline(v = quantile(perms,.95), col = "red", lty=2)
summary(s1, perms = perms, pvalues=T)

plot(s1, main = "one-way QTL scan with significance threshold")
add.threshold(out=s1, perms=perms, alpha = 0.05, col = "red", lty=2)

## ------------------------------------------------------------------------
maxS1 <- max(s1)
s1chr <- as.numeric(as.character(maxS1$chr))
s1pos <- as.numeric(as.character(maxS1$pos))

## ------------------------------------------------------------------------
qtl1 <- makeqtl(cross, chr=s1chr, pos=s1pos, what="prob")
plot(qtl1)

## ------------------------------------------------------------------------
form1 <- "y ~ Q1"
s2 <- addqtl(cross, qtl=qtl1, formula=form1, method="hk")
plot(s2)
s2chr <- max(s2)$chr
s2pos <- max(s2)$pos
qtl2 <- addtoqtl(cross, qtl=qtl1, chr=s2chr, pos=s2pos)

## ------------------------------------------------------------------------
form2 <- "y ~ Q1 + Q2"
qtl2 <- refineqtl(cross, formula=form2, qtl=qtl2, method="hk", verbose=F)
plotLodProfile(qtl2)

## ------------------------------------------------------------------------
plot(s3 <- addqtl(cross, qtl=qtl2, formula=form2, method="hk"))
plot(qtl3 <- addtoqtl(cross, qtl=qtl2, chr=max(s3)$chr, pos=max(s3)$pos))
form3 <- "y ~ Q1 + Q2 + Q3"
plotLodProfile(qtl.model <- refineqtl(cross, formula=form3, qtl=qtl3, method="hk", verbose=F))

## ------------------------------------------------------------------------
step.model<-stepwiseqtl(cross, max.qtl = 3, method = "hk", additive.only = T)
plotLodProfile(step.model)
plot(step.model)

## ------------------------------------------------------------------------
print(fit<-summary(fitqtl(cross, qtl = step.model, formula = formula(step.model), 
                    dropone = T, get.ests = T, method = "hk")))

## ------------------------------------------------------------------------
lodint(step.model, qtl.index = 1, drop = 1.5) # LOD drop of 1.5
lodint(step.model, qtl.index = 1, drop = 4) # LOD drop of 4
bayesint(step.model, qtl.index = 1, prob = .80)
bayesint(step.model, qtl.index = 1, prob = .99)

## ------------------------------------------------------------------------
print(cis<-calcCis(step.model, lodint = TRUE, drop = 1.5))

## ------------------------------------------------------------------------
cis$phenotype<-"fake"
segmentsOnMap(cross=cross, phe=cis$phenotype, 
              chr=cis$chr, 
              l = cis$lowposition, 
              h =cis$highposition, 
              lwd = 5, legendPosition = "right", leg.inset=.1)

## ------------------------------------------------------------------------
qtlStats(cross,  pheno.col = "phenotype",
         form = formula(step.model), 
         mod = step.model)

## ------------------------------------------------------------------------
library(lsmeans)
alllsms<-lsmeans4qtl(cross, 
                     pheno.col = "phenotype",
                     form = "y ~ Q1 + Q2 + Q3", 
                     mod = step.model, 
                     covar=NULL)
print(alllsms)

## ------------------------------------------------------------------------
alllsms<-lsmeans4qtl(cross, 
                     pheno.col = "phenotype",
                     form = "y ~ Q1 + Q2 + Q3 + Q1*Q3", 
                     mod = step.model, 
                     covar=NULL)
lsms<-alllsms[!is.na(alllsms$Q1) & 
                !is.na(alllsms$Q3),]
lsms<-lsms[,c("Q1","Q3","lsmean","SE","mean","sem")]

library(ggplot2)
pos<-position_dodge(.1)
ggplot(lsms, aes(x = Q1, y = lsmean, shape = Q3,
   color = Q3, group = Q3))+
   geom_point(position = pos)+
   geom_line(position = pos)+
   theme_jtl()+
   geom_errorbar(aes(ymin = lsmean - SE, ymax = lsmean+SE), width = .1,position = pos)+
   ggtitle("sas-style LSMeans")

