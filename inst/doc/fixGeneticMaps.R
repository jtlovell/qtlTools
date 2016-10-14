## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.width = 5, fig.height = 5, warning = FALSE)
library(knitr)
library(qtl)
library(qtlTools)

## ----get qtlTools, eval = FALSE------------------------------------------
#  library(qtl)
#  
#  library(devtools)
#  install_github("jtlovell/qtlTools")
#  library(qtlTools)

## ----make a messy map----------------------------------------------------
set.seed(42)
map<-sim.map(len = c(50,50,20,30), n.mar = c(25, 10, 10, 50), include.x=F)
plot(map)
cross0<-sim.cross(map, n.ind=50, type="f2",
          error.prob=0.001, missing.prob=0.001, map.function="kosambi")
plot.rf(cross0)

##########
jitterMarkerOrder<-function(cross, chr){
  mars<-1:nmar(cross)[chr]
  set.seed(42)
  badorder<-order(jitter(mars, factor = 10))
  cross<-switch.order(cross, chr = chr, order = badorder, 
                      error.prob=0.001, map.function="kosambi")
}
##########
cross1<-cross0
for(i in 1:nchr(cross1))   cross1<-jitterMarkerOrder(cross=cross1, chr = i)
plot.rf(cross1)

newmap<-est.map(cross1, error.prob=0.001, map.function="kosambi")
cross1 <- replace.map(cross1, newmap)
cross1<-est.rf(cross1)

## ----findSimilar markers-------------------------------------------------
cross2<-dropSimilarMarkers(cross1, rf.threshold = 0.03)
plot.map(cross1, cross2, main = "comparison of full and culled maps")

## ----repRipple markers---------------------------------------------------
cross3<-repRipple(cross2, error.prob=0.001, map.function="kosambi",window = 6)
plot.rf(cross2, main = "recombination fractions before ripple")
plot.rf(cross3, main = "recombination fractions after ripple")

