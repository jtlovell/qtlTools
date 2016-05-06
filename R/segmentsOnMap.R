#' @title Visualize the positions of QTL on a genetic map.
#'
#' @description
#' \code{segmentsOnMap} A basic function to plot QTL confidence intervals.
#' Useful for <10 traits. All inputs (excluding cross) must have the same length,
#' since each element represents a single segment on the map.
#'
#' @param cross The qtl cross object with marker names that need to be changed.
#' @param phe Character or numeric vector indicating the phenotype to be tested
#' @param chr Vector of chromosome ids - will be coerced to numeric
#' @param l The lower confidence interval bound for each qtl
#' @param h The upper confidence interval bound for each qtl
#' @param seqSpread How far apart (x axis) the semgents are
#' @param legendPosition Where to place the legend. This is passed to the first argument of
#' the legend function. See legend documentation
#' @param legendCex The character expansion for the legend. This is passed to the cex argument
#' of legend.
#' @param chr.adj How far should the first segment be from the map drawing? Defaults to scale with the
#' number of chromosomes
#' @param ... Other arguments passed to segments
#'
#' @details Pass output from bayesint, lodint, or another confidence
#' interval estimation program to visualize this.
#'
#' @return The plot
#'
#' @import qtl
#' @export

segmentsOnMap<-function(cross, phe, chr, l, h, segSpread = 0.15,
                        legendPosition = "bottom", legendCex = 1, chr.adj=NULL, ...){

  ### Plot the map ###
  if(class(cross)[1]=="4way"){
    chrlens<-chrlen(cross)[1,]
    map<-pull.map(cross, as.table=T)[,1:2]
    colnames(map)<-c("chr","pos")
  }else{
    chrlens<-chrlen(cross)
    map<-pull.map(cross, as.table=T)
  }
  map<-as.data.frame(map, stringsAsFactors=F)
  map$chr <- as.numeric(map$chr)
  chrns<-unique(map$chr)

  plot(chrns, rep(0, nchr(cross)), bty="n",type="n",
       ylim=c(max(chrlens),0),
       xlab="chromosome", ylab = "mapping position (cM)",
       xlim=c(min(chrns), max(chrns)+1),
       xaxt = "n")
  axis(1, at = chrns)
  segments(x0=chrns, x1=chrns, y0=rep(0, nchr(cross)), y1=chrlens)
  scl<-length(chrns)/100
  for(i in chrns){
    dat<-map[map$chr == i,]
    segments(x0 = i-scl, x1= i+scl, y0= dat$pos, y1=dat$pos)
  }

  ### Generate the color distributions
  jColors <- data.frame(phe = unique(phe),
                        color = rainbow(length(unique(phe))))
  cols<- jColors$color[match(phe, jColors$phe)]

  temp<-data.frame(phe=as.factor(phe), chr, l, h)
  temp$qname<-paste(temp$phe, temp$chr, temp$l, sep="_")
  temp<-merge(temp, jColors, by="phe")
  ### Add confidence interval segments
  for(i in chrns){
    if(i %in% chr){
      tem<-temp[temp$chr == i,]
      if(nrow(tem)==1){
       x=0
      }else{
        tem$x<-0
        tem$phe<-as.factor(as.character(tem$phe))
        seqs<-seq(from=min(map$pos[map$chr==i]), to = max(map$pos[map$chr==i]), by = 1)
        out<-seqs
        cmat<-sapply(levels(tem$phe), function(x) {
          for(y in 1:nrow(tem[tem$phe==x,])){
            tem2<-tem[tem$phe==x,][y,]
            out<-ifelse(out>=floor(tem2$l) & out<=ceiling(tem2$h),99999,out)
          }
          out<-ifelse(out==99999,1,0)
        })
        poss<-data.frame(t(apply(cmat, 1, cumsum)))
        poss$index<-seqs
        tem$phecols<-as.numeric(tem$phe)
        tem$x<-sapply(1:nrow(tem), function(x) {
          tem1<-tem[x,]
          (max(with(tem1,poss[poss$index>=floor(l) & poss$index<=ceiling(h),phecols]))-1)*segSpread
        })
      }
      if(is.null(chr.adj)){
        chr.adj<-nchr(cross)*.02
      }
      with(tem, segments(x0 = x+i+chr.adj, x1=x+i+chr.adj, y0 = l, y1=h, col = tem$color, ...))
    }
  }

  ### Add legend
  legend(legendPosition, legend=jColors$phe, col=jColors$color, pch=19, cex = legendCex, bty="n")
}




