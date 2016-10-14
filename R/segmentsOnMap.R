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
#' @param jColors A dataframe with two columns, "phe" and "color". The phe column contains a all unique
#' phenotypes and the color is the associated color. The default sets the colors to the rainbow color
#' palette.
#' @param ... Other arguments passed to segments
#'
#' @details Pass output from bayesint, lodint, or another confidence
#' interval estimation program to visualize this.
#'
#' @examples
#' data(multitrait)
#'
#' cross <- multitrait
#' cross <- calc.genoprob(cross)
#' phes <- phenames(cross)
#' s1 <- scanone(cross, pheno.col = phes, method = "hk")
#' perms <- scanone(cross, pheno.col = phes, method = "hk",
#'    n.perm = 100, verbose=F)
#' cis<-calcCis(s1.output = s1, perm.output = perms)
#' with(cis, segmentsOnMap(cross, phe = pheno, chr = chr, l = lowposition,
#'    h = highposition, legendCex = .5))
#' segmentsOnMap(cross, calcCisResults = cis, legendCex = .5)
#' segmentsOnMap(cross, calcCisResults = cis, legendCex = .5, orderByLOD = FALSE)
#' segmentsOnMap(cross, calcCisResults = cis, legendCex = .5, col = terrain.colors(length(unique(cis$phe))))
#' segmentsOnMap(cross, calcCisResults = cis, legendCex = .5, max.lwd=6, min.lwd=.5)
#' segmentsOnMap(cross, calcCisResults = cis, legendCex = .5, max.lwd=4, min.lwd=2)
#' segmentsOnMap(cross, calcCisResults = cis, legendCex = .5, lwd = 2)
#' segmentsOnMap(cross, calcCisResults = cis, legendCex = .5, palette = rainbow)
#'
#' @return The plot
#'
#' @import qtl
#' @export

segmentsOnMap<-function(cross, phe, chr, l, h, peaklod = NA, calcCisResults=NULL,
                        legendPosition = "bottom", legendCex = 0.8, col = NULL,
                        palette = highContrastColors, lwd = "byLod", max.lwd = 5, min.lwd = 1,
                        leg.inset = 0.01, orderByLOD = TRUE, ...){

  ############
  # 1. Combine the results into a dataframe
  if(!is.null(calcCisResults)){
    dat<-calcCisResults[,c("pheno","chr","maxLod","lowposition","highposition")]
    colnames(dat)<-c("phe","chr","lod","l","h")
  }else{
    dat<-data.frame(phe = phe, chr = chr, lod = peaklod, l = l, h = h, stringsAsFactors=F)
  }
  dat$phenonum<-as.numeric(as.factor(dat$phe))
  dat$col<-NA

  for(i in c("chr","l","lod","h")) dat[,i]<-as.numeric(as.character(dat[,i]))

  ############
  # 2. Get the colors in order
  if(!is.null(col)){
    if(!any(length(col) == max(dat$phenonum), length(col) == 1))
      stop("provide a color vector of length 1 or with the same length as unique phenotypes\n")
    if(length(col) == 1) col<-rep(col, nrow(dat))
    for(i in 1:max(dat$phenonum)){
      dat$col[dat$phenonum == i]<-col[i]
    }
  }else{
    cols<-palette(max(dat$phenonum))
    for(i in 1:max(dat$phenonum)){
      dat$col[dat$phenonum == i]<-cols[i]
    }
  }

  ############
  # 3. Get the line weights in order
  if(!any(length(lwd) == 1 & is.numeric(lwd) | lwd == "byLod"))
    stop("lwd, if specified, must be a numeric vector of length 1\n or with the same length as the numbe of unique phenotypes\n")
  if(lwd == "byLod"){
    dat$lwd<-log2(dat$lod)
    dat$lwd[dat$lwd>5]<-5
    dat$lwd<-dat$lwd-min(dat$lwd)
    dat$lwd<-dat$lwd/(max(dat$lwd))
    dat$lwd<-dat$lwd*(max.lwd-min.lwd)
    dat$lwd<-dat$lwd+min.lwd
  }else{
    if(length(lwd) == 1){
      dat$lwd<-rep(lwd, nrow(dat))
    }
  }
  dat.ci<-dat

  ############
  # 4. Plot the map
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

  ## 5. Figure out how tightly to pack the segments
  max.nq<-max(table(dat.ci$chr))
  min.st<-scl
  max.st<-1-sepdist-scl
  compress<-max.nq/(max.st-min.st)

  ### Add confidence interval segments
  for(i in chrns){
    if(i %in% dat.ci$chr){
      tem<-dat.ci[dat.ci$chr == i,]
      if(nrow(tem)==1){
        tem$x=0
      }else{
        if(orderByLOD){
          tem<-tem[order(-tem$lod),]
        }else{
          tem<-tem[order(tem$l),]
        }
        tem$x<-c(0,sapply(2:nrow(tem),function(x){
          a<-tem[x,]
          o<-tem[1:(x-1),]
          ho<-o$h>=a$l & o$h<=a$h
          lo<-o$l>=a$l & o$l<=a$h
          bo<-o$l<=a$l & o$h>=a$h
          no<-o$l>=a$l & o$h<=a$h
          sum(sapply(1:length(ho), function(y){
            any(ho[y],lo[y],no[y],bo[y])
          }))
        }))
      }
      xs<-0:max(tem$x)
      if(!all(xs %in% tem$x)){
        wh<-xs[-which(xs %in% tem$x)]
        wh2<-which(tem$x>max(wh))
        tem$x[wh2]<-tem$x[wh2]-length(wh)
      }
      tem$x<-(tem$x/compress)+i+min.st
      with(tem, segments(x0 = x, x1=x, y0 = l, y1=h, col = col, lwd=lwd))
    }
  }

  ### Add legend
  if(!is.null(legendPosition)){
    leg.dat<-dat.ci[!duplicated(dat.ci$phe),]
    leg.dat<-leg.dat[order(leg.dat$phenonum),c("phe","col","lwd")]
    leg.dat$lwd = leg.lwd
    if(lwd == "byLod"){
      leg.dat<-rbind(leg.dat,
                     with(dat.ci,
                          data.frame(phe = paste0("LOD = ",
                                                  floor(c(min(lod),
                                                          mean(lod),
                                                          max(lod)))),
                                     col = "black",
                                     lwd = c(min(lwd),
                                             mean(lwd),
                                             max(lwd)))))
    }
    with(leg.dat,
         legend(legendPosition, legend=phe, col=col, lty = 1,lwd =lwd, cex = legendCex, bty="n",
                inset = leg.inset))
  }
}
