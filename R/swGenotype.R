
#' @title Call markers based on allelic imbalance.
#'
#' @description
#' \code{swGenotype} Takes a 2-dimensional dataset (markers x libraries) with
#' the proportion of reference allele calls (for bi-allelic markers) and returns a matrix of genotype calls
#' across a specified sliding window width. Requires 'zoo' to do the sliding window.
#'
#' @param reference.prop The dataset containing ratios of reads to reference genotype.
#' All data must be numeric.
#' @param chr A vector (matching rows in reference.prop) with the chromosome identifier
#' @param pos A vector (matching rows in reference.prop) with the marker position
#' @param marker.id A vector (matching rows in reference.prop) with the marker identifier
#' @param width The width of the window within which calls are made.
#' @param prop2call The threshold of majority voting to make soft calls.
#' If the proportion of genotypes in agreement are < prop2call, NA is returned
#' @param plot.diagnostics Should histograms of reference.prop be returned for each library.
#' Helps to look at this if you are concerned about a few library's quality
#' @param verbose Should updates be printed?
#' @param calcMeans Should mean values for each window be returned? Helps with diagnostics.
#' @param makeSoftCalls Should soft-calls be made? See prop2call.
#' @param makeHardCalls Should hard-calls be made? Hard calls are the majority vote for each
#' window, disregarding prop2call.
#'
#' @details The approach operates library-by-library and ignores among library correlations.
#' For each library, prior thresholds for each genotype are calculated based on the least common
#' observed ratio above 60% and below 40% of the reference genotype. Markers with proportions
#' betweent these thresholds are initially called as A/B (0.5), below the threshold are A/A (0.0)
#' and above the upper threshold is B/B (1.0). Given these initial calls, a sliding window
#' is applied counting the proportion of each. Mean, soft- and hard-calls are made within each
#' window
#'
#' @return  A list containing the desired output.
#'
#'
#' @import qtl
#' @import zoo
#' @export
swGenotype<-function(reference.prop, chr, pos, marker.id,
                     width = 20, prop2call = .95, plot.diagnostics = FALSE, verbose = TRUE,
                     calcMeans = TRUE, makeSoftCalls = TRUE, makeHardCalls = TRUE){
  if(any(duplicated(marker.id))){
    stop("all marker.id's must be unique\n")
  }
  if(length(chr) != nrow(reference.prop) | length(pos) != nrow(reference.prop) | length(marker.id) != nrow(reference.prop)){
    stop("chr, pos and marker id must have the same length as the number of rows in reference.prop\n")
  }
  if(verbose) cat("running sliding-window allele-calling for",nrow(reference.prop),"markers and",ncol(reference.prop),"libraries:\n")
  if(verbose) cat("\t -window width = ",width,"\n")
  if(verbose & makeSoftCalls) cat("\t -% markers in agreement to make call = ",prop2call*100,"\n...\n")
  dat<-data.frame(chr = chr, pos = pos, reference.prop, stringsAsFactors = F)
  rownames(dat)<-marker.id
  byChr.list<-split(dat,as.factor(dat$chr))
  if(verbose) cat("calculating allelic imbalance for each locus\n")
  thresh<-apply(reference.prop,2, function(x){
    his<-data.frame(hist(x, breaks = 100, plot = F)[c("counts","mids")])
    his1<-his[his$mids<0.4,]
    his2<-his[his$mids>0.6,]
    thr1<-his1$mids[which.min(his1$counts)]
    thr2<-his2$mids[which.min(his2$counts)]
    if(plot.diagnostics){
      hist(x, breaks = 100, plot = T, xlab = "proportion of calls", main = "")
      abline(v = c(thr1,thr2), col = "red",lty =2)
      text(y = max(his[,1]), x = thr1, adj = c(0,1), label = "fewest obs <.4")
      text(y = max(his[,1])*.75, x = thr2, adj = c(1,0), label = "fewest obs >.6")
    }
    return(c(thr1,thr2))
  })

  ##
  softcall<-function(z, thrs){
    ifelse(sum(z>thrs[2])/width>=prop2call,1,
           ifelse(sum(z<thrs[1])/width>=prop2call,0,
                  ifelse(sum(z>thrs[1] & z<thrs[2])/width>=prop2call,.5,NA)))
  }
  hardcall<-function(z,thrs){
    c(.5,1,0)[which.max(c(sum(z>=thrs[1] & z<=thrs[2], na.rm = T),
                          sum(z>thrs[2], na.rm = T),
                          sum(z<thrs[1], na.rm = T)))[1]]
  }
  ##

  out<-list()
  if(calcMeans){
    if(verbose) cat("sliding window analysis of mean allelic imbalance for each locus\n")
    meansByChr<-do.call(rbind, lapply(byChr.list, function(x)
      sapply(colnames(reference.prop), USE.NAMES = T, function(y)
        rollapply(x[,y], width = width,  partial = T, na.rm =T, mean))))
    means<-data.frame(marker.id, chr, pos, meansByChr, stringsAsFactors = F)
    out[["means"]]<-means
  }

  if(makeSoftCalls){
    if(verbose) cat("sliding window analysis, making soft-calls\n")
    callsByChr<-do.call(rbind, lapply(byChr.list, function(x)
      sapply(colnames(reference.prop), USE.NAMES = T, function(y)
        rollapply(x[,y], width = width, partial = T, softcall, thrs = thresh[,y]))))
    softCalls<-data.frame(marker.id, chr, pos, callsByChr, stringsAsFactors = F)
    out[["softCalls"]]<-softCalls
  }

  if(makeHardCalls){
    if(verbose) cat("sliding window analysis, making hard-calls\n")
    hardCallsByChr<-do.call(rbind, lapply(byChr.list, function(x)
      sapply(colnames(reference.prop), USE.NAMES = T, function(y)
        rollapply(x[,y], width = width, partial = T, hardcall, thrs = thresh[,y]))))
    hardCalls<-data.frame(marker.id, chr, pos, hardCallsByChr, stringsAsFactors = F)
    out[["hardCalls"]]<-hardCalls
  }
  return(out)
  cat("done\n")
}
