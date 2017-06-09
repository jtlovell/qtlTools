#' @title Method to improve a genetic map.
#'
#' @description
#' \code{dropSimilarMarkers} finds markers that have a small recombination fraction and
#' drops the one with combined greater segregation distortion and/or missing data. ***Note:
#' if the cross object has many markers (>1000), avoid running this function on more than one
#' chromosome at a time. Also make sure to run est.rf first and use re.est.map = FALSE***
#'
#' @param cross The qtl cross object to search
#' @param chr numeric, Should the analysis be restricted to a set of chromosomes.
#' Defaults to all chromosomes in the cross object.
#' @param rf.threshold The recombination fraction threshold to drop a marker. If est.rf has
#' not been run on cross, it will be done so automatically. See qtl::est.rf for details
#' @param sd.weight The weighting of segregation distortion rank in dropping a marker.
#' Higher values relative to na.weight increase the weight of the sd rank. Setting a value
#' of 0 removes sd as a factor in choosing the best marker.
#' @param na.weight Same as sd.weight, but for the number of NAs.
#' @param keepEnds Logical, should markers on the ends of the chromosomes always be retained?
#' @param doNotDrop Character vector of markers to retain no matter their rfs.
#' @param verbose Logical, should updates be printed?
#' @param blockSize If not NULL, do an initial culling by splitting markers into
#' blocks of this size. Smaller blocks run more quickly than large blocks, but when the total
#' number of blocks excedes ~ 2000, it can take a very long time to parse the cross object
#' into blocks.
#' @param byChr Should the procedure be run chromosome-by-chromosome. If blocksize != NULL,
#' this procedure is run following block-wise culling. If there are many thousands of markers
#' it is recommended to run multiple block-wise calls prior to whole-chromosome procedures. In
#' general, chromosomes with > 1k markers should first be culled using blockSize != NULL.
#' @param runFullMatrix should the full matrix ever be assessed?
#'
#' @param ... if recombination fractions are not included in the cross object,
#' pass on additional arguments to est.rf.
#'
#' @return A new cross object with the similar markers dropped.
#'
#' @examples
#' set.seed(42)
#' map<-sim.map(len = c(50,20), n.mar = c(20,30), include.x=FALSE)
#' cross0<-sim.cross(map, n.ind=50, type="f2", map.function="kosambi",
#'    error.prob=.01, missing.prob = .05)
#' cross0<-est.rf(cross0)
#' cross1<-dropSimilarMarkers(cross0)
#' cross2<-dropSimilarMarkers(cross0, keepEnds=TRUE)
#' par(mfrow=c(2,1))
#' plot.map(cross0, cross1, main = "comparison of full and culled maps")
#' plot.map(cross0, cross2, main = "comparison of full and culled maps")
#'
#' @import qtl
#' @export

dropSimilarMarkers<-function(cross,
                             chr = NULL,
                             rf.threshold=0.01,
                             sd.weight=1,
                             na.weight=1,
                             keepEnds = FALSE,
                             doNotDrop = NULL,
                             verbose=TRUE,
                             blockSize = 100,
                             byChr = TRUE,
                             runFullMatrix = FALSE,
                             ...){
  loadNamespace("qtl")
  dsm<-function(cross,
                rf.threshold=0.01,
                sd.weight=1,
                na.weight=1,
                keepEnds = FALSE,
                doNotDrop = NULL,
                verbose=TRUE){
    # 1. Get the rfs, geno table and chromosomes in order
    gt<-geno.table(cross)
    if(!"rf" %in% names(cross)){
      cross<-est.rf(cross)
    }
    rf<-pull.rf(cross, what = "rf")
    rf[!upper.tri(rf)]<-1
    diag(rf)<-1

    # 1.1 Fancy calculation of sd.weight
    if(class(cross)[1]=="4way"){
      gt.names<-colnames(gt)[3:6]
      gt$P.value<-apply(gt[,gt.names],1,function(x)
        min(x, na.rm = T)/nind(cross))
    }
    gt$rank.p<-with(gt, rank(rank(-P.value, ties.method = "min")*sd.weight))
    gt$rank.sd<-with(gt, rank(rank(missing, ties.method = "min")*na.weight))
    gt$rank<-with(gt, rank(rank.p+rank.sd))

    # 2. drop the markers to retain from the matrix
    if(!is.null(doNotDrop)){
      dnd.index<-which(colnames(rf) %in% doNotDrop)
      rf<-rf[-dnd.index, -dnd.index]
    }
    if(keepEnds){
      tokeep<-as.character(
        unlist(
          lapply(pull.map(cross),function(x)
            c(names(x)[1], names(x)[length(x)]))))
      ends.index<-which(colnames(rf) %in% tokeep)
      rf<-rf[-ends.index, -ends.index]
    }

    # 3. Loop through the rf matrix, dropping one of the two markers with the lowest
    # recombination fraction.
    nmarstart<-sum(nmar(cross))
    while(min(rf)<rf.threshold){
      worst<-colnames(rf)[which(rf == min(rf, na.rm=TRUE), arr.ind=T)[1,]]
      gtm<-gt[worst,]
      badmars<-rownames(gtm)[which.max(gtm$rank)[1]]

      which.todrop<-which(colnames(rf) == badmars)
      rf<-rf[-which.todrop,-which.todrop]

      cross<-drop.markers(cross, markers = badmars)
    }
    nmarend<-sum(nmar(cross))
    return(cross)
  }

  if(!is.null(chr)) {
    cross = subset(cross, chr = chr)
  }
  if(!is.null(blockSize)){
    if(verbose) cat("initial n markers =", totmar(cross),"\n")
    spl<-split(markernames(cross), ceiling(seq_along(markernames(cross))/blockSize))
    if(length(spl)>2000) warning("breaking cross into > 2k blocks can be very slow\n")
    if(verbose) cat("parsing cross object into blocks\n")
    temp.cross<-newLG(cross, markerList=spl)
    if(verbose) cat("running on", length(spl), "blocks of", blockSize, "markers ... \nblock: ")
    goodMars<-lapply(chrnames(temp.cross), function(x){
      ctp<-ifelse(length(spl)>1000, 100, ifelse(length(spl)>100,10, ifelse(length(spl)>50,5,1)))
      if(which(chrnames(temp.cross) == x) %% ctp == 0) cat(x,"")
      cr<-subset(temp.cross, chr = x)
      if(totmar(cr)>1){
        cr<-dsm(cr, rf.threshold = rf.threshold,
                sd.weight = sd.weight,verbose = FALSE,
                keepEnds = keepEnds,
                doNotDrop = doNotDrop)
      }
      return(markernames(cr))
    })
    if(verbose) cat("\n")
    toKeep<-unlist(goodMars)
    toDrop<-markernames(cross)[!markernames(cross) %in% toKeep]
    cross<-drop.markers(cross, markers = toDrop)
    if(verbose) cat("n markers after block-wise culling:", totmar(cross),"\n")
  }

  if(byChr){
    if(verbose) cat("running on each chromosome\ninitial n markers:", nmar(cross),"\n")
    if(verbose) cat("Chromosome:")
    goodMars<-lapply(chrnames(cross), function(x){
      if(verbose) cat(x,"")
      cr<-subset(cross, chr = x)
      if(totmar(cr)>1){
        cr<-dsm(cr, rf.threshold = rf.threshold,
                sd.weight = sd.weight,verbose = FALSE,
                keepEnds = keepEnds,
                doNotDrop = doNotDrop)
      }
      return(markernames(cr))
    })
    toKeep<-unlist(goodMars)
    toDrop<-markernames(cross)[!markernames(cross) %in% toKeep]
    cross<-drop.markers(cross, markers = toDrop)
    if(verbose) cat("\nn markers after chromosome-wise culling:", nmar(cross),"\n")
  }
  if(runFullMatrix){
    if(verbose) cat("running for the whole matrix of",totmar(cross),"markers\n")
    cross<-dsm(cross, rf.threshold = rf.threshold,
               sd.weight = sd.weight,verbose = FALSE,
               keepEnds = keepEnds,
               doNotDrop = doNotDrop)
    if(verbose) cat("n markers after whole-map culling:",totmar(cross),"\n")
  }
  return(cross)
}
