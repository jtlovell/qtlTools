#' @title Orient chromosomes to match some previously found marker order.
#'
#' @description
#' \code{matchMarkerOrder} Use previously known marker orders (e.g. from a
#' annotation) to inform the marker orientation along a chromosome. Does not re-order
#' markers within a linkage group, but only flips orders of chromosomes that seem to be
#' inverted relative to the annotation.
#'
#' @param cross The QTL cross object. If marker names
#' are formated as physicalChromosome_physicalBasePairPosition ("_" separated), this
#' is the only required argument. *Note, the physical Chromosome must match exactly
#' the chromosome names of the cross.*
#' If there are additional characters (e.g. "Chr01" = "1"), specify the character
#' string to sub stripped in marker.char.rm. In this example, marker.char.rm = "Chr0".
#' @param marker.chr.bp If marker names are not as specified above, provide chromosome
#' identities of each marker here. Must match chromosome names of the cross exactly.
#' See qtlTools::newLG to rename the chromosomes in the cross object.
#' @param marker.pos.bp If marker names are not as specified above, provide the physical
#' positions (e.g. bp) of each marker.
#' @param marker.sep If markers are named as chr and pos, but the separating character is
#' not "_", specify here.
#' @param marker.char.rm See above. The character string to drop from the marker name so that
#' it matches the chromosome names exactly.
#' @param plotit If the ggplot2 package is available, should plots be made? Increases
#' run time substantially.
#'
#' @details For each chromosome, fits a linear model and asks if the slope < 0. If so,
#' flips the maker order on that chromosome.
#'
#' @return A cross object with flipped marker order on chromosomes that seem to be in the
#' reverse orientation of the annotation.
#'
#' @examples
#' library(qtlTools)
#' data(fake.f2)
#' cross<-fake.f2
#' \dontrun{
#' fake.f2<-est.rf(fake.f2)
#' cross<-fake.f2
#' #Perturb the marker order and chromosome names
#' markerlist<-lapply(chrnames(cross), function(x) sample(markernames(cross, chr =x)))
#' names(markerlist)<-as.character(chrnames(cross))
#' cross2<-newLG(cross, markerList = markerlist)
#' library(TSP)
#' plot.rf(cross2)
#' cross3<-cross3<-tspOrder(cross = cross2,
#'   hamiltonian = T,
#'   method="nn") # change to your path
#' cross3<-cross3<-tspOrder(cross = cross2,
#'   hamiltonian = T,
#'   method="concorde",
#'   concorde_path = "/Users/John/Documents/concorde/TSP") # change to your path
#' plot.rf(cross3)
#'
#' orig.marker.info<-pull.map(cross, as.table = T)
#' orig.marker.info$marker<-rownames(orig.marker.info)
#' orig.marker.info<-orig.marker.info[match(markernames(cross3),orig.marker.info$marker),]
#'
#' cross4<-matchMarkerOrder(cross3,
#'   marker.chr.bp = orig.marker.info$chr,
#'   marker.pos.bp = orig.marker.info$pos,
#'   plotit=T)
#' }
#' @import qtl
#' @export

matchMarkerOrder<-function(cross, marker.chr.bp = NULL, marker.pos.bp = NULL,
                           marker.sep = "_", marker.char.rm = "Chr0", plotit = F){

  if(!requireNamespace("ggplot2", quietly = TRUE) & plotit){
    warning("install the ggplot2 package to use plotit\n")
    plotit = FALSE
  }else{
    require("ggplot2", quietly = TRUE)
  }

  map<-pullMap(cross)
  map$marker.name<-gsub(marker.char.rm,"", map$marker.name)

  if(is.null(marker.chr.bp)){
    map$chr.orig<-sapply(map$marker.name,
                         function(x) strsplit(x, "_")[[1]][1])
  }else{
    map$chr.orig<-marker.chr.bp
  }
  if(is.null(marker.pos.bp)){
    map$pos.orig<-as.numeric(sapply(map$marker.name,
                                    function(x) strsplit(x, "_")[[1]][2]))
  }else{
    map$pos.orig<-marker.pos.bp
  }
  if(class(cross)[1]=="4way"){
    colnames(map)[3]<-"pos"
    map$pos.male<-NULL
  }

  if(any(!unique(map$chr) %in% unique(map$chr.orig)) |
     any(!unique(map$chr.orig) %in% unique(map$chr)))
    stop("some marker chromosome names do not match the annotation chromosome names\n")

  if(plotit){
    p1<-ggplot(map, aes(x = pos, y = pos.orig, col = chr))+
      geom_point(size=.3)+
      theme(axis.text=element_blank(),
            axis.ticks=element_blank())+
      facet_grid(chr.orig~chr, as.table=F, scale = "free", space="free")+
      labs(x = "mapping position of marker (cM)",
           y = "annotation position of marker (bp)",
           title = "initial position comparison")
  }

  map.same<-map[map$chr == map$chr.orig,]
  flipIt<-sapply(unique(map.same$chr), function(i){
    out<-lm(pos ~ pos.orig,
            data = map.same[map.same$chr==i,])$coefficients["pos.orig"]
    ifelse(out>0,FALSE,TRUE)
  })
  toflip<-chrnames(cross)[flipIt]
  for(i in toflip) cross<-flip.order(cross, chr = i)

  map<-pullMap(cross)
  map$marker.name<-gsub(marker.char.rm,"", map$marker.name)

  if(is.null(marker.chr.bp)){
    map$chr.orig<-sapply(map$marker.name,
                         function(x) strsplit(x, "_")[[1]][1])
  }else{
    map$chr.orig<-marker.chr.bp
  }
  if(is.null(marker.pos.bp)){
    map$pos.orig<-as.numeric(sapply(map$marker.name,
                                    function(x) strsplit(x, "_")[[1]][2]))
  }else{
    map$pos.orig<-marker.pos.bp
  }
  if(class(cross)[1]=="4way"){
    colnames(map)[3]<-"pos"
    map$pos.male<-NULL
  }
  if(plotit){
    p2<-ggplot(map, aes(x = pos, y = pos.orig, col = chr))+
      geom_point(size=.3)+
      theme(axis.text=element_blank(),
            axis.ticks=element_blank())+
      facet_grid(chr.orig~chr, as.table=F, scale = "free", space="free")+
      labs(x = "mapping position of marker (cM)",
           y = "annotation position of marker (bp)",
           title = "marker position comparison, following order matching")
    print(p1)
    print(p2)
  }
  return(cross)
}
