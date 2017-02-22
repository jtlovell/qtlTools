#' @title Drop markers that significantly expand the genetic map.
#'
#' @description
#' \code{dropByDropone} Leverage the results of qtl::droponemarker to find markers
#' that are likely misplaced or erroneous.
#'
#' @param cross The QTL cross object, with an estimated genetic map
#' @param droponeRes Results from qtl::droponemarker
#' @param which.map If the cross is a 4-way, which sex-specific map should be used?
#' See the column names of the dropone results. Defaults to Ldiff.female, but can be
#' set to Ldiff.male. Ignored if cross is not a 4way.
#' @param midMarkerThresh The ldiff (length increase of map) that causes a marker
#' to be dropped from the map.
#' @param endMarkerThresh The ldiff (length increase of map) that causes a marker
#' to be dropped from the map if the marker is at the end/start of a chromosome.
#'
#' @param re.est.map Should the genetic map be re-estimated?
#' @param map.function If re.est.map = TRUE, the map functions to employ
#' @param sex.sp If re.est.map = TRUE, should a sex-specific map be estimated?
#' @param ... Additional arguemnts passed on to est.map if re.est.map = TRUE.
#' @details A simple function to parse output from droponemarker postion.
#'
#' @return A cross object without markers that caused map expansion > thresholds.
#'
#' @examples
#' library(qtlTools)
#' data(fake.f2)
#' cross<-fake.f2
#' \dontrun{
#' dropone<-droponemarker(cross, map.function = "kosambi")
#' plot(dropone, lodcolumn = 2)
#'
#' cross1<-dropByDropone(cross = cross, droponeRes = dropone, endMarkerThresh = 50, re.est.map=T)
#' plot.map(cross,cross1)
#' }
#' @import qtl
#' @export

dropByDropone<-function(cross, droponeRes,
                        endMarkerThresh = 12, midMarkerThresh = 4,
                        which.map = NULL,
                        map.function = "kosambi", re.est.map = F, sex.sp=F,
                        ...){

  if(class(cross)[1] == "4way"){
    if(is.null(which.map)){
      which.map<-"Ldiff.female"
    }
    ldCol = which(colnames(droponeRes)==which.map)
    map<-lapply(pull.map(cross), colnames)
  }else{
    ldCol = which(colnames(droponeRes)=="Ldiff")
    map<-lapply(pull.map(cross), names)
  }

  endMars<-unlist(lapply(map, function(x) c(x[1],x[length(x)])))

  todrop<-rownames(
    droponeRes[droponeRes[,ldCol]>midMarkerThresh,])
  todrop<-todrop[!todrop %in% endMars]
  tail2drop<-rownames(
    droponeRes[rownames(droponeRes) %in% endMars &
                 droponeRes[,ldCol]>endMarkerThresh,])

  cross<-drop.markers(cross, c(todrop,tail2drop))

  if(re.est.map){
    map<-est.map(cross, map.function = "kosambi", sex.sp=F, ...)
    cross<-replace.map(cross, map)
  }

  return(cross)
}
