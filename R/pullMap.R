#' @title Pull the genetic map and marker names
#'
#' @description
#' \code{pullMap} An extension of qtl::pull.map
#' @param cross The qtl cross object.
#'
#' @import qtl
#' @export
pullMap<-function(cross){
  map <- pull.map(cross, as.table = TRUE)
  map <- data.frame(marker.name = row.names(map),
                  map, stringsAsFactors = FALSE)
  rownames(map) <- NULL
  return(map)
}
