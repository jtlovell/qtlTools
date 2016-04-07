#' @title Change marker names
#'
#' @description
#' \code{renameMarkers} Take a cross object and the corresponding genotype matrix
#' and rename a set of specified markers.
#'
#' @param cross The qtl cross object with marker names that need to be changed.
#' @param crossfile The location of the .csv file containing the genotype matrix. Formats accepted
#' are only the .csv and .csvs r/qtl read.cross formats. Other formats may run without error, but
#' do not work.
#' @param oldnames A character vector of the names of markers that need to be changed
#' @param newnames A character vector of the names of markers to replace the current ones
#' @param outputName A character string with the name of the new genotype matrix with new names
#'
#' @return A genotype matrix is written to file in the current working directory, unless a path is
#' specificed in the "crossfile" and "outputName" arguments. The new matrix will have marker names
#' that have been replaced by the newnames provided.
#'
#' @examples
#' library(qtl)
#' data(fake.f2)
#' write.cross(fake.f2, format = "csv", filestem = "test") ### Must be "csv" or "csvs" formats.
#' old.x <- markernames(fake.f2)[grep("X", markernames(fake.f2))]
#' new.x <- c("X1", "X2", "X3")
#' renameMarkers(cross = fake.f2, crossfile = "test.csv",
#'    oldnames = old.x, newnames = new.x, outputName = "new.fake.f2.csv")
#' new.fake.f2 <- read.cross(file = "new.fake.f2.csv", format = "csv")
#' cbind(markernames(fake.f2), markernames(new.fake.f2))
#' @import qtl
#' @export
renameMarkers<-function(cross, crossfile, oldnames, newnames, outputName){
  cat("renaming markers and writing to a new file:", outputName, "\n")
  m <- markernames(cross)
  c <- read.csv(crossfile, header = T, quote = "", na.strings = "NA")
  n <- data.frame(oldnames,newnames)
  for(x in oldnames){
    if(x %in% m){
      colnames(c)[which(colnames(c) == x)] <- as.character(n$newnames[n$oldnames == x])
    }
  }
  c[1:2, ][is.na(c[1:2, ])] <- ""
  write.csv(c, file = outputName, row.names = FALSE, quote = FALSE)
}
