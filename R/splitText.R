#' @title Simplify strsplit functionality
#'
#' @description
#' \code{splitText} Loops along a string splitting text by sep and returning
#' the element num of that split text
#'
#' @param x A character vector, if not character, function coerces it to character.
#' @param sep The character that separates values
#' @param num The index of the split vector to return.
#' @export
splitText<-function(x, sep = "_", num = 1){
  x<-as.character(x)
  sapply(x, function(y) strsplit(x, sep)[[1]][num])
}
