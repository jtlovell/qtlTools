pal_jtl<-c("black","darkred","cornflowerblue","gold",
           "darkorchid3","lightseagreen","orangered","thistle","blue4","grey60",
           "chartreuse","orange","turquoise1","darkolivegreen3","violet")

pal_jtl3<-colors()[c(418,204,24,
                     100,553,536,
                     621,502,
                     652,513,
                     494,48,81,
                     72,70,431,
                     128,26,73,
                     468,625,
                     509,456,645)]

pal_jtl4<-c(
  RColorBrewer::brewer.pal(9,"Greys")[c(3,6)],
  colors()[c(536,404)],RColorBrewer::brewer.pal(9,"Reds")[c(6,8)],
  colors()[c(622,92)],
  colors()[c(411)],"gold",
  colors()[c(259,48)],
  colors()[c(139,81)],
  colors()[c(69,114)],
  colors()[c(27,73)],
  colors()[c(468,549)],
  colors()[c(509,456)]
)

pal_jtl4_pairs<-function(ncols,
                         nogrey = FALSE,
                         orderCols = TRUE){
  x<-split(pal_jtl4,rep(1:11,each = 2))
  ord <- c(1,2,9,5,8,7,10,3,6,4,11)
  if(nogrey){
    ord<-ord[-1]
  }
  tokeep<-ord[1:ncols]
  names(x)<-1:11
  if(orderCols){
    return(unlist(x[names(x) %in% tokeep]))
  }else{
    return(unlist(x[tokeep]))
  }
}
