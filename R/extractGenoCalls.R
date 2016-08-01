extractGenoCalls<-function(cross, model, covar = NULL, threshold = 0){
  temp<-model[[1]]
  for(i in 1:length(temp)) temp[[i]][temp[[i]]<=threshold]<-NA

  gp <- lapply(temp, function(x){
    apply(x, 1, function(y) {
      ifelse(sum(is.na(y)) == length(y),NA, which.max(y))
    })
  })
  gp2 <- data.frame(do.call(cbind, gp))
  colnames(gp2) <- mod$altname
  if(!is.null(covar)){
    gp2<-cbind(covar, gp2)
  }
  return(gp2)
}
