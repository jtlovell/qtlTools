# eQTLPipe_cisTrans<-function(cross, phe, perms, covar = NULL, intcovar = NULL, method = "cisTrans", alpha = 0.05){
#   permThresh<-quantile(perms, 1-alpha)
#   cis<-makeCiseqtl()
#   cis.stats<-getCiseqtlStats()
#   trans<-scan4trans()
#   trans.stats<-ciFromScan4trans(, threshold = )
#   return(merge(cis.stats, trans.stats, by=phe))
# }
