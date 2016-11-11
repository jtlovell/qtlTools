#' @title Remove close QTL from stepwiseqtl model selection
#'
#' @description
#' \code{dropCloseQTL}  Assists in simplifying models following stepwise model selection
#' @import qtl
#' @export
dropCloseQTL<-function(cross, qtl, formula, covar, pheno.col, method = "hk",
                       min.dist = 30){
  anyClose<-any(unlist(sapply(chrnames(cross), function(i){
    pos.chr = qtl$pos[qtl$chr == i]
    if(length(pos.chr)>1){
      min(dist(pos.chr))<min.dist
    }
  })))
  while(anyClose){
    fout<-data.frame(summary(fitqtl(cross, qtl = qtl, formula = formula(qtl),
                                    covar = covar, method = "hk",
                                    pheno.col = pheno.col))$result.drop)
    fout$id = rownames(fout)
    for(i in 1:length(qtl$name)) fout$id = gsub(qtl$name[i],qtl$altname[i],fout$id, fixed = T)
    bads<-unlist(sapply(chrnames(cross), function(i){
      pos.chr = qtl$pos[qtl$chr == i]
      if(length(pos.chr)>1){
        if(min(dist(pos.chr))<min.dist){
          todrop = qtl$altname[qtl$chr == i]
          todrop<-todrop[which.min(fout$LOD[fout$id %in% todrop])]
          return(todrop)
        }
      }
    }))

    if(length(bads)>0){
      new.qtl<-dropfromqtl(qtl, qtl.name = bads)

      dat<-merge(data.frame(name = qtl$altname,
                            name2 = qtl$name, stringsAsFactors=F),
                 data.frame(new.name = new.qtl$altname,
                            name2 = new.qtl$name, stringsAsFactors=F),
                 by = "name2")

      fout$new.id<-row.names(fout)
      for(i in 1:nrow(dat)) fout$new.id<-gsub(dat$name2[i],dat$new.name[i], fout$new.id, fixed = T)
      fout<-fout[!grepl("@",fout$new.id, fixed = T),]

      newform2 = paste("y ~ ", paste(fout$new.id, collapse = " + "))
      attr(new.qtl, "formula")<-newform2
      qtl = refineqtl(cross, pheno.col = pheno.col, qtl = new.qtl, formula = newform2, method = "hk", covar = covar,
                      verbose = F)
    }
    anyClose<-any(unlist(sapply(chrnames(cross), function(i){
      pos.chr = qtl$pos[qtl$chr == i]
      if(length(pos.chr)>1){
        min(dist(pos.chr))<min.dist
      }
    })))
  }

  return(qtl)
}
