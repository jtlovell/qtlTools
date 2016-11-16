#' @title Conduct mixed effect QTL mapping.
#'
#' @description
#' \code{scanone.lmer} Use lmer linear mixed effect models
#' as the engine for scanone
#'
#' @param cross A cross object
#' @param base.formula The formula of a model that serves as the baseline. Must
#' contain at least one random effect per lmer's functionality. All terms in the
#' formula must be contained in the `cross` phenotype matrix slot.
#' @param test.formula The formula containing an additional fixed effect to test.
#' The word "QTL" must be contained here, and may also occur in the base.formula.
#' This term is a placeholder and is replaced by the genotypes at each postion in the
#' conditional probability matrix when lmer is called.
#' @param ... Additional arguments to pass on to lmer.
#'
#' @details This function is still under development and is too
#' slow for permutations as it stands now. An overly conservative approach is to
#' adjust the p-values. However, since each marker is not independent, this results in a
#' severe reduction in power.
#' @return An object of class scanone
#' @examples
#' \dontrun{
#' data(fake.f2)
#' cross<-fake.f2
#' cross<-subset(cross, ind = c(1:nind(cross),1:nind(cross)))
#' cross$pheno<-data.frame(cross$pheno, block = rep(1:(nind(cross)/2),2))
#' cross$pheno<-data.frame(cross$pheno, pheno2 = jitter(cross$pheno$phenotype))
#' sex<-data.frame(sex = pull.pheno(cross, pheno.col = "sex"))
#' cross = calc.genoprob(cross)
#'
#' base.formula = "pheno2 ~ sex + (1|block)"
#' qtl.formula = "pheno2 ~ sex + QTL + sex*QTL + (1|block)"
#'
#' s1<-scanone.lmer(cross, base.formula = base.formula, test.formula = qtl.formula)
#' par(mfrow = c(2,1))
#' plot(s1, main = "linear mixed effect model")
#' plot(scanone(cross, pheno.col = "pheno2", intcovar = sex), col = "grey", main = "scanone")
#' }
#' @export
scanone.lmer<-function(cross, base.formula,
                       test.formula, ...){
  if(!"id" %in% tolower(phenames(cross))){
    cross$pheno<-data.frame(id = 1:nind(cross), cross$pheno)
  }
  if("prob" %in% names(cross$geno[[1]])){
    atr<-attributes(cross$geno[[1]]$prob)
  }else{
    if("draws" %in% names(cross$geno[[1]])){
      atr<-attributes(cross$geno[[1]]$draws)
    }else{
      stop("run either calc.genoprob or sim.geno first.\n")
    }
  }
  ag<-argmax.geno(cross, step = atr$step, error.prob = atr$error.prob,
                     off.end = atr$off.end, map.function = atr$map.function,
                     stepwidth = atr$stepwidth)
  s1<-scanone(cross,method = ifelse("prob" %in% names(cross$geno[[1]]),"hk","imp"))[,-3]

  gen<-pull.argmaxgeno(ag, include.pos.info=F)
  if(class(cross)[1] == "4way"){
    for(i in colnames(gen)) gen[,i]<-as.character(gen[,i])
  }

  mars<-rownames(s1)
  dat<-cbind(pull.pheno(cross), gen)

  s1$negLog10P<-sapply(mars, function(x){
    form.base <- as.formula(gsub("QTL",x,base.formula))
    no.marker <- lmer(form.base, data = dat, REML = FALSE, ...)
    form.full<-as.formula(gsub("QTL",x,test.formula))
    w.marker <- lmer(form.full, data = dat, REML = FALSE, ...)
    return(-log10(anova(no.marker, w.marker)[["Pr(>Chisq)"]][2]))
  })
  return(s1)
}
