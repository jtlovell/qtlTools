% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scanone.GWERk.R
\name{scanone.GWERk}
\alias{scanone.GWERk}
\title{Scanone qtl analysis permitting GWERk analysis.}
\usage{
scanone.GWERk(cross, chr, pheno.col = 1, model = c("normal", "binary",
  "2part", "np"), method = c("em", "imp", "hk", "ehk", "mr", "mr-imp",
  "mr-argmax"), addcovar = NULL, intcovar = NULL, weights = NULL,
  use = c("all.obs", "complete.obs"), upper = FALSE,
  ties.random = FALSE, start = NULL, maxit = 4000, tol = 1e-04,
  n.perm, perm.Xsp = FALSE, perm.strata = NULL, verbose,
  batchsize = 250, n.cluster = 1, ind.noqtl, GWERk = 0)
}
\description{
\code{scanone.GWERk}  see qtl::scanone for details. GWERk = 0 runs normal scanone
permutations. Otherwise GWERk perms are used.
}
