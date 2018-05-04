#' Simultaneous posterior bands for HFM model fits
#'
#' This function is a wrapper for a Cpp Gibbs Sampler. Region-specific regression is modeled through a varying coefficient functions.
#'
#' @param fit An object returned from \code{hfm}
#' @param beta A matrix of mcmc samples from \code{hfm}
#' @param contrast A matrix of contrasts of interests
#' @param alpha Value between 0 and 1 for the (1-alpha)\% Credible Band
#' @param LR Logical, Left-Right difference (specific to DA)
#' @return
#' An array with simultaneous bands for the varying coefficients in fit
#' \itemize{
#' \item{mean:}{ Posterior mean}
#' \item{low:}{ Lower bound of the simultaneous confidence band}
#' \item{up:}{ Upper bound of the simultaneous confidence band}
#' }
#'
#' @details
#' For details please point to the references in our manuscript...
#'
#'
#' @export
bands <- function(fit, beta, contrast=NULL, alpha=NULL, LR = FALSE){
  #
  ## Construct Bands in C -----------------------------------------------------
  #
  if(is.null(alpha)) alpha = 0.05;
  if(is.null(contrast)) contrast = diag(dim(fit$coef)[3]);
  #
  if(!LR){
    band1  <- bands_cpp(fit$coef, fit$coef2, fit$Bs, beta, contrast, alpha)
  } else{
    band1  <- da_bands_cpp(fit$coef, fit$coef2, fit$Bs, beta, contrast, alpha)
  }
  #
  return(list(mean=band1$mean, low=band1$low, up=band1$up))
}

#' @useDynLib HFM
#' @importFrom Rcpp sourceCpp
NULL
