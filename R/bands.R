#' Simultaneous posterior bands from hierarchical regionally referenced varying coefficient model
#'
#' This function is a wrapper for a Cpp Gibbs Sampler. Region-specific regression is modeled through a varying coefficient functions.
#'
#' @param fit An object returned from \code{hfm}
#' @param beta A matrix of mcmc samples from \code{hfm}
#' @param alpha Value between 0 and 1 for the (1-alpha)\% Credible Band.
#' @return
#' An array with simultaneous bands for the varying coefficients in fit
#'
#' @details
#' For region r, subject i, let y_i(t) ...
#'
#' @examples
#'
#'
#' @export
bands <- function(fit, beta, alpha=NULL){
  #
  ## Construct Bands in C -----------------------------------------------------
  #
  if(is.null(alpha)) alpha = 0.05;
  #
  band1  <- bands_cpp(fit$coef, fit$coef2, fit$Bs, beta, alpha)
  #
  return(list(mean=band1$mean, low=band1$low, up=band1$up))
}

#' @useDynLib HFM
#' @importFrom Rcpp sourceCpp
NULL
