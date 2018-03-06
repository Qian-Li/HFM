#' Multivariate Hierarchical Functional Models: Matrix Normal Implementation (NB prior)
#'
#' This function is a wrapper for a Cpp Gibbs Sampler. Multivariate functional observations, for examples t time points observed over p regions simultaneously, are projected onto B-spline bases, and further modeled by means of a matrix normal prior.
#'
#' @param y A matrix with vectorized region measurements in columns.
#' @param t A vector with vectorized observed time points or segments.
#' @param X A regression design matrix (subject by predictor).
#' @param nsegs An integer vector with the number of observed segments per subject.
#' @param knots Number of spline knots or a vector of spline knots
#' @param mcmc A list of mcmc specifications. Default is (burnin = 1000, nsim = 2000, thin = 1)
#' @return
#' A list of with model fit summaries list(fit, coef, y, t, Bs).
#' \itemize{
#' \item{fit:}{ Posterior mean fit to y in a 3d array [subject, region, segmenr]}
#' \item{coef:}{ Array containing mean posterior regression coefficients [time, predictor, region]}
#' \item{coef2:}{ Array containing second moment of posterior regression coefficients [time, predictor, region]}
#' \item{y:}{ Reashaped data array [subject, region, segment]}
#' \item{t:}{ Matrix of observed segments [subject, segments]}
#' \item{Bs:}{ Matrix of spline polynomials used in the model}
#' \item{X:}{ Covariate matrix}
#' \item{nsim:}{ Number of MCMC samples}
#' }
#'
#' @details
#' For region r, subject i, let y_i(t) ...
#'
#' @examples
#'
#'
#' @export
hfm.NB <- function(y, t, X, nsegs, knots=NULL, mcmc=NULL, spatial = FALSE){
  #
  # Normalize longitudinal time scale to be between 0 and 1 -------------------
  #
  tmin <- min(t, na.rm=TRUE)
  tmax <- max(t, na.rm=TRUE)
  tn    <- (t - tmin)/(tmax - tmin)
  if(!is.null(knots)) knots <- (knots - tmin)/(tmax - tmin)
  #
  ## Get basic data summaries -------------------------------------------------
  #
  nreg <- ncol(y)            # number of regions
  seg1 <- sort(unique(tn))   # unique segments or time-points
  seg2 <- sort(unique(t))
  ns1  <- length(seg1)       # number of unique segments
  nsub <- nrow(X)            # number of subjects
  #
  # Put Observed Segments in square Matrix and pad NA for unobserved -------------
  #
  i1 <- cumsum(c(1, nsegs))[1:nsub]
  i2 <- cumsum(c(1, nsegs))[2:(nsub+1)] - 1
  #
  x1 <- matrix(NA, nrow=nsub, ncol=ns1)
  xx <- matrix(NA, nrow=nsub, ncol=ns1)
  y3 <- array(NA, c(nsub, nreg, ns1))
  #
  for(i in 1:nsub){
    index         <- (seg1 %in% tn[i1[i]:i2[i]])
    x1[i,index]   <- seg1[index]
    xx[i,index]   <- seg2[index]
    y3[i,,index]  <- t(y[i1[i]:i2[i],])
  }
  #
  # To be edited out -------------
  #
  # # Remove margins if too many NA's
  # yy <- y3[,1,]
  # nn <- rep(0, ncol(yy))
  # for(j in 1:ncol(yy)){
  #     yj    <- yy[,j]
  #     nn[j] <- length(yj[is.na(yj)])/ncol(yy)
  # }
  # # Leave aside for now.
  # # Trim segments with too many NA's
  # # y3   <- y3[ , , nn < 0.20]
  # # x1   <- x1[ , nn < 0.20]
  # # seg1 <- seg1[nn < 0.20]
  # # ns1  <- length(seg1)       # number of unique segments
  #
  y3_b <- y3
  y3[is.na(y3)] <- 12345.12345

  ## Get B-Spline Matrix ------------------------------------------------------
  if(is.null(knots)){
    ns1 <- round(length(seg1)/5)
    knots1 <- seq(0,seg1[ns1],length=(ns1+2))[2:(ns1+1)]
  } else if(length(knots)==1){
    knots1 = seq(0,seg1[ns1],length=(knots+2))[2:(knots+1)]
  } else if(length(knots)>1) {
    knots1 <- (knots - min(seg1))/(max(seg1) - min(seg1))  ## renormalize knots in case outta bound
  } 
  BSX    <- splines::bs(seg1, knots=knots1, intercept=T)
  nbase  <- ncol(BSX)
  #
  # Create Output directory ---------------------------------------------------
  files <- list.files();
  dir   <- files == "post_NB";
  if(length(files[dir]) == 0) dir.create("post_NB");
  #
  ## Get MCMC summaries -------------------------------------------------------
  if(is.null(mcmc)){
    mcmc$burnin <- 1000
    mcmc$nsim   <- 2000
    mcmc$thin   <- 1
  }
  ## Run MCMC -----------------------------------------------------------------
  #
  fit1 <- NBmix_mcmc(y3, X, BSX, mcmc$burnin, mcmc$nsim, mcmc$thin, spatial)
  #
  ## END MCMC -----------------------------------------------------------------
  y3[is.na(y3_b)] <- NA;
  return(list(fit=fit1$fit, coef=fit1$coef, coef2=fit1$coef2,covs=fit1$covcoef, y=y3, t=x1, orig_t = xx, Bs=BSX, nsim=mcmc$nsim))
}

#' @useDynLib HFM
#' @importFrom Rcpp sourceCpp
NULL
