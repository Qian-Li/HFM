#' Multivariate Hierarchical Functional Models: Latent Factors Implementation (SS and NS prior)
#'
#' This function is a wrapper for a Cpp Gibbs Sampler. Multivariate functional observations, for examples t time points observed over p regions simultaneously, are projected onto B-spline bases, and further modeled by means of a latent factors encoded hierarchical priors. We implement both a Strongly Separable and Non-Separable priors.
#'
#' @param y A matrix with vectorized region measurements in columns.
#' @param t A vector with vectorized observed time points or segments.
#' @param X A regression design matrix (subject by predictor).
#' @param nsegs An integer vector with the number of observed segments per subject.
#' @param nLF An integer (Non-separable Prior) or vector of two integers (Separable Prior)
#' @param knots Number of spline knots or a vector of spline knots
#' @param mcmc A list of mcmc specifications. Default is \code{burnin = 1000, nsim = 2000, thin = 1}
#' @param display Logical, whether to display the sampling progress bar. Default is \code{true}
#' @return
#' A list of with model fit summaries list(fit, coef, coef2, cov, facs, fac_vars, y, t, orig_t, Bs, nsim).
#' \itemize{
#' \item{fit:}{ Posterior mean fit to y in a 3d array [subject, region, segments]}
#' \item{coef:}{ Array containing posterior mean of regressional coefficients [segments, region, covariates]}
#' \item{coef2:}{ Array containing posterior mean of the second moments of regressional coefficients [segments, region, covariates]}
#' \item{cov:}{ Matrix of the posterior mean of cov(Theta), size of \code{nBS*nregion} square}
#' \item{facs:}{ Posterior mean of latent factors per individual}
#' \item{fac_vars:}{Posterior vairance of latent factors per individual}
#' \item{y:}{ Reshaped data array (NA filled) [subject, region, segments]}
#' \item{t:}{ Reshaped segment array (NA filled) [subject, segments]}
#' \item{orig_t:}{ Time as input}
#' \item{Bs:}{ Matrix of spline polynomials used in the model [segments, nBS]}
#' \item{nsim:}{ Number of MCMC samples}
#' }
#'
#' @details
#' For region r, subject i, let y_i(t) ... (Please refer to our manuscript...)
#'
#'
#'
#' @export
hfm.LF <- function(y, t, X, nsegs, nLF = NULL, knots=NULL, mcmc=NULL, display=TRUE){
  #
  # Normalize longitudinal time scale to be between 0 and 1 -------------------
  #
  tmin <- min(t, na.rm=TRUE)
  tmax <- max(t, na.rm=TRUE)
  tn   <- (t - tmin)/(tmax - tmin)
  if(!is.null(knots)) knots <- (knots - tmin)/(tmax - tmin)
  #
  ## Get basic data summaries -------------------------------------------------
  #
  nreg <- ncol(y)            # number of regions
  seg1 <- sort(unique(tn))    # unique segments or time-points
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
  # # To be edited out -------------
  # #
  # # Remove margins if too many NA's
  # yy <- y3[,1,]
  # nn <- rep(0, ncol(yy))
  # for(j in 1:ncol(yy)){
  #   yj    <- yy[,j]
  #   nn[j] <- length(yj[is.na(yj)])/ncol(yy)
  # }
  # # Leave aside for now.
  # # Trim segments with too many NA's
  # # y3   <- y3[ , , nn < 0.20]
  # # x1   <- x1[ , nn < 0.20]
  # # seg1 <- seg1[nn < 0.20]
  # # ns1  <- length(seg1)       # number of unique segments
  # #
  # #
  y3_b <- y3
  y3[is.na(y3)] <- 12345.12345
  
  ## Get B-Spline Matrix ------------------------------------------------------
  if(is.null(knots)){
    ns1 <- round(length(seg1)/5)
    knots1 <- seq(0,seg1[ns1],length=(ns1+2))[2:(ns1+1)]
  } else if(length(knots)==1){
    knots1 = seq(0,seg1[ns1],length=(knots+2))[2:(knots+1)]
  } else if(length(knots)>1) {
    knots1 <- (knots - min(seg1))/(max(seg1) - min(seg1))  ## renormalize knots
  } 
  BSX    <- splines::bs(seg1, knots=knots1, intercept=T)
  nbase  <- ncol(BSX)
  #
  ## Get MCMC summaries -------------------------------------------------------
  if(is.null(mcmc)){
    mcmc$burnin <- 1000
    mcmc$nsim   <- 2000
    mcmc$thin   <- 1
  }
  if(length(nLF) == 2){
    # Create Output directory ---------------------------------------------------
    files <- list.files();
    dir   <- files == "post_SS";
    if(length(files[dir]) == 0) dir.create("post_SS");
    nLF = ceiling(nLF)
    ## sandwich HFMM- MCMC---
    cppfit = SSmix_mcmc(y3, X, BSX, nLF[1], nLF[2], mcmc$burnin, mcmc$nsim, mcmc$thin, display)
  } else if(length(nLF == 1)){
    # Create Output directory ---------------------------------------------------
    files <- list.files();
    dir   <- files == "post_NS";
    if(length(files[dir]) == 0) dir.create("post_NS");
    nLF = ceiling(nLF)
    ## vectorized HFMM- MCMC---
    cppfit <- NSmix_mcmc(y3, X, BSX, nLF,mcmc$burnin, mcmc$nsim, mcmc$thin, display)
  } else {
    # Create Output directory ---------------------------------------------------
    files <- list.files();
    dir   <- files == "post_NS";
    if(length(files[dir]) == 0) dir.create("post_NS");
    nLF = ceiling(nreg * nbase / 2.0);
    warning('Number of Latent Factors not given, default nLF = ',nLF)
    ## vectorized HFMM- MCMC---
    cppfit <- NSmix_mcmc(y3, X, BSX, nLF,mcmc$burnin, mcmc$nsim, mcmc$thin, display)
  }
  ## END MCMC -----------------------------------------------------------------
  y3[is.na(y3_b)] <- NA;
  return(list(fit=cppfit$fit, coef=cppfit$coef, coef2=cppfit$coef2, 
              covs=cppfit$covcoef, facs = cppfit$LF, fac_vars = cppfit$LF2, 
              y=y3, t=x1, orig_t = xx, Bs=BSX, nsim=mcmc$nsim))
}
#' @useDynLib HFM
#' @importFrom Rcpp sourceCpp
NULL
