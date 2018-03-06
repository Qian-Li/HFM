#' Model Generated Simulations
#' 
#' This simulation function generates functional observations from the HFMM models (separable and non-separable).
#' 
#' @export
#' @importFrom MASS mvrnorm
sim.HFM <- function(nsub = 10, nreg = 6, tt = NULL, type = "sep", pr.miss = 0.2, SNR = 1, 
                    nLF = 5, nbase = 10, sp_rho = runif(1)*0.9
                    ){
  ## Simulated time
  if(is.null(tt)){tt <- seq(0,30, by = .5)}
  if(pr.miss >= 1 | pr.miss < 0){stop("Missing percentage must be in (0,1)!")}
  p <- rep(pr.miss, nsub * length(tt)); p <- (runif(length(p)) >= p)
  p <- split(p, rep(1:nsub, rep(length(tt),nsub)))
  time  <- unlist(lapply(p, function(x) tt[x]))
  nsegs <- unname(unlist(lapply(p,sum)))
  ## Spline project matrix
  Bs    <- splines::bs(tt, df = nbase, intercept = TRUE)
  ## 
  sp_cov <- sp_rho^abs(outer(1:nreg,1:nreg,"-")); 
  signal <- list();
  time <- unname(unlist(lapply(p, function(x) tt[x])))
  if(type == "nonsep"){
    ## Generative from vectorized model
    pen = rgamma(nLF, 2, 1); tau = cumprod(pen);
    loading <- c()
    for(base in 1:nbase){
      nu = nreg + 10;
      chol_sp <- chol(solve(rWishart(1,nu,solve(sp_cov)/nu)[,,1]))
      loading <- rbind(loading, 
                       t(chol_sp) %*% MASS::mvrnorm(nreg, rep(0,nLF), diag(1/tau)))
    }
    complete_data <- list();
    bs_err = 1/kronecker(rgamma(nbase,10,10), rgamma(nreg,10,10));
    true_cov = diag(bs_err) + loading %*% t(loading);
    true_fac = c()
    for(i in 1:nsub){
      factors = rnorm(nLF)
      true_fac= rbind(true_fac, factors)
      coefs = matrix(loading %*% factors, nrow = nreg, ncol = nbase)
      coefs = coefs + matrix(rnorm(nreg*nbase), nrow = nreg, ncol = nbase) *
        matrix(sqrt(bs_err), nrow = nreg, ncol = nbase);
      sig = mean(diag(var(Bs %*% t(coefs))))
      sim_vec <- Bs %*% t(coefs) + t(MASS::mvrnorm(nreg, rep(0, length(tt)),
                                             Sigma = sig/SNR*diag(length(tt))))
      complete_data[[i]] <- sim_vec[p[[i]],]
      signal[[i]] <- Bs %*% t(coefs);
    }
    results = list(data = do.call(rbind, lapply(complete_data, as.matrix)),
                   sig  = signal,
                   time = time,
                   nsegs=nsegs,
                   snr= SNR,
                   nbase = nbase,
                   nLF = nLF,
                   true_cov = true_cov,
                   factors = true_fac)
  } else {
    # ## Individual curves smoothness ++ for smoother
    ind.smooth = 2
    # ## Gaussian Process covs
    knots <- default.knots(tt,num.knots = nbase )
    bs_cov = calcSigma(knots, knots, l = ind.smooth)
    # bs_cov = bs_rho^(abs(outer(1:nbase,1:nbase,"-")))
    true_cov = kronecker(bs_cov, sp_cov)
    ## Regional AR-1 dependency
    chol_sp <- chol(sp_cov);
    ## Sample data from GP
    complete_data <- list()
    for(sub in 1:nsub){
      bij <- MASS::mvrnorm(nreg, rep(0, nbase), Sigma = bs_cov) %*% t(Bs) ## true signal for individuals
      sim_vec <- t(bij) %*% chol_sp + 
        t(MASS::mvrnorm(nreg, rep(0, length(tt)), Sigma = 1.0/SNR*diag(length(tt))))
      complete_data[[sub]] <- sim_vec[p[[sub]],]
      signal[[sub]] <- t(bij) %*% chol_sp;
    }
    results = list(data = do.call(rbind, lapply(complete_data, as.matrix)),
                   sig  = signal,
                   time = time,
                   nsegs=nsegs,
                   snr= SNR,
                   nbase = nbase,
                   nLF = nLF,
                   true_cov = true_cov)
  }
  return(results)
}
#' @importFrom splines bs
#' @importFrom stats rWishart rgamma rnorm runif var
NULL