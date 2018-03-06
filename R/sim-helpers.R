## Simulation Helpers

#' Simulation Helper: quantile knots
#'
#' Getting default Knots
#'
#' @export
default.knots <- function(x,num.knots)
{
  # Delete repeated values from x

  x <- unique(x)

  # Work out the default number of knots

  if (missing(num.knots))
  {
    n <- length(x)
    d <- max(4,floor(n/35))
    num.knots <- floor(n/d - 1)
  }

  nx <- names(x)
  x <- as.vector(x)
  nax <- is.na(x)
  if(nas <- any(nax))
    x <- x[!nax]

  knots <- seq(0,1,length=num.knots+2)[-c(1,num.knots+2)]
  knots <- quantile(x,knots)

  names(knots) <- NULL

  return(knots)
}

## Gaussian Process Kernal:

#' Simulation Helper: GP-kernel
#'
#' Calculates the Gaussian Process kernel for simulation
#' @export
calcSigma<-function(X1,X2,l=1){
  # Sigma<-matrix(rep(0,length(X1)*length(X2)),nrow=length(X1))
  # for(i in 1:nrow(Sigma)){
  #   for (j in 1:ncol(Sigma)) Sigma[i,j]<-exp(-1/2*(abs(X1[i]-X2[j])/l)^2)
  # }
  sig <- outer(X1, X2, "-"); Sigma <- exp(-1/2 * (abs(sig)/l)^2)
  return(Sigma)
}
# Importation of ourside functions
#' @importFrom stats quantile
#' @importFrom splines bs
NULL