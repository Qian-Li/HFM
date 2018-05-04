## Simulation Helpers
##
#' Simulation Helper: quantile-based knots
#'
#' @param x A vector of observations
#' @param num.knots An integer, Number of knots.
#' @return A vector of \code{num.knots} length as quantile-based knots for \code{x}
#'
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
#' 
#' @param X1 A vector of input
#' @param X2 A vector of input
#' @param l A number of kernel window size
#' 
#' @export
calcSigma<-function(X1,X2,l=1){
  ## Simplified code
  sig <- outer(X1, X2, "-"); Sigma <- exp(-1/2 * (abs(sig)/l)^2)
  return(Sigma)
}
# Importation of ourside functions
#' @importFrom stats quantile
#' @importFrom splines bs
NULL