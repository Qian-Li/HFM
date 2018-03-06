#include "HFM.h"

/***********************************************************************************/
/* STRUCTURES & GLOBAL VARIABLES ***************************************************/
/***********************************************************************************/
int np, nt, nr, nbase;
//
cube low1;         // Lower bands
cube up1;          // Upper bands

/***********************************************************************************/
/* bands()                                                                         */
/***********************************************************************************/
// [[Rcpp::export]]
List bands_cpp(arma::cube const& mbeta,
               arma::cube const& mbeta2,
               arma::mat const& Bs,
               arma::mat const &beta,
               double &alpha)
{

  int  i, r, j, p, k, nsim;
  cube betasd, maxf, betaj;
  mat  mb1r, mb2r, betar, normb, Malpha, bandrj;
  vec  maxrj, qmax;

  // Get initialization summaries ---------------------------------------------
  np    = mbeta.n_slices;     // number of varying coefficients
  nt    = Bs.n_rows;          // number of evaluation points
  nr    = mbeta.n_cols;       // number of regions
  nbase = Bs.n_cols;          // number of basis functions
  nsim  = beta.n_rows;        // number of mcmc samples

  // Memory allocation --------------------------------------------------------
  betaj  = randu<cube>(nbase, np, nr);
  betasd = randu<cube>(nt, np, nr);
  low1   = randu<cube>(nt, np, nr);
  up1    = randu<cube>(nt, np, nr);
  maxf   = randu<cube>(np, nr, nsim);
  qmax   = randu<vec>(nsim);
  bandrj = randu<mat>(nt, np);
  Malpha = randu<mat>(np, nr);

  // Compute Posterior SD -----------------------------------------------------
  for(r=0; r<nr; r++){
    // mb1r = mbeta.slice(r);
    mb1r = longslice(mbeta, r); //nt*np
    // mb2r = mbeta2.slice(r);
    mb2r = longslice(mbeta2, r);
    betasd.slice(r) = sqrt(mb2r - mb1r % mb1r);
  }

  // Read posterior by row ----------------------------------------------------
  for(i=0; i<nsim; i++){
    k = 0;
  // Loop over MCMC samples ---------------------------------------------------
    for(p=0; p<np; p++){
      for(r=0; r<nr; r++){
        for(j=0; j<nbase; j++){
          betaj(j,p,r) = beta(i,k);    // Read beta sample into cube
          k++;
        }
      }
    }
  // Compute simoultaneous posterior bands -----------------------------------
   for(r=0; r<nr; r++){
     betar = betaj.slice(r);      // np x nbase
     normb = abs((Bs * betar - longslice(mbeta, r))/betasd.slice(r));  // nt x np
     for(j=0; j<np; j++){
        maxf(j,r,i) = max(normb.col(j));
     }
    }
  }
  for(r=0; r<nr; r++){
    for(j=0; j<np; j++){
      maxrj       = maxf.tube(j,r);
      qmax        = sort(maxrj);
      Malpha(j,r) = qmax(nsim*(1.0 - alpha - 0.00000001));
    }
  }
  // Form Simultaneous bands --------------------------------------------------
  for(r=0; r<nr; r++){
      // mb1r = mbeta.slice(r);  // np x nt
      mb1r = longslice(mbeta, r);
      mb2r = betasd.slice(r); // nt * np
      for(j=0; j<np; j++){
        bandrj.col(j) = mb1r.col(j) + Malpha(j,r)*mb2r.col(j);
      }
      up1.slice(r)  = bandrj;
      for(j=0; j<np; j++){
        bandrj.col(j) = mb1r.col(j) - Malpha(j,r)*mb2r.col(j);
      }
      low1.slice(r) = bandrj;
  }

  // Output to R --------------------------------------------------------------
  return List::create(
    Named("low")  = low1,
    Named("mean") = mbeta,
    Named("up")   = up1);
}

// END bands ------------------------------------------------------------------