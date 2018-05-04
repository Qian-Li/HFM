#include "HFM.h"

/***********************************************************************************/
/* STRUCTURES & GLOBAL VARIABLES ***************************************************/
/***********************************************************************************/
int np, nt, nr, nc, nbase;
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
               arma::mat const &XX,
               double &alpha)
{

  int  i, r, j, p, k, nsim;
  cube betasd, maxf, betaj;
  mat  mb1r, mb2r, betar, normb, Malpha, bandrj;
  vec  maxrj, qmax;

  // Get initialization summaries ---------------------------------------------
  np    = mbeta.n_slices;     // number of fixed effects parameters
  nt    = Bs.n_rows;          // number of time points
  nr    = mbeta.n_cols;       // number of regions
  nbase = Bs.n_cols;          // number of basis functions
  nsim  = beta.n_rows;        // number of mcmc samples
  nc    = XX.n_rows;          // number of contrasts

  // Memory allocation --------------------------------------------------------
  betaj  = randu<cube>(nbase, np, nr);
  betasd = randu<cube>(nt, np, nr);
  low1   = randu<cube>(nt, nc, nr);
  up1    = randu<cube>(nt, nc, nr);
  maxf   = randu<cube>(nc, nr, nsim);
  qmax   = randu<vec>(nsim);
  bandrj = randu<mat>(nt, nc);
  Malpha = randu<mat>(nc, nr);

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
     normb = abs((Bs * betar - longslice(mbeta, r))/betasd.slice(r) * XX.t());  // nt x nc
     for(j=0; j<nc; j++){
        maxf(j,r,i) = max(normb.col(j));
     }
    }
  }
  for(r=0; r<nr; r++){
    for(j=0; j<nc; j++){
      maxrj       = maxf.tube(j,r);
      qmax        = sort(maxrj);
      Malpha(j,r) = qmax(nsim*(1.0 - alpha - 0.00000001));
    }
  }
  // Form Simultaneous bands --------------------------------------------------
  for(r=0; r<nr; r++){
      // mb1r = mbeta.slice(r);  // np x nt
      mb1r = longslice(mbeta, r) * XX.t();
      mb2r = sqrt(betasd.slice(r) % betasd.slice(r) * (XX.t() % XX.t()));  // nt * nc
      for(j=0; j<nc; j++){
        bandrj.col(j) = mb1r.col(j) + Malpha(j,r)*mb2r.col(j);
      }
      up1.slice(r)  = bandrj;
      for(j=0; j<nc; j++){
        bandrj.col(j) = mb1r.col(j) - Malpha(j,r)*mb2r.col(j);
      }
      low1.slice(r) = bandrj;
  }
  cube out = cubexmat(mbeta,XX);

  // Output to R --------------------------------------------------------------
  return List::create(
    Named("low")  = low1,
    Named("mean") = out,
    Named("up")   = up1);
}

/***********************************************************************************/
/* DAbands()                                                                       */
/***********************************************************************************/
// [[Rcpp::export]]
List da_bands_cpp(arma::cube const& mbeta,
                  arma::cube const& mbeta2,
                  arma::mat const& Bs,
                  arma::mat const &beta,
                  arma::mat const &XX,
                  double &alpha)
{
  
  int  i, r, j, p, k, nsim;
  cube betasd, maxf, betaj,mmbeta;
  mat  mb1r, mb2r, betar, normb, Malpha, bandrj, bandrg;
  vec  maxrj, qmax;
  
  // Get initialization summaries ---------------------------------------------
  np    = mbeta.n_slices;     // number of fixed effects parameters
  nt    = Bs.n_rows;          // number of time points
  nr    = mbeta.n_cols;       // number of regions
  nbase = Bs.n_cols;          // number of basis functions
  nsim  = beta.n_rows;        // number of mcmc samples
  nc    = XX.n_rows;          // number of contrasts
  int nn=4;                   // number of region pairs, Left-Right
  
  // Memory allocation --------------------------------------------------------
  betaj  = randu<cube>(nbase, np, nr);
  betasd = randu<cube>(nt, np, nr);
  low1   = randu<cube>(nt, nc, nn);
  up1    = randu<cube>(nt, nc, nn);
  maxf   = randu<cube>(nc, nn, nsim);
  qmax   = randu<vec>(nsim);
  bandrj = randu<mat>(nt, nc);
  bandrg = randu<mat>(nt, nc);
  Malpha = randu<mat>(nc, nr);
  
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
    // pair1
    normb = (Bs*betaj.slice(0) - longslice(mbeta, 0))/betasd.slice(0)*XX.t() - 
      (Bs*betaj.slice(1) - longslice(mbeta, 1))/betasd.slice(1)*XX.t();
    for(j=0; j<nc; j++){
      maxf(j,0,i) = max(normb.col(j));
    }
    // pair2
    normb = (Bs*betaj.slice(2) - longslice(mbeta, 2))/betasd.slice(2)*XX.t() - 
      (Bs*betaj.slice(3) - longslice(mbeta, 3))/betasd.slice(3)*XX.t();
    for(j=0; j<nc; j++){
      maxf(j,1,i) = max(normb.col(j));
    }
    // pair3
    normb = (Bs*betaj.slice(5) - longslice(mbeta, 5))/betasd.slice(5)*XX.t() - 
      (Bs*betaj.slice(7) - longslice(mbeta, 7))/betasd.slice(7)*XX.t();
    for(j=0; j<nc; j++){
      maxf(j,2,i) = max(normb.col(j));
    }
    // pair4
    normb = (Bs*betaj.slice(8) - longslice(mbeta, 8))/betasd.slice(8)*XX.t() - 
      (Bs*betaj.slice(10) - longslice(mbeta, 10))/betasd.slice(10)*XX.t();
    for(j=0; j<nc; j++){
      maxf(j,3,i) = max(normb.col(j));
    }
  //
  }
  for(r=0; r<nn; r++){
    for(j=0; j<nc; j++){
      maxrj       = maxf.tube(j,r);
      qmax        = sort(maxrj);
      Malpha(j,r) = qmax(nsim*(1.0 - alpha - 0.00000001));
    }
  }
  // Form Simultaneous bands --------------------------------------------------
  mat mb1rr, mb2rr;
  //pair 1
  mb1rr = (longslice(mbeta,0) - longslice(mbeta,1)) * XX.t();
  mb2rr = sqrt((betasd.slice(0) % betasd.slice(0) + betasd.slice(1) % betasd.slice(1)) * (XX.t() % XX.t()));
  for(j=0; j<nc; j++){
    bandrj.col(j) = mb1rr.col(j) + Malpha(j,0)*mb2rr.col(j);
    bandrg.col(j) = mb1rr.col(j) - Malpha(j,0)*mb2rr.col(j);
  }
  up1.slice(0) = bandrj; low1.slice(0) = bandrg;
  //
  //pair 2
  mb1rr = (longslice(mbeta,2) - longslice(mbeta,3)) * XX.t();
  mb2rr = sqrt((betasd.slice(2) % betasd.slice(2) + betasd.slice(3) % betasd.slice(3)) * (XX.t() % XX.t()));
  for(j=0; j<nc; j++){
    bandrj.col(j) = mb1rr.col(j) + Malpha(j,1)*mb2rr.col(j);
    bandrg.col(j) = mb1rr.col(j) - Malpha(j,1)*mb2rr.col(j);
  }
  up1.slice(1) = bandrj; low1.slice(1) = bandrg;
  //
  //pair 3
  mb1rr = (longslice(mbeta,5) - longslice(mbeta,7)) * XX.t();
  mb2rr = sqrt((betasd.slice(7) % betasd.slice(7) + betasd.slice(5) % betasd.slice(5)) * (XX.t() % XX.t()));
  for(j=0; j<nc; j++){
    bandrj.col(j) = mb1rr.col(j) + Malpha(j,2)*mb2rr.col(j);
    bandrg.col(j) = mb1rr.col(j) - Malpha(j,2)*mb2rr.col(j);
  }
  up1.slice(2) = bandrj; low1.slice(2) = bandrg;
  //
  //pair 4
  mb1rr = (longslice(mbeta,8) - longslice(mbeta,10)) * XX.t();
  mb2rr = sqrt((betasd.slice(8) % betasd.slice(8) + betasd.slice(10) % betasd.slice(10)) * (XX.t() % XX.t()));
  for(j=0; j<nc; j++){
    bandrj.col(j) = mb1rr.col(j) + Malpha(j,3)*mb2rr.col(j);
    bandrg.col(j) = mb1rr.col(j) - Malpha(j,3)*mb2rr.col(j);
  }
  up1.slice(3) = bandrj; low1.slice(3) = bandrg;
  //
  cube out = cubexmat(mbeta,XX);
  
  // Output to R --------------------------------------------------------------
  return List::create(
    Named("low")  = low1,
    Named("mean") = out,
    Named("up")   = up1);
}

// END bands ------------------------------------------------------------------