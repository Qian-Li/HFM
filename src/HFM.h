#ifndef __HFM_H__
#define __HFM_H__

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <fstream>
//
using namespace std;
using namespace arma;
using namespace Rcpp;
//
//
/************************************************************************************/
/* Data Structure                                                                   */
/************************************************************************************/
//
struct dta{
  icube miss; // Indicator cube for missing observations
  int ns;     // Number of subjects
  int nt;     // Number of time segments
  int np;     // Number of covariates
  int nr;     // Number of brain regions
  vec nmiss;  // Number of missing observations
};
//
/************************************************************************************/
/* Parameters Structure                                                             */
/************************************************************************************/
//
struct pars{
  // Common Model Parameters---------------------------------------------------------
  cube bi;         // Individual coefficient (fixed + random)
  cube mbi;        // Individual mean coefficient (fixed)
  cube beta;       // Coefficients mean;
  mat He;          // Region-specific Residual precision, diagonal
  mat taup;        // coef residual precision, diagonal
  mat tauq;        // coef residual precision, diagonal
  // Static Summaries --------------------------------------------------------------
  int nbase;       // Number of basis functions;
  cube fit;        // Current sub-level fit
  mat BtB;         // Bs'Bs;
  mat iBtB;        // pinv(Bs'Bs);
  mat XtXi;        // pinv(X'X)
  mat X;           // design matrix X
};
//
struct parNB{
  // Naive Bayesian Implementation Parameters-------------------------------------
  mat iS;         // Spatial Precision
  double lbeta;   // Population smoothing
  double lbi;     // Subject-level smoothing
  // Static Summaries -------------------------------------------------------------
  mat Df;         // Smoothing Matrix
  mat cholXX;     // Cholesky of X'X
  mat cholDf;     // Cholesky of Df
};
//
struct parNS{
  // Non-Separable Implementation Parameters--------------------------------------
  cube Ld;         // Latent Factor loading
  mat H;           // Latent Factors
  cube Omega;      // Spatial covariance, free
  vec pen;         // Column covariance penalty on Loadings
  mat HtH;         // H'H
  // static
  int nfac;        // Number of factors;
};
//
struct parSS{
  // Strongly Separable Implementation Parameters---------------------------------
  mat LL;          // Latent Factor loading, Left
  mat LR;          // Latent Factor loading, Right
  cube H;          // Latent Factors
  mat tauL;        // Precision, Left
  mat tauR;        // Precision, Right
  vec penL;        // Penalty, Left
  vec penR;        // Penalty, Right
  // Static Summaries --------------------------------------------------------------
  int nfacL;       // Number of factors;
  int nfacR;       // Number of factors;  
};
/************************************************************************************/
/* Prior Structure                                                                  */
/************************************************************************************/
//
struct prNB{
  // Priors for Naive Bayesian Priors;
  // Error variance
  double ae;
  double be;
  // subject level smoothing
  double ab;
  double bb;
  // spatial smoothing
  double aS;
  mat W;
};
//
struct prNS{
  // Prior for Non-Separable Priors;
  // Error variance
  double ae;
  double be;
  // Coefficients variances (Sigp and Sigq)
  double athe;
  double bthe;
  // Penalties shapes, the first and the rest
  double ap1;
  double ap2;
  // spatial covs (IW)
  double aOm;
  mat Om0;
  // Latent factors;
  vec eta0;
};
//
struct prSS{
  // Prior for Strongly Separable Priors
  // Error variance
  double ae;
  double be;
  // Coefficients variances (Sigp and Sigq)
  double athe;
  double bthe;
  // Penalties shapes, the first and the rest
  double ap1;
  double ap2;
  // Loading precision
  double nuL;
  double nuR;
  // Latent factors;
  vec eta0;
};
/************************************************************************************/
/* Estimate Structure                                                               */
/************************************************************************************/
//
struct ergodics{
   cube fit;    // Posterior mean fit (bi*Bs')
   cube beta;   // Posterior mean population regression functions
   cube beta2;  // Second posterior moment of population regression functions
   double n;    // number of samples used in mc estimators
   // More can be added here;
   mat covcoef; // Covariance of the coefficient matrices
   mat LF;      // Latent factors; zeros if not exists;
   mat LF2;     // latent factors; second moments;
};
/************************************************************************************/
/* Utility Functions                                                                */
/************************************************************************************/
//
// Bayesian Regression and Triangular Systems ------------------------------------
vec BayesRegL(vec const &b, mat const &L);
void trisolve(vec &solution, mat const &L, vec const &x, bool const &lower);
//
// Multivariate Samplers
mat rWish(mat const &S, double v);
mat rMN(mat const &m, mat const &S, mat const& V);
// Other utilities
mat cov2cor(mat S);
double meanNA(vec &y);
mat vec2mat(vec &v, int const& nr, int const& nc);
void vec2long(cube &a, vec const &v, int const& j);
cube cubexmat(cube &c, mat const&m);
cube matxcube(mat const& L, cube const&C, mat const&R);
mat cube2mat(cube const &c);
cube mat2cube(mat const& m, int const& nr, int const& nc);
//
// slicing of Armadillo cubes ---------------------------
mat longslice(cube const &A, int r);
mat flatslice(cube const &A, int r);
//
#endif