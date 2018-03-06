#include "HFM.h"
//
/***********************************************************************************/
/* OUTPUT FILES      ***************************************************************/
/***********************************************************************************/
#define fileOutput1 "post_NB/beta.txt"             //fixed-effects related coeffs
#define fileOutput2 "post_NB/precision.txt"
#define fileOutput3 "post_NB/spatial.txt"
#define fileOutput4 "post_NS/beta.txt"             //fixed-effects related coeffs  -V
#define fileOutput5 "post_NS/precision.txt"        //sigP,sigq and sige^-1         -V
#define fileOutput6 "post_NS/factors.txt"          //factors and penalties         -V
#define fileOutput7 "post_SS/beta.txt"             //fixed-effects related coeffs  -S
#define fileOutput8 "post_SS/precision.txt"        //sigP,sigq and sige^-1         -S
#define fileOutput9 "post_SS/factors.txt"          //factors and penalties         -S
//
// Output files ------------------------------------------------------------------
ofstream out1;
ofstream out2;
ofstream out3;
//
/***********************************************************************************/
/* STRUCTURES & GLOBAL VARIABLES ***************************************************/
/***********************************************************************************/
dta         dat;                   // Data summaries
pars        th;                    // Common Parameters
ergodics    hat;                   // Ergodic summaries
parNB       NB;                    // Naive Bayesian parameters
parNS       NS;                    // Non-Separable parameters
parSS       SS;                    // Strongly Separable parameters
prNB        pNB;                   // NB Prior;
prNS        pNS;                   // NS Prior;
prSS        pSS;                   // SS Prior;
//
vec yit(dat.nt), bir, vbi;         // Global var: nt*1; nbase*1; vec(bi) pq*1 vectors
mat Sinv;                          // Sigma^{-1} of size pq*pq;
//
// Utility constants
double w1, w2;
double g = 100.0;                  // Mean Prior constant, set large;
bool   ch_ok;                      // Whether cholesky failed;
int    i,r,j,k,l,t,s,iter,p;
//
/***********************************************************************************/
/* completeY()                                                                    */
/***********************************************************************************/
void completeY(cube &y)
{
  vec yn(1);
  //
  for(r=0; r<dat.nr; r++){
    for(i=0; i<dat.ns; i++){
      for(j=0; j<dat.nt; j++){
        if(dat.miss(i,r,j) == 1){
          yn.randn();
          y(i,r,j) = th.fit(i,r,j) + yn[0] * std::pow(th.He(r,r), -0.5);
        }
      }
    }
  }
}
//
/***********************************************************************************/
/* samplePrecisions()                                                              */
/***********************************************************************************/
void samplePrecisions(cube &y, bool const& spatial)
{ // Feb 1st update: Regional error precisions
  double    a, b, v;
  mat       resE, S;
  a = 0.5*(dat.ns*dat.nt - accu(dat.nmiss)) + pNB.ae;
  // Error precision ------------------------------------------------
  for(r=0; r<dat.nr; r++){
    resE  =   longslice(y, r) - longslice(th.fit, r);
    // b correction
    b     =   pNB.be;
    for(i=0; i<dat.ns; i++){
      for(t=0; t<dat.nt; t++){
        b += 0.5 * ((dat.miss(i,r,t)==0)? std::pow(resE(i,t),2.0) : 0.0);
      }
    }
    // Residual sampler;
    th.He(r, r) = Rf_rgamma(a, 1.0/b);
  }
  S = 0.0 * zeros<mat>(dat.nr, dat.nr);
  for(i=0; i<dat.ns; i++){
    S += th.bi.slice(i) * NB.Df * trans(th.bi.slice(i));
  }
  if(spatial){
    // Spatial -----------------------------
    v = dat.ns*th.nbase + pNB.aS;
    mat iS = pinv(S + inv(pNB.W));
    //
    NB.iS = rWish(iS, v);                                   // Free Regional Precision                                     
  } else {
    a     = pNB.ab + th.nbase * dat.ns * 0.5;
    for(r=0; r<dat.nr; r++){
      NB.iS(r,r) = Rf_rgamma(a, 1.0/(pNB.bb + S(r,r)*0.5));  // Diag Regional Precision
    }
  }
}
//
/***********************************************************************************/
/* samplebi()                                                                      */
/***********************************************************************************/
// 
void samplebi(mat const& Bs, cube &y)
{ //
  vec sbi;
  mat L, V, fit;
  // Random-effects cov:
  V     = kron(NB.iS, NB.Df) + kron(th.He, th.BtB);
  ch_ok = chol(L, V, "lower");
  if(!ch_ok){V.print("V");stop("bi sampler failed, chol() Error on V at iter " + std::to_string(iter));}
  for(i=0; i<dat.ns; i++){
    vbi = vectorise(trans(Bs) * trans(flatslice(y,i)) * th.He);
    vbi += vectorise(NB.Df * th.mbi.slice(i).t() * NB.iS);
    sbi = BayesRegL(vbi, L);
    th.bi.slice(i) = trans(vec2mat(sbi, th.nbase, dat.nr));   // Theta-i, not centered
    fit = th.bi.slice(i) * trans(Bs); //p*nt
    for(r=0; r<dat.nr; r++){
      th.fit.tube(i,r) = conv_to<vec>::from(fit.row(r));
    }
    th.bi.slice(i) = th.bi.slice(i) - th.mbi.slice(i);        // normalized
  }
}
//
/***********************************************************************************/
/* sampleBetas()                                                                   */
/***********************************************************************************/
void sampleBetas(double const& mode)
{ // Sample M = B'xi; B| Sigma, Theta, X
  // Mode: 1 NB; 2 NS; 3 SS; 
  mat S, V, Theta, Bn, Bsamp, Lmat;
  if(mode == 2.0){
    Lmat = cube2mat(NS.Ld);          //pq*m
  }
  Theta = cube2mat(th.mbi+th.bi);     //pq*n
  //
  S = g / (g+1.0) * th.XtXi;
  if(mode > 1.0){
    V = kron(pinv(th.tauq), pinv(th.taup));
    if(mode == 2.0){
      V = V + Lmat * trans(Lmat);
    } else {
      V = V + kron(SS.LR*SS.LR.t(), SS.LL * SS.LL.t());
    }
  } else {
    V = kron(pinv(NB.Df), pinv(NB.iS));
  }
  //
  Bn = S * trans(th.X) * trans(Theta); 
  //sample population means: beta;-------------------------------------------------/
  Bsamp = rMN(Bn, S, V);                                // pq*d
  th.beta = mat2cube(trans(Bsamp), dat.nr, th.nbase);   // p*q*d
  // Update mbi based on beta------------------------------------------------------/
  th.mbi = cubexmat(th.beta, th.X);                     // p*q*n
}
//
void sampleBeta0(){
  mat mbeta, betas, fit;
  for(r=0; r<dat.nr;r++){
    mbeta = th.XtXi * trans(th.X) * (trans(flatslice(th.bi,r))+trans(flatslice(th.mbi,r))); // d * q
    betas = rMN(mbeta, th.XtXi, pinv(NB.iS(r,r)*NB.Df));       // d * q
    for(j=0; j<th.nbase; j++){
      th.beta.tube(r,j) = conv_to<vec>::from(betas.col(j));
    }
  }
  th.mbi = cubexmat(th.beta, th.X);
}
/***********************************************************************************/
/* sampleLoadingsV()                                                                */
/***********************************************************************************/
void sampleLoadingsV()
{// Sample Lambda, Omega, delta | Theta, M, SigP, SigQ, Fac
  vec     samp, tau, b, copy;
  mat     Sp, sampM, W, T, sch;
  // cumsum for penalty update
  vec     bg = 1.0 * ones<vec>(NS.nfac);
  //
  double  nu = pNS.aOm + NS.nfac;
  // tau = prod pen's
  tau     = cumprod(NS.pen);
  T       = diagmat(tau);
  for(k=0; k<th.nbase; k++){
    vbi   = vectorise(longslice(th.bi, k));
    Sp    = th.taup * th.tauq(k,k);
    Sinv  = kron(T, NS.Omega.slice(k));
    Sinv  = Sinv + kron(NS.HtH, Sp);
    b     = kron(trans(NS.H), Sp) * vbi;
    //
    // Sample loadings, pm * 1
    //---------------------Safe Cholesky-------------------------//
    // sch   = safechol(Sinv);
    ch_ok   = chol(sch, Sinv, "lower"); 
    if(!ch_ok){Sinv.print("Sinv");stop("Loading Matrix sampler failed, chol() Error at iter " + std::to_string(iter));}
    samp  = BayesRegL(b, sch);
    // Matricize; p*m
    sampM = vec2mat(samp, dat.nr, NS.nfac);
    // Update Loading matrix ---------------------------------------------------- /
    vec2long(NS.Ld, samp, k);
    // Update Omega --------------------------------------------------------------/
    W     = sampM * T * trans(sampM);
    W     = pinv(W + pinv(pNS.Om0));
    NS.Omega.slice(k) = rWish(W, nu);
    // Update cumsum for penalties
    for(l=0; l<NS.nfac; l++){
      copy = NS.pen;     copy[l] = 1.0;
      tau = cumprod(copy);    T = diagmat(tau);
      bg[l] += 0.5 * trace(T.submat(l,l,NS.nfac-1,NS.nfac-1) * 
        trans(sampM.cols(l, NS.nfac-1)) * 
        NS.Omega.slice(k) * sampM.cols(l, NS.nfac-1));
    }
  }
  // shrinkage(penalty) update
  // delta0;----------------------------------------------------------------------/
  {
    double a1   = pNS.ap1 + dat.nr * th.nbase * NS.nfac * 0.5;
    NS.pen[0]  = Rf_rgamma(a1, 1.0/bg[0]);
  }
  // delta1+;---------------------------------------------------------------------/
  for(l=1; l<NS.nfac; l++){
    double a1   = pNS.ap2 + dat.nr * th.nbase * (NS.nfac-l) * 0.5;
    NS.pen[l]  = Rf_rgamma(a1, 1.0/bg[l]);
  }
}
//
/***********************************************************************************/
/* sampleLoadingsS()                                                                */
/***********************************************************************************/
void sampleLoadingsS()
{ // Sample LL LR, tauL tauR, penL penR | bi, taup tauq 
  // Variable declaration
  vec     pL(SS.nfacL), pR(SS.nfacR), sampL, sampR, ssL(SS.nfacL), ssR(SS.nfacR);
  mat     TL, TR, schL, schR;
  double  ga, gb;
  // Tau as prod of pen's
  pL    = cumprod(SS.penL); TL = diagmat(pL);
  pR    = cumprod(SS.penR); TR = diagmat(pR);
  // ------ Left Loading ------
  for(j=0; j<dat.nr; j++){
    vec temp = pL % conv_to<vec>::from(SS.tauL.row(j));
    vec b = zeros<vec>(SS.nfacL);
    mat S = zeros<mat>(SS.nfacL, SS.nfacL);
    for(i=0; i<dat.ns; i++){
      b += SS.H.slice(i) * SS.LR.t() * th.tauq * trans(th.bi.slice(i).row(j));
      S += SS.H.slice(i) * SS.LR.t() * th.tauq * SS.LR * SS.H.slice(i).t();
    }
    b = b * th.taup(j,j);
    S = S * th.taup(j,j) + diagmat(temp);
    //---------------------Safe Cholesky-------------------------//
    // schL   = safechol(S);
    ch_ok   = chol(schL, S, "lower"); 
    if(!ch_ok){S.print("S");stop("Loading Matrix sampler failed, chol() Error on LL at iter " + std::to_string(iter));}
    sampL = BayesRegL(b, schL);
    // Update loading matrix by rows;
    SS.LL.row(j) = conv_to<rowvec>::from(sampL);
    // Update loading precisions;
    for(r=0; r<SS.nfacL; r++){
      ga = 0.5 * (pSS.nuL + 1.0);
      gb = 0.5 * (pSS.nuL + pL[r] * std::pow(SS.LL(j,r), 2.0));
      SS.tauL(j,r) = Rf_rgamma(ga, 1.0/gb);
    }
  }
  // ------ Left Penalties ------
  ssL = conv_to<vec>::from(sum((SS.tauL % SS.tauL),0));
  for(r=0; r<SS.nfacL; r++){
    ga = (r==0) ? pSS.ap1 : pSS.ap2;
    ga += 0.5 * dat.nr * (SS.nfacL - r);
    vec tp = SS.penL;  tp[r] = 1.0;
    pL     = cumprod(tp);
    gb = 1.0;
    for(k=r; k<SS.nfacL; k++){
      gb += 0.5 * pL[k] * ssL[k];
    }
    SS.penL[r] = Rf_rgamma(ga, 1.0/gb);
  }
  //
  // ------ Right Loading ------
  for(j=0; j<th.nbase; j++){
    vec temp = pR % conv_to<vec>::from(SS.tauR.row(j));
    vec b = zeros<vec>(SS.nfacR);
    mat S = zeros<mat>(SS.nfacR, SS.nfacR);
    for(i=0; i<dat.ns; i++){
      b += SS.H.slice(i).t() * SS.LL.t() * th.taup * th.bi.slice(i).col(j);
      S += SS.H.slice(i).t() * SS.LL.t() * th.taup * SS.LL * SS.H.slice(i);
    }
    b = b * th.tauq(j,j);
    S = S * th.tauq(j,j) + diagmat(temp);
    //---------------------Safe Cholesky-------------------------//
    // schR   = safechol(S);
    ch_ok   = chol(schR, S, "lower"); 
    if(!ch_ok){S.print("S");stop("Loading Matrix sampler failed, chol() Error on RL at iter " + std::to_string(iter));}
    sampR  = BayesRegL(b, schR);
    // Update loading matrix by rows;
    SS.LR.row(j) = conv_to<rowvec>::from(sampR);
    // Update loading precisions;
    for(r=0; r<SS.nfacR; r++){
      ga = 0.5 * (pSS.nuR + 1.0);
      gb = 0.5 * (pSS.nuR + pR[r] * std::pow(SS.LR(j,r), 2.0));
      SS.tauR(j,r) = Rf_rgamma(ga, 1.0/gb);
    }
  }
  // ------ Reft Penalties ------
  ssR = conv_to<vec>::from(sum((SS.tauR % SS.tauR),0));
  for(r=0; r<SS.nfacR; r++){
    ga = (r==0) ? pSS.ap1 : pSS.ap2;
    ga += 0.5 * dat.nr * (SS.nfacR - r);
    vec tp = SS.penR;  tp[r] = 1.0;
    pR     = cumprod(tp);
    gb = 1.0;
    for(k=r; k<SS.nfacR; k++){
      gb += 0.5 * pR[k] * ssR[k];
    }
    SS.penR[r] = Rf_rgamma(ga, 1.0/gb);
  }
}
//
/***********************************************************************************/
/* sampleFacCovs()                                                                 */
/***********************************************************************************/
void sampleFacCovs(bool const& v)
{// Sample SigQ, Sigp | Theta, Lambda H
  double  ap, aq;
  cube    res;
  vec     bp(dat.nr), bq(th.nbase);
  //
  if(v){
    res = th.bi - cubexmat(NS.Ld, NS.H);               //Theta_i-M_i-mat(Lambda*H)
  } else {
    res = th.bi - matxcube(SS.LL, SS.H, SS.LR.t());
  }
  {
    ap = pNS.athe + 0.5 * dat.ns * th.nbase;
    aq = pNS.athe + 0.5 * dat.ns * dat.nr;
    bp.fill(pNS.bthe);
    bq.fill(pNS.bthe);
    mat sump, sumq;
    for(i=0; i<dat.ns; i++){
      sump = res.slice(i) * th.tauq * trans(res.slice(i));
      sumq = trans(res.slice(i)) * th.taup * res.slice(i);
      bp += 0.5 * diagvec(sump);
      bq += 0.5 * diagvec(sumq);
    }
  }
  // Sample col covs SigP;---------------------------------------------------------/
  for(j=0; j<dat.nr; j++){th.taup(j,j) = Rf_rgamma(ap, 1.0/bp[j]);}
  // Sample row covs SigQ;---------------------------------------------------------/
  for(k=0; k<th.nbase; k++){th.tauq(k,k) = Rf_rgamma(aq, 1.0/bq[k]);}
}
//
/***********************************************************************************/
/* sampleErrCovs()                                                                 */
/***********************************************************************************/
void sampleErrCovs(cube &y)
{// Sample sige | Yi, Theta
  double a, b;
  vec  h1;
  mat  resE, bi, mi, betar, S;
  
  // a    = pNS.ae + (dat.ns*dat.nt*dat.nr - accu(dat.nmiss)) * 0.5;
  a = (dat.ns*dat.nt - accu(dat.nmiss))*0.5 + pNS.ae;
  
  // Error sum of squares ---------------------
  for(r=0; r<dat.nr; r++){
    // Residuals
    resE = longslice(y, r) - longslice(th.fit, r);
    // b
    b    = pNS.be;
    for(i=0; i<dat.ns; i++){
      for(t=0; t<dat.nt; t++){
        b += 0.5 * ((dat.miss(i,r,t)==0)? std::pow(resE(i,t),2.0) : 0.0);
      }
    }
    //sample residual precision;  --------------------------------------------------/
    th.He(r,r) = Rf_rgamma(a, 1.0/b);
  }
  // th.taue = Rf_rgamma(a, 1.0/b);
}
// /***********************************************************************************/
// /* ErgodicAverages()                                                               */
// /***********************************************************************************/
// void ErgodicAverages(mat const& Bs)
// {
//   int p;
//   mat beta, betat;                                                    // nt*np
//   // output as nt * nr * np
//   w1 = 1.0/(1.0 + hat.n);
//   w2 = 1.0 - w1;
//   hat.n++;
//   // Average fit ------------------------------------------
//   hat.fit  = w1*th.fit + w2*hat.fit;
//   // Average regression coefficients ----------------------
//   for(p=0; p<dat.np; p++){
//     beta  = flatslice(th.beta, p);                                     // nbase * nr
//     betat = Bs * beta;                                                 // nt x nr
//     hat.beta.slice(p)  = w1*betat + w2*hat.beta.slice(p);              // nt x nr
//     hat.beta2.slice(p) = w1*betat%betat + w2*hat.beta2.slice(p);       // nt x nr
//   }
// }
// //
/***********************************************************************************/
/* sampleFactors()                                                                 */
/***********************************************************************************/
void sampleFactors(cube &y, mat const& Bs, bool const& v)
{// Sample factors | everything else;
  mat     lmat, res, Q, Stinv, Sipq, fit, qch, sch;
  cube    LFcube;
  vec     b, bi, thi;
  int     ind;
  // -- common variables:
  Sinv =  kron(pinv(th.tauq), pinv(th.taup));
  // Sinv += kron((th.iBtB), 1.0/th.taue*eye<mat>(dat.nr, dat.nr));
  Sinv += kron(th.iBtB, pinv(th.He));
  mat pSinv = pinv(Sinv);
  if(v){
    lmat    = cube2mat(NS.Ld);               //loading matrix
    LFcube  = cubexmat(NS.Ld, NS.H);        //latent factor
    Q = trans(lmat) * pSinv * lmat + 1.0 * eye<mat>(NS.nfac, NS.nfac);
  } else{
    LFcube  = matxcube(SS.LL, SS.H, SS.LR.t());//latent factor
    lmat    = zeros<mat>(dat.nr*th.nbase, SS.nfacL*SS.nfacR);
    for(r=0; r<SS.nfacL;r++){
      for(s=0; s<SS.nfacR; s++){
        ind = r*SS.nfacR + s;
        lmat.col(ind) = vectorise(SS.LL.col(r)*SS.LR.col(s).t());
      }
    }
    Q = trans(lmat) * pSinv * lmat + 1.0 * eye<mat>(SS.nfacL*SS.nfacR, 
              SS.nfacL*SS.nfacR);
  }
  // qch     = safechol(Q);
  ch_ok   = chol(qch, Q, "lower"); 
  if(!ch_ok){Q.print("Q");Rcpp::stop("Factors sampler failed, chol() Error at iter " + std::to_string(iter));}
  Sipq    = kron(th.taup, th.tauq);
  // Stinv   = Sipq + kron(th.taue * eye<mat>(dat.nr, dat.nr), th.BtB);
  Stinv   = Sipq + kron(th.He, th.BtB);
  // sch     = safechol(Stinv);
  ch_ok   = chol(sch, Stinv, "lower"); 
  if(!ch_ok){Stinv.print("Stinv");stop("Factors sampler failed, chol() Error on theta_i at iter " + std::to_string(iter));}
  //
  for(i=0; i<dat.ns; i++){
    //residual of pq
    res = (flatslice(y,i) * Bs * th.iBtB - th.mbi.slice(i));
    b   = (v) ? pNS.eta0 : pSS.eta0;
    b   = trans(lmat) * pSinv * vectorise(res);
    // Sample latent Factors;--------------------------------------------------------/
    vec bsamp   = BayesRegL(b,qch);
    if(v){
      NS.H.row(i) = conv_to<rowvec>::from(bsamp);
    } else {
      SS.H.slice(i) = vec2mat(bsamp, SS.nfacR, SS.nfacL).t();
    }
    // residual of Theta_i
    // bi  = vectorise(th.taue * trans(Bs) * trans(flatslice(y,i)));
    bi  = vectorise(trans(Bs) * trans(flatslice(y,i)) * th.He);
    bi  = bi + Sipq * vectorise(trans(th.mbi.slice(i) + LFcube.slice(i)));
    // Sample coefficient matrices --------------------------------------------------/
    thi = BayesRegL(bi, sch);
    // Update bi = Theta_i - mbi
    th.bi.slice(i) = trans(vec2mat(thi, th.nbase, dat.nr)) - th.mbi.slice(i);
    fit = (th.bi.slice(i) + th.mbi.slice(i)) * trans(Bs);       //p*nt      
    // Update fit_i = Theta_i * Bs'
    for(r=0; r<dat.nr; r++){
      for(t=0; t<dat.nt; t++){
        th.fit(i,r,t) = fit(r,t);
      }
    }
  }
  // Update H'H;
  NS.HtH = trans(NS.H) * NS.H;
}
//
/***********************************************************************************/
/* OutputSample()                                                                  */
/***********************************************************************************/
void OutputSample(double const& mode)
{
  { // Population Splines --------------------------------
    std::stringstream tempstr;
    for(p=0; p<dat.np; p++){
      for(r=0; r<dat.nr; r++){
        for(j=0; j<th.nbase; j++){
          tempstr << th.beta(r,j,p) << " ";
        }
      }
    }
    tempstr << "\n"; out1 << tempstr.str();
  }
  { // Precisions ----------------------------------------
    std::stringstream tempstr;
    for(r=0; r<dat.nr; r++){
      tempstr << th.He(r,r) << " "; // Error precision
    }
    if(mode > 1.0){
      for(r=0; r<dat.nr; r++){
        tempstr << th.taup(r,r) << " ";   // Spatial res covs;
      }
      for(j=0; j<th.nbase; j++){
        tempstr << th.tauq(j,j) << " ";   // Temp res covs;
      }
    }
    tempstr << "\n"; out2 << tempstr.str();
  }
  { // Spatial or Latent Factors -------------------------------------------
    std::stringstream tempstr;
    if(mode == 1.0){
      mat S  = (pinv(NB.iS+0.0));
      mat SS = cov2cor(S);
      for(r=0; r<dat.nr; r++){
        for(j=r; j<dat.nr; j++){
          tempstr << SS(r,j) << " "; // Spatial matrix
        }
      }
    } else if (mode == 2.0){
      for(j=0; j<NS.nfac; j++){
        tempstr << NS.H(r,j) << " "; // LF scores
      }
    } else if (mode == 3.0){
      for(j=0; j<SS.nfacR; j++){
        for(k=0; k<SS.nfacL; k++){
          tempstr << SS.H(k,j,r) << " "; // LF-s scores
        }
      }
    }
    tempstr << "\n"; out3 << tempstr.str();
  }
}
//
/***********************************************************************************/
/* ErgodicAverages()                                                               */
/***********************************************************************************/
void ErgodicAverages(mat const& Bs, double const& mode)
{ //
  mat betat; //, beta2;
  mat V, Lmat;
  //
  w1 = 1.0/(1.0 + hat.n);
  w2 = 1.0 - w1;
  hat.n++;
  //
  // Average fit ------------------------------------------
  hat.fit  = w1 * th.fit + w2 * hat.fit;                   //p*nt individual fit
  // Average Coeffients covariances;
  if(mode == 1.0) {
    hat.covcoef = w2 * hat.covcoef + w1 * kron(pinv(NB.Df),pinv(NB.iS));
  } else {
    V = kron(pinv(th.tauq), pinv(th.taup));
    if(mode == 2.0){
      Lmat = cube2mat(NS.Ld);
      V += Lmat * Lmat.t();
      hat.LF = w2 * hat.LF + w1 * NS.H;
      hat.LF2 = w2 * hat.LF2 + w1 * NS.H%NS.H;
    } else if (mode == 3.0){
      V += kron(SS.LR * SS.LR.t(), SS.LL * SS.LL.t());
      mat LF = cube2mat(SS.H);
      hat.LF = w2 * hat.LF + w1 * LF;
      hat.LF2 = w2 * hat.LF2 + w1 * LF % LF;
    }
    hat.covcoef = w2 * hat.covcoef + w1 * V;
  }
  //
  for(i=0; i<dat.np; i++){
    betat     = Bs * th.beta.slice(i).t();
    hat.beta.slice(i) = w1 * betat + w2 * hat.beta.slice(i);        //nt * nr * nvars
    hat.beta2.slice(i)= w1 * betat%betat + w2 * hat.beta2.slice(i); //nt * nr * nvars
  }
}
//
/***********************************************************************************/
/* init_HFM()                                                                      */
/***********************************************************************************/
//
void init_HFM(arma::cube y,
               arma::mat const& X,
               arma::mat const& Bs)
{
  // Get initialization summaries ---------------------------------------------------
  dat.ns = y.n_rows;           // number of subjects
  dat.nr = y.n_cols;           // number of regions
  dat.nt = y.n_slices;         // number of segments
  dat.np = X.n_cols;           // number of covariates (groups)
  //
  th.nbase = Bs.n_cols;        // number of basis functions
  //
  // Initialize missing data ----------------------------------------
  {
    vec rn1;
    dat.miss = icube(dat.ns, dat.nr, dat.nt);
    dat.nmiss= vec(dat.ns);
    yit = randu<vec>(dat.ns);
    for(i=0; i<dat.ns; i++){
      for(r=0; r<dat.nr; r++){
        yit = y.tube(i,r);
        for(j=0; j<dat.nt; j++){
          dat.miss(i,r,j) = 0;
          if( yit(j)==12345.12345 ){
            dat.miss(i,r,j) = 1;
            rn1 = randn(1)*0.5;
            y(i,r,j) = meanNA(yit) + rn1(0);
          }
        }
      }
      dat.nmiss[i] = accu(dat.miss.tube(0,0));
    }
  }
  // Initialize and decompose regression matrices -------------------
  th.X     = X;                         // X
  th.XtXi  = pinv(trans(X)*X);          // pinv(X'X)
  // Random initialization;
  th.bi   = randu<cube>(dat.nr, th.nbase, dat.ns);    // basis functions
  th.mbi  = randu<cube>(dat.nr, th.nbase, dat.ns);    // mean bases funcs
  th.fit  = randu<cube>(dat.ns, dat.nr, dat.nt);       // current fit
  th.beta = randu<cube>(dat.nr, th.nbase, dat.np);    // Group means
  NS.Omega = randu<cube>(dat.nr, dat.nr, th.nbase);    // Spatial correlations
  th.BtB = trans(Bs) * Bs;
  th.iBtB= pinv(th.BtB);
  th.taup= 1.0 * eye(dat.nr, dat.nr);
  th.tauq= 1.0 * eye(th.nbase, th.nbase);
  // th.taue= Rf_rgamma(ae, 1.0/be);                         // Residual vars
  th.He = 0.1 * eye(dat.nr, dat.nr);
  // Meaningful initialization:                       // ----------------------
  bir    = randu<vec>(th.nbase);
  for(i=0; i<dat.ns; i++){
    for(r=0; r<dat.nr; r++){
      yit = y.tube(i,r);
      bir = solve((th.BtB + th.He(r,r) * eye<mat>(th.nbase, th.nbase)), trans(Bs) * yit);
      for(j=0; j<th.nbase; j++){
        th.bi(r,j,i) = bir[j];
      }
      th.fit.tube(i,r) = vectorise(Bs*bir);
    }
  }
  for(i=0; i<dat.ns; i++){th.mbi.slice(i) = mean(th.bi,2);}
  int nn1 = 1.0 * dat.ns* dat.nt * dat.nr;
  for(r=0; r<dat.nr; r++){
    double err = 0.0;
    for(i=0; i<dat.ns; i++){
      for(j=0; j<dat.nt; j++){
        err += pow(y(i,r,j) - th.fit(i,r,j), 2.0)/nn1;
      } // j
    } //i
    th.He(r,r) = 3.0/err;
  }
  // Initialize ergodics --------------------------------------------
  hat.n     = 0.0;
  hat.beta  = randu<cube>(dat.nt, dat.nr, dat.np);
  hat.beta2 = randu<cube>(dat.nt, dat.nr, dat.np);
  hat.fit   = randu<cube>(dat.ns, dat.nr, dat.nt);
  hat.covcoef=randu<mat>(dat.nr*th.nbase,dat.nr*th.nbase);
}
/***********************************************************************************/
/* NBmix_mcmc()                                                                     */
/***********************************************************************************/

// [[Rcpp::export]]
List NBmix_mcmc(arma::cube y,
               arma::mat const& X,
               arma::mat const& Bs,
               int  const& burnin, int const& nsim, int const& thin,
               bool const& spatial = true) {
  //
  mat W, D, D2;               // AR matrices
  vec rho;                    // regularization eigenvalue
  /* ----------------------------------------------------------------------------- */
  /* Initialize Chain and static summaries                                         */
  /* ----------------------------------------------------------------------------- */
  init_HFM(y, X, Bs);
  // Prior parameters -----------------------------------------------
  pNB.ae = 0.0001;         pNB.be = 0.0001;               // Error variance
  pNB.ab = 0.1;  pNB.bb = 0.1;                            // Spline-coefs precision
  pNB.aS = 0.0 + dat.nr;
  pNB.W = eye<mat>(dat.nr, dat.nr) * 1.0/pNB.aS;          // spatial smoothing
  //
  // Smoothing Matrix (AR-1) ----------------------------------------
  W  = eye<mat>(th.nbase, th.nbase)*1.0;
  D  = eye<mat>(th.nbase, th.nbase)*3.0;
  D2 = eye<mat>(th.nbase, th.nbase)*1.0;
  // AR1 adjecency -----------------------------------
  for(i=1; i<th.nbase; i++){ W(i,(i-1)) = 1.0; }
  for(i=0; i<(th.nbase-1); i++){ W(i,(i+1)) = 1.0; }
  // AR1 no. neighbors -------------------------------
  D(0,0) = D((th.nbase-1),(th.nbase-1)) = 2.0;
  // Regularizing eigenvalue -------------------------
  for(i=0; i<th.nbase; i++) D2(i,i) = std::pow(D(i,i), -0.5);
  rho = randu<vec>(th.nbase);
  eig_sym(rho, D2*W*D2);
  // Final penalization matrix -----------------------
  NB.Df = (D - rho[(th.nbase-2)]*W)*10.0;
  // th.Df = eye<mat>(th.nbase, th.nbase)*1.0;
  //
  // Individual functional means and smoothing ----------------------
  NB.lbi = 1.0;                                   // precision
  NB.iS = eye<mat>(dat.nr, dat.nr)*1.0;
  // Population means means and smoothing ---------------------------
  NB.lbeta = 0.01/(dat.ns*dat.nr);
  // Initialize and decompose regression matrices -------------------
  NB.cholXX = chol(th.XtXi);          // Cholesky of (X'X)^-1 (upper)
  NB.cholDf = chol(pinv(NB.Df));      // Cholesky of Df^-1 (upper)
  //
  // Open output files ----------------------------------------------
  out1.open(fileOutput1);
  out2.open(fileOutput2);
  out3.open(fileOutput3);
  //
  /* ----------------------------------------------------------------------------- */
  /* Gibbs Sampling                                                                */
  /* ----------------------------------------------------------------------------- */
  //
  for(iter=0; iter<(nsim+burnin+1); iter++){
    completeY(y);
    samplePrecisions(y, spatial);
    samplebi(Bs, y);
    sampleBetas(1.0);  //Unified sampler
    // sampleBeta0();   //old sampler
    // Store mcmc samples --------------
    if( (iter > burnin) && (iter%thin == 0) ){
      OutputSample(1.0);           // Write mcmc samples to file
      ErgodicAverages(Bs,1.0);      // Update ergodic averages
    }
  }
  //
  // Close output files ---------------------------------------------
  out1.close(); out2.close(); out3.close();
  //
  // Output to R --------------------------------------------------------------------
  return List::create(
    Named("fit")      = hat.fit,
    Named("coef")     = hat.beta,
    Named("coef2")    = hat.beta2,
    Named("covcoef")  = hat.covcoef);
}
// Vectorized HFM
/***********************************************************************************/
/* NSmix_mcmc()                                                                     */
/***********************************************************************************/
// [[Rcpp::export]]
List NSmix_mcmc(arma::cube y, 
               arma::mat const& X, 
               arma::mat const& Bs,
               int const& nfac,
               int const& burnin, int const& nsim, int const& thin
)
{ 
  // number of factors in v;
  NS.nfac  = nfac;             // number of latent factors
  /* ----------------------------------------------------------------------------- */
  /* Initialize Chain and static summaries                                         */
  /* ----------------------------------------------------------------------------- */
  init_HFM(y, X, Bs);
  // Prior parameters -----------------------------------------------
  pNS.ae = 0.0001;          pNS.be = 0.0001;              // Error variance
  pNS.athe = 0.001;         pNS.bthe = 0.001;             // coeff res row/col var
  pNS.ap1 = 1.5;            pNS.ap2 = 1.5;                // penalties
  pNS.aOm = 0.0 + dat.nr;                                // --------------could be less informative
  pNS.Om0 = eye<mat>(dat.nr, dat.nr) * 1.0/pNS.aOm;     // spatial smoothing
  pNS.eta0 = zeros<vec>(NS.nfac);                        // Latent factor priors
  //
  // RANDOM initialization                              // ----------------------
  NS.H   = randu<mat>(dat.ns, nfac);                    // Latent factors
  NS.Ld  = randu<cube>(dat.nr, th.nbase, nfac);         // Loading matrix
  NS.HtH = trans(NS.H) * NS.H;
  NS.pen = randu<vec>(nfac);
  for(i=0; i<nfac; i++){NS.pen[i] = Rf_rgamma(1*pNS.ap1, 1.0/pNS.ap2);}
  for(j=0; j<th.nbase; j++){NS.Omega.slice(j) = 1.0*eye<mat>(dat.nr, dat.nr);}
  hat.LF = NS.H;
  hat.LF2= NS.H;
  //
  // Open output files ----------------------------------------------
  out1.open(fileOutput4);
  out2.open(fileOutput5);
  out3.open(fileOutput6);
  /* ----------------------------------------------------------------------------- */
  /* Gibbs Sampling                                                                */
  /* ----------------------------------------------------------------------------- */
  for(iter=0; iter<(nsim+burnin+1); iter++){
    // mcmc samples --------------------
    completeY(y);
    sampleLoadingsV();
    sampleFacCovs(true);
    sampleErrCovs(y);
    sampleBetas(2.0);
    sampleFactors(y, Bs, true);
    // Store mcmc samples --------------
    if( (iter > burnin) && (iter % thin == 0) ){
      OutputSample(2.0);             // Write mcmc samples to file
      ErgodicAverages(Bs,2.0);      // Update ergodic averages
    }
  }
  // Close output files ---------------------------------------------
  out1.close(); out2.close(); out3.close(); 
  // Output to R --------------------------------------------------------------------
  return List::create(
    Named("fit")      = hat.fit,
    Named("coef")     = hat.beta,
    Named("coef2")    = hat.beta2,
    Named("covcoef")  = hat.covcoef,
    Named("LF")       = hat.LF,
    Named("LF2")      = hat.LF2
  );
}
// Sandwich HFMM
/***********************************************************************************/
/* SSmix_mcmc()                                                                     */
/***********************************************************************************/
// [[Rcpp::export]]
List SSmix_mcmc(arma::cube y, 
               arma::mat const& X, 
               arma::mat const& Bs,
               int const& nfacL, int const& nfacR,
               int const& burnin, int const& nsim, int const& thin
)
{
  // number of factors in v;
  SS.nfacL  = nfacL; SS.nfacR  = nfacR;                // number of latent factors
  /* ----------------------------------------------------------------------------- */
  /* Initialize Chain and static summaries                                         */
  /* ----------------------------------------------------------------------------- */
  init_HFM(y, X, Bs);
  // Prior parameters -----------------------------------------------
  pSS.ae = 0.0001;          pSS.be = 0.0001;              // Error variance
  pSS.athe = 0.001;         pSS.bthe = 0.001;             // coeff res row/col var
  pSS.ap1 = 1.5;            pSS.ap2 = 1.5;                // penalties
  pSS.nuL = 5.0;            pSS.nuR = 5.0;                // LF precision
  pSS.eta0 = zeros<vec>(nfacL*nfacR);                     // Latent factor priors
  //
  // RANDOM initialization                                // ----------------------
  SS.H   = randu<cube>(nfacL, nfacR, dat.ns);           // Latent factors
  SS.LL  = randu<mat>(dat.nr, nfacL);                   // Loading matrix (left)
  SS.LR  = randu<mat>(th.nbase, nfacR);                // Loading matrix (right)
  SS.penL = randu<vec>(nfacL);
  for(i=0; i<nfacL; i++){SS.penL[i] = Rf_rgamma(pSS.ap1, 1.0/pSS.ap2);}
  SS.penR = randu<vec>(nfacR);
  for(r=0; r<nfacR; r++){SS.penR[r] = Rf_rgamma(pSS.ap1, 1.0/pSS.ap2);}
  SS.tauL= randu<mat>(dat.nr, nfacL);
  SS.tauR= randu<mat>(th.nbase, nfacR);
  for(i=0; i<SS.tauL.n_elem; i++){SS.tauL[i] = Rf_rgamma(pSS.nuL, 1.0/pSS.nuL);}
  for(j=0; j<SS.tauR.n_elem; j++){SS.tauR[j] = Rf_rgamma(pSS.nuR, 1.0/pSS.nuR);}
  hat.LF  = randu<mat>(nfacL*nfacR, dat.ns);
  hat.LF2 = hat.LF;
  //
  // Open output files ----------------------------------------------
  out1.open(fileOutput7);
  out2.open(fileOutput8);
  out3.open(fileOutput9);
  /* ----------------------------------------------------------------------------- */
  /* Gibbs Sampling                                                                */
  /* ----------------------------------------------------------------------------- */
  for(int rep=0; rep<(nsim+burnin+1); rep++){
    // mcmc samples --------------------
    completeY(y);
    sampleLoadingsS();
    sampleFacCovs(false);
    sampleErrCovs(y);
    sampleBetas(3.0);
    sampleFactors(y, Bs, false);
    // Store mcmc samples --------------
    if( (rep > burnin) && (rep%thin == 0) ){
      OutputSample(3.0);             // Write mcmc samples to file
      ErgodicAverages(Bs, 3.0);      // Update ergodic averages
    }
  }
  // Close output files ---------------------------------------------
  out1.close(); out2.close(); out3.close();
  
  // Output to R --------------------------------------------------------------------
  return List::create(
    Named("fit")      = hat.fit,
    Named("coef")     = hat.beta,
    Named("coef2")    = hat.beta2,
    Named("covcoef")  = hat.covcoef,
    Named("LF")       = hat.LF,
    Named("LF2")      = hat.LF2
  );
}