#include "HFM.h"
// disable printing of error messages
std::ostream nullstream(0);
//
/***********************************************************************************/
/* meanNA()                                                                        */
/***********************************************************************************/
double meanNA(vec &y)
{
  int i, n, ny;
  double m;
  ny = y.n_elem;
  m  = 0.0;
  n  = 0;
  for(i=0; i<ny; i++){
    if(y(i) != 12345.12345){
      m += y(i);
      n ++;
    }
  }
  return m/n;
}


/***********************************************************************************/
/* longslice()                                                                     */
/***********************************************************************************/
mat longslice(cube const &A, int r)
{
  int i, j;
  mat b(A.n_rows, A.n_slices); b.fill(0.0);
  for(i=0; i<A.n_rows; i++){
    for(j=0; j<A.n_slices; j++){
      b(i,j) = A(i,r,j);
    }
  }
  return b;
}
/***********************************************************************************/
/* flatslice()                                                                     */
/***********************************************************************************/
mat flatslice(cube const &A, int i)
{
  int r, j;
  mat b(A.n_cols, A.n_slices); b.fill(0.0);
  for(j=0; j<A.n_cols; j++){
    for(r=0; r<A.n_slices; r++){
      b(j,r) = A(i,j,r);
    }
  }
  return b;
}


/***********************************************************************************/
/* trisolve()                                                                      */
/* Sole triangular system Lx = b                                                   */
/***********************************************************************************/
void trisolve(vec &solution, mat const &L, vec const &x, bool const &lower)
{
  int i,j;
  double sum;
  int n = x.n_elem;
  //
  if(lower){ // Lower (forward) ---------------------
    solution(0) = x(0)/L(0,0);
    for(i=1; i<n; i++){
      sum = 0.0;
      for(j=0; j<i; j++)  sum +=  L(i,j)*solution(j);
      solution(i) = (x(i) - sum)/L(i,i);
    }
  }
  if(!lower){ // Upper (backward) --------------------
    solution((n-1)) = x((n-1))/L((n-1),(n-1));
    for(i=(n-2); i>=0; i--) {
      sum = 0.0;
      for(j=i+1; j<n; j++) sum +=  L(i,j)*solution(j);
      solution(i) = (x(i) - sum)/L(i,i);
    }
  }
}
// /***********************************************************************************/
// /* BayesReg()                                                                      */
// /* Sample x ~ N_p(Q^-1b, Q^-1)                                                     */
// /***********************************************************************************/
// vec BayesReg(vec const &b, mat const &Q)
// {
//   vec z, x, m;
//   int p = b.n_elem;
//   x = zeros<vec>(p); m = zeros<vec>(p);
//   //
//   mat U = chol(Q, "upper");       // Q = U'U : U is upper triang
//   trisolve(x, trans(U), b, true); // U'(Um) = b; Um = x
//   trisolve(m, U, x, false);       // m = Q^-1 b
//   z = randn<vec>(p);              // z ~ N_p(0, I)
//   trisolve(x, U, z, false);       // x ~ N_p(0, Q^-1)
//   z = x + m;                      // z ~ N_p(b, Q^-1)
//   return z;
// }
/***********************************************************************************/
/* BayesRegL()                                                                     */
/* Sample x ~ N_p((LL')^-1b, (LL')^-1)                                             */
/***********************************************************************************/
vec BayesRegL(vec const &b, mat const &L)
{ // L is lower triang
  vec z, x, m;
  int p = b.n_elem;
  //
  x = zeros<vec>(p); m = zeros<vec>(p);
  //
  trisolve(x, L, b, true);            //
  trisolve(m, trans(L), x, false);    // m = Q^-1 b
  z = randn<vec>(p);                  // z ~ N_p(0, I)
  trisolve(x, trans(L), z, false);    // x ~ N_p(0, Q^-1)
  z = x + m;                          // z ~ N_p(b, Q^-1)
  return z;
}
/***********************************************************************************/
/* rWish()                                                                         */
/* Sample W ~ Wishart(S, v)                                                        */
/* Mean: vS; Mode (v-p-1)S                                                         */
/***********************************************************************************/
mat rWish(mat const &V, double v)
{
  int m = V.n_rows;
  mat T = zeros<mat>(m,m);
  // Barlett decomposition;
  for(int i=0; i<m; i++) {
    T(i,i) = std::pow(rchisq(1,v-i)[0], 0.5);
    for(int j=i+1; j<m; j++){
      T(i,j) = rnorm(1)[0];
    }
  }
  mat C = T * chol(V, "upper");
  return trans(C) * C ;
}
/***********************************************************************************/
/* rMN()                                                                           */
/* Sample X ~ MN(m, S, V)                                                          */
/***********************************************************************************/
mat rMN(mat const &m, mat const &S, mat const& V)
{
  int p, k;
  mat z, Sc, Vc;
  //
  p  = S.n_cols;
  k  = V.n_cols;
  //
  Sc = chol(S, "lower");                    // Cholesky S = ScSc'
  Vc = chol(V, "upper");                    // Cholesky V = Vc'Vc
  //
  z  = randn<mat>(p, k);                    // MN(0, I_p, I_k)
  return m + Sc*z*Vc;                       // MN(m, Sc'Sc, Vc'Vc)
}
/***********************************************************************************/
/* cov2cor()                                                                       */
/* Transform a covariance matrix W into a correlation matrix R                     */
/***********************************************************************************/
mat cov2cor(mat S)
{
  int m = S.n_cols;
  mat D = eye<mat>(m, m);
  for(int i=0; i<m; i++) D(i,i) = pow(S(i,i), -0.5);
  return D * S * D;
}
/***********************************************************************************/
/* vec2mat()                                                                       */
/* column-wise matricization                                                       */
/***********************************************************************************/
mat vec2mat(vec &v, int const& nr, int const& nc)
{
  int i, j, index;
  mat out(nr,nc);
  for(j=0; j<nc; j++){
    for(i=0; i<nr; i++){
      index = nr*j+i;
      out(i,j) = v[index];
    }
  }
  return out;
}
/***********************************************************************************/
/* vec2long()                                                                      */
/* feeding vector to longslice of cube a(i,j,k)                                    */
/***********************************************************************************/
void vec2long(cube &a, vec const &v, int const& j)
{
  int i, k, index;
  int nr = a.n_rows;
  int ns = a.n_slices;
  for(k=0; k<ns; k++){
    for(i=0; i<nr; i++){
      index = nr*k+i;
      a(i,j,k) = v[index];
    }
  }
}
/***********************************************************************************/
/* cubexmat()                                                                      */
/* cube(a,b,c) multiplies mat(d,c)                                                 */
/***********************************************************************************/
cube cubexmat(cube const&c, mat const&m)
{
  int nr,nc,ns,i,j,k;
  mat tmp, bmat;
  nr = c.n_rows;
  nc = c.n_cols;
  ns = m.n_rows;
  cube out(nr,nc,ns);
  for(j=0;j<nc;j++){
    bmat = longslice(c, j); // a*c
    tmp  = bmat * trans(m); // a*d
    for(i=0;i<nr;i++){
      for(k=0;k<ns;k++){
        out(i,j,k) = tmp(i,k);
      }
    }
  }
  return out;
}
/***********************************************************************************/
/* matxcube()                                                                      */
/* cube(a,b,c) multiplies mat(d,c)                                                 */
/***********************************************************************************/
cube matxcube(mat const& L, cube const&C, mat const&R)
{
  int nr = L.n_rows;
  int nc = R.n_cols;
  int ns = C.n_slices;
  cube out(nr,nc,ns);
  for(int i=0; i<ns; i++){
    out.slice(i) = L * C.slice(i) * R;
  }
  return out;
}
/***********************************************************************************/
/* cube2mat()                                                                      */
/* cube(a,b,c) into ab * c                                                         */
/***********************************************************************************/
mat cube2mat(cube const &c)
{
  int nr = c.n_rows;
  int nc = c.n_cols;
  int ns = c.n_slices;
  mat out(nr*nc, ns);
  vec temp(nr*nc);
  for(int i=0; i<ns; i++){
    temp = vectorise(c.slice(i));
    out.col(i) = conv_to< colvec >::from(temp);
  }
  return out;
}
/***********************************************************************************/
/* mat2cube()                                                                      */
/* mat(a*b,c) into cube(a,b,c)                                                     */
/***********************************************************************************/
cube mat2cube(mat const& m, int const& nr, int const& nc)
{ 
  int ns = m.n_cols;
  cube out(nr,nc,ns);
  vec temp;
  for(int s=0; s<ns; s++){
    temp = vectorise(m.col(s));
    out.slice(s) = vec2mat(temp, nr, nc);
  }
  return out;
}
/***********************************************************************************/
/* debugpinv()                                                                     */
/* pseudo inverse with exceptions safely                                           */
/***********************************************************************************/
mat debugpinv(mat m)
{ // Try inv_sympd first
  bool success = false;
  mat out;
  while(!success){
    set_cerr_stream(nullstream);
    success = inv_sympd(out, m);
    if(!success){
      Rcpp::warning("inv_sympd failed!!!!");
      if(m.has_inf()){
        // Inf-containing, replace
        m.replace(datum::inf, 1e8);
      } else if (m.has_nan()){
        // NaN-containing, relace
        m.replace(datum::nan, 1e-8);
      } else {
        // Singularity, perturb the diagonal
        m += (eye(size(m)) * 1e-4);
      }
    }
  }
  return out;
}
/***********************************************************************************/
/* diaginv()                                                                      */
/* Inverting a diagonal matrix                                                    */
/***********************************************************************************/
mat diaginv(mat const&m)
{ // Diagonal processing (Faster than inv_sympd)
  return(diagmat(1.0/diagvec(m)));
}
//