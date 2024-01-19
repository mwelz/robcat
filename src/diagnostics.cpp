#include <Rcpp.h>
#include <math.h> // for std::sqrt() and pow()
#include "matrixalgebra.h" // relevant helpers
#include "polycor_variance.h"
using namespace Rcpp;

// right quantile of standard normal distribution
//[[Rcpp::export]]
double pnorm_right(double x)
{
  Rcpp::Environment stats = Rcpp::Environment::namespace_env("stats");
  Rcpp::Function f = stats["pnorm"];
  NumericVector p = f(x, 0., 1., false, false);
  return p[0];
}


// p-value of two-sided test
//[[Rcpp::export]]
double pval_twosided(double z)
{
  return 2.0 * pnorm_right(std::fabs(z));
}


// p-value of one-sided test with right-hand alternative H_0: Z > z
//[[Rcpp::export]]
double pval_right(double z)
{
  return pnorm_right(z);
}


// test if a cell us outlying
//[[Rcpp::export]]
List celltest_cpp(
  double rho, 
  NumericVector thresX,
  NumericVector thresY,
  int Kx, int Ky,
  NumericMatrix probs,
  NumericMatrix f,
  NumericMatrix sigma,
  int N,
  bool twosided)
{
  int x, y;
  NumericMatrix stderr(Kx,Ky), teststat(Kx,Ky), pval(Kx,Ky);
  double sigma2, stderr_xy, z, p;
  double sqrtN = std::sqrt((double)N);
  
  // declare function pointer
  double (*p_fun)(double);
  
  // depending if test is 1- or 2-sided, call correct p-value function
  if(twosided)
  {
    p_fun = &pval_twosided;
  } else{
    p_fun = &pval_right;
  } // IF

  
  for(int i = 0; i < Kx; ++i)
  {
    x = i + 1;
    for(int j = 0; j < Ky; ++j)
    {
      y = j + 1;
      
      // calculate gradient at cell
      NumericVector grad = pk_prime_theta(x, y, thresX, thresY, rho, Kx, Ky);
      
      // cast to matrix
      NumericMatrix gradm = mat2vec(grad);
      
      /// delta rule variance (scalar)
      NumericMatrix tmp = mmult(transpose(gradm), sigma);
      sigma2 = mmult(tmp, gradm)(0,0);
      stderr_xy = std::sqrt(sigma2);
      
      // test statistic at cell
      z = sqrtN * (probs(i,j) - f(i,j)) / stderr_xy;
      
      // p-value of test statistic (unadjusted for multiple comparisons)
      p = p_fun(z);
      
      // collect results
      stderr(i,j) = stderr_xy;
      teststat(i,j) = z;
      pval(i,j) = p;
    } // i
  } // j
  
  return List::create(
    Rcpp::Named("teststat") = teststat,
    Rcpp::Named("stderr") = stderr,
    Rcpp::Named("pval") = pval
  );
}
