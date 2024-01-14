#include <Rcpp.h>
#include "matrixalgebra.h"
using namespace Rcpp;


// outer product of a vector
// [[Rcpp::export]]
NumericMatrix outer_vec(NumericVector v)
{
  // initialize
  int n = v.size();
  NumericMatrix m(n);
  double z;
  
  for(int i = 0; i < n; i++)
  {
    for(int j = 0; j < n; j++)
    {
      z = v[i] * v[j];
      m(i,j) = z;
      m(j,i) = z;
    }
  }
  return m;
}

// cbind() in C++
// we don't ever use this function; consider dropping in!
// taken from https://stackoverflow.com/a/31921185, thanks to Tony Tonov + Romain Francois!
// NumericMatrix cbind_cpp(NumericMatrix a, NumericMatrix b)
// {
//   int acoln = a.ncol();
//   int bcoln = b.ncol();
//   NumericMatrix out = no_init_matrix(a.nrow(), acoln + bcoln);
//   for (int j = 0; j < acoln + bcoln; j++) {
//     if (j < acoln) {
//       out(_, j) = a(_, j);
//     } else {
//       out(_, j) = b(_, j - acoln);
//     }
//   }
//   return out;
// }


// matrix addition: A + B
// A and B are of same dimension and square
// [[Rcpp::export]]
NumericMatrix matplus(NumericMatrix A, NumericMatrix B)
{
  int n = A.ncol();
  NumericMatrix m(n);
  
  for(int i = 0; i < n; i++)
  {
    for(int j = 0; j < n; j++)
    {
      m(i,j) = A(i,j) + B(i,j);
    }
  }
  return m;
}


// matrix multiplication (faster than with Armadillo)
//[[Rcpp::export]]
NumericMatrix mmult(NumericMatrix A, NumericMatrix B)
{
  Rcpp::Environment base = Rcpp::Environment::namespace_env("base");
  Rcpp::Function f = base["%*%"];
  return f(A, B);
}


// invert matrix
//[[Rcpp::export]]
NumericMatrix inv(NumericMatrix A)
{
  Rcpp::Environment base = Rcpp::Environment::namespace_env("base");
  Rcpp::Function f = base["solve"];
  return f(A);
}


// cast vector to one-column matrix
//[[Rcpp::export]]
NumericMatrix mat2vec(NumericVector v)
{
  int n = v.length();
  NumericMatrix m(n, 1);
  for(int i = 0; i < n; ++i)
  {
    m(i,0) = v[i];
  }
  return m;
}
