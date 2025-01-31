#include <Rcpp.h>
using namespace Rcpp;

double pmvnorm_cpp(NumericVector,
                   NumericVector,
                   NumericVector,
                   NumericMatrix);

double pk_theta(int, int,      
                NumericMatrix, // correlation matrix; unit diagonal
                NumericVector, // X-thresholds; first and last values are ∓∞
                NumericVector, // Y-thresholds; first and last values are ∓∞
                NumericVector);

NumericMatrix make_cormat(double);
  
