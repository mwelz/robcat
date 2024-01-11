#include <Rcpp.h>
#include "otherfun.h"
using namespace Rcpp;


// [[Rcpp::export]]
double internalfun(double x) {
  return x * 2;
}

