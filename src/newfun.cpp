#include <Rcpp.h>
#include <cmath>
#include "otherfun.h"
using namespace Rcpp;



//' Multiplies two doubles
//'
//' @param v1 First value
//' @param v2 Second value
//' @return Product of v1 and v2
//' @export
 // [[Rcpp::export]]
double mult( double v1, double v2 ) {return v1*v2;}


// helper function, shall not be exported to user: note the absence pof @export
// [[Rcpp::export]]
double add(double x, double y){return x+y;}


// [[Rcpp::export]]
double callothercpp(double x)
{
	double foo = std::log(internalfun(x));
	return foo;	
}

