#include <Rcpp.h>
#include <cmath> // for std::log(); the natural logarithm
#include "polycor.h"
using namespace Rcpp;

// empirical probability mass function for a bivariate sample
// [[Rcpp::export]]
NumericVector fhat(NumericVector x, 
                   NumericVector y,
                   int Kx,
                   int Ky)
{
  
  // initialize
  int K = Kx * Ky;
  int n = x.length();
  NumericVector out(K);
  int k = 0;
  
  for(int i = 0; i < Kx; ++i)
  {
    for(int j = 0; j < Ky; ++j)
    {
      
      // given responses (accounts for 0-indexing)
      int s = i + 1;
      int r = j + 1;
      double count = Rcpp::sum(x == s & y == r); // double to avoid integer division
      out[k] = count / n;
      k++;
    } // FOR r
  } // FOR s
  
  return out;
}


// generate correlation matrix with correlation rho
// [[Rcpp::export]]
NumericMatrix make_cormat(double rho)
{ 
  NumericMatrix m(2);
  m(0,0) = 1.0;
  m(1,0) = rho;
  m(0,1) = rho;
  m(1,1) = 1.0;
  return m;
}

// helper function for C++ equivalent to any(x < 10)
// taken from https://gallery.rcpp.org/articles/sugar-any/
// [[Rcpp::export]]
bool any_sug(LogicalVector x){
  // Note the use of is_true to return a bool type
  return is_true(any(x == TRUE));
}

// Rcpp implementation of mvtnorm::pmvnorm()
// partially taken from https://stackoverflow.com/q/51290014
// [[Rcpp::export]]
double pmvnorm_cpp(NumericVector lower,
                   NumericVector upper,
                   NumericVector mean,
                   NumericMatrix corr)
{
  Rcpp::Environment mvtnorm = Rcpp::Environment::namespace_env("mvtnorm");
  Rcpp::Function f = mvtnorm["pmvnorm"];
  
  // NumericVector mean(lower.length());
  NumericVector prob(1);
  prob = f(lower, upper, mean, corr, Named("keepAttr", false)); 
  
  // it may happen that the computed probability is slightly negative
  // due to being computationally zero, so truncate to avoid issues 
  double out;
  if(prob[0] < 0.0)
  {
    out = 0.0;
  } else{
    out = prob[0];
  } // IF
  
  return out;
} // FUN


// p_k(theta), at given responses (x,y)
// the given responses (x,y) \in [Kx]x[Ky] are 1-based
// README: Whenever variances aren't unit, we MUST pass sigma (not corr) to mvtnorm::pmvnorm() [ipv corr].
// In the polychoric model, we can pass either because variances under this model are unit, but we pass corr
// because it's faster to compute.
// However, when variances aren't unit---like in a contaminating distribution---it's crucial that we pass sigma 
// and not corr. That's why we have separate functions for the model distribution and contaminating distribution
// [[Rcpp::export]]
double pk_theta(int x, int y,      
                NumericMatrix cormat, // correlation matrix; unit diagonal
                NumericVector thresX, // X-thresholds; first and last values are ∓∞
                NumericVector thresY, // Y-thresholds; first and last values are ∓∞
                NumericVector mean)   // = {0,0} at polycor model
{
  // define bounds of integral
  NumericVector lower = {thresX[x-1], thresY[y-1]};
  NumericVector upper = {thresX[x], thresY[y]};

  // calculate the probability p_k(theta)
  double prob = pmvnorm_cpp(lower, upper, mean, cormat);
  return prob;
}


// p_k(theta) for p_k a contaminating distribution, at given responses (x,y)
// the given responses (x,y) \in [Kx]x[Ky] are 1-based
// README: Whenever variances aren't unit, we MUST pass sigma (not corr) to mvtnorm::pmvnorm() [ipv corr].
// In the polychoric model, we can pass either because variances under this model are unit. 
// However, when variances aren't unit---like in a contaminating distribution---it's crucial that we pass sigma 
// and not corr. That's why we have separate functions for the model distribution and contaminating distribution
// [[Rcpp::export]]
double pk_theta_contam(int x, int y,      
                       NumericMatrix sigma,  // covariance matrix
                       NumericVector thresX, // X-thresholds; first and last values are ∓∞
                       NumericVector thresY, // Y-thresholds; first and last values are ∓∞
                       NumericVector mean)  
{
  // define bounds of integral
  NumericVector lower = {thresX[x-1], thresY[y-1]};
  NumericVector upper = {thresX[x], thresY[y]};
  
  // call mvtnorm::pmvnorm()
  Rcpp::Environment mvtnorm = Rcpp::Environment::namespace_env("mvtnorm");
  Rcpp::Function f = mvtnorm["pmvnorm"];
  
  NumericVector prob(1);
  prob = f(lower, upper, mean, Named("sigma", sigma), Named("keepAttr", false)); 
  
  // it may happen that the computed probability is slightly negative
  // due to being computationally zero, so truncate to avoid issues 
  double out;
  if(prob[0] < 0.0)
  {
    out = 0.0;
  } else{
    out = prob[0];
  } // IF
  
  return out;
  
}

// loss of Ruckstuhl + Welsh (2001, AoS), with slope parts predefined
// [[Rcpp::export]]
double rho_fun_fast(double x,
                    double c1,
                    double c2,
                    double logc1p1, // log(c1) + 1
                    double logc2p1) // log(c2) + 1
{
  double out;
  if(x == 0)
  {
    // case 1: x is 0, so return 0
    out = 0.0;
  } else
  {
    // case 2: x \neq 0
    if(x < c1)
    {
      out = x * logc1p1 - c1;
    } else if(x > c2)
    {
      out = x * logc2p1 - c2;
    } else
    {
      // occurs iff c1 <= x && x <= c2
      out = x * std::log(x); 
    } // IF
  } // IF

  return out;
} // FUN


// loss of Ruckstuhl + Welsh (2001, AoS)
// [[Rcpp::export]]
double rho_fun_cpp(double x,
                   double c1,
                   double c2)
{
  double out;
  double logc1p1 = std::log(c1) + 1.0;
  double logc2p1 = std::log(c2) + 1.0;
  out = rho_fun_fast(x, c1, c2, logc1p1, logc2p1);
  return out;
}

//[[Rcpp::export]]
double objective_cpp_fast(
    double rho, 
    NumericVector f,      // density function
    NumericVector thresX, // X-thresholds; first and last values are ∓∞
    NumericVector thresY, // Y-thresholds; first and last values are ∓∞
    double c1,
    double c2,
    int Kx,
    int Ky,
    int K, // = Kx * Ky
    double logc1p1, // = std::log(c1) + 1.0
    double logc2p1, // = std::log(c2) + 1.0
    NumericVector mean, // = {0,0}
    double maxcor) // = 0.999
{
  double OUT;
  NumericVector out(K);
  NumericVector probs(K);
  NumericMatrix cormat = make_cormat(rho);
  
  if(rho < maxcor * (-1) || rho > maxcor )
  {
    OUT = 100000.; // TODO: come up with better way for implicit constraints
  } else
  {
    // initialize
    int k = 0;
    
    // loop over responses
    for(int x = 1; x <= Kx; ++x)
    {
      for(int y = 1; y <= Ky; ++y)
      {
        // x and y are the given responses at (x,y)
        double prob_k = pk_theta(x, y, cormat, thresX, thresY, mean);
        double x = f[k] / prob_k; // pearson residual
        out[k] = rho_fun_fast(x, c1, c2, logc1p1, logc2p1) * prob_k;
        probs[k] = prob_k;
        k++;
      } // FOR y
    } // FOR x
    
    // penalize values of theta that lead to zero class probabilities
    if(any_sug(probs == 0))
    {
      // avoid a stepwise pattern, so whenever zero probs assumption is violated, return same large penalty value
      OUT = 100000.; // TODO: make better 
    } else
    {
      // assign zero to empirically empty classes
      out[f == 0] = 0;
      OUT = sum(out);
    } // IF penalization
  } // IF maxcor
  
  return OUT;
}


//[[Rcpp::export]]
double objective_cpp(
    double rho, 
    NumericVector f,      // density function
    NumericVector thresX, // X-thresholds; first and last values are ∓∞
    NumericVector thresY, // Y-thresholds; first and last values are ∓∞
    double c1,
    double c2,
    NumericVector mean, // = {0,0}
    double maxcor = 0.999)
{
  int Kx = thresX.length() - 1;
  int Ky = thresY.length() - 1;
  int K = Kx * Ky;
  double logc1p1 = std::log(c1) + 1.0;
  double logc2p1 = std::log(c2) + 1.0;

  double out = objective_cpp_fast(rho, f, 
                                  thresX, thresY,
                                  c1, c2,
                                  Kx, Ky, K,
                                  logc1p1, logc2p1, 
                                  mean,
                                  maxcor);
  return out;
}

// MLE loss as in paper
//[[Rcpp::export]]
double objective_mle_cpp(
  double rho,
  NumericVector freq,   // frequency vector
  NumericVector thresX, // X-thresholds; first and last values are ∓∞
  NumericVector thresY, // Y-thresholds; first and last values are ∓∞
  int Kx,
  int Ky,
  int K, // = Kx * Ky
  NumericVector mean) // = {0,0}
{
  double OUT;
  NumericVector out(K);
  NumericVector probs(K);
  NumericMatrix cormat = make_cormat(rho);
  
  // initialize
  int k = 0;
  
  // loop over responses
  for(int x = 1; x <= Kx; ++x)
  {
    for(int y = 1; y <= Ky; ++y)
    {
      // x and y are the given responses at (x,y)
      double prob_k = pk_theta(x, y, cormat, thresX, thresY, mean);
      out[k] = freq[k] * std::log(prob_k);
      probs[k] = prob_k;
      k++;
    } // FOR y
  } // FOR x
  
  
  // penalize values of rho that lead to zero class probabilities
  if(any_sug(probs == 0))
  {
    // avoid a stepwise pattern, so whenever zero probs assumption is violated, return same large penalty value
    OUT = 100000.; // TODO: make better 
  } else
  {
    // assign zero to empirically empty classes
    out[freq == 0] = 0;
    OUT = -sum(out); // make negative to make it a minimization problem
  } // IF penalization
  
  return OUT;
}



// probability mass function under contaminated model
//[[Rcpp::export]]
NumericVector feps_cpp(
    NumericVector thresX,
    NumericVector thresY,
    double eps, // contamination fraction
    NumericMatrix covmat_true, // = cormat_true since polycor model has unit variance
    NumericMatrix covmat_contam, 
    NumericVector mean_true,
    NumericVector mean_contam)
{
  // initialize
  int Kx = thresX.length() - 1;
  int Ky = thresY.length() - 1;
  int K  = Kx * Ky;
  NumericVector out(K);
  double meps = 1.0 - eps;
  int k = 0;
  
  // loop over the given responses
  for(int x = 1; x <= Kx; ++x)
  {
    for(int y = 1; y <= Ky; ++y)
    {
      // i and j are the given responses at (i,j)
      double fk_true   = pk_theta(x, y, covmat_true,   thresX, thresY, mean_true);
      double fk_contam = pk_theta_contam(x, y, covmat_contam, thresX, thresY, mean_contam);
      out[k] = meps * fk_true + eps * fk_contam;
      k++;
    } // FOR x
  } // FOR y
  
  return out;
}

// calculate theoretical probabilities
//[[Rcpp::export]]
NumericVector model_probabilities(
    double rho,
    NumericVector thresX,
    NumericVector thresY,
    int Kx, int Ky)
{
  
  NumericMatrix cormat = make_cormat(rho);
  NumericVector mean = {0., 0.};
  NumericVector out(Kx * Ky);
  int k = 0;
  
  for(int x = 1; x <= Kx; ++x)
  {
    for(int y = 1; y <= Ky; ++y)
    {
      // x and y are the given responses at (x,y)
      double prob_k = pk_theta(x, y, cormat, thresX, thresY, mean);
      out[k] = prob_k;
      k++;
    }
  }
  return out;
}

