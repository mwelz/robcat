#include <Rcpp.h>
#include <math.h> // for std::sqrt() and pow()
#include <limits.h> // for infinity
#include <vector> // for vector class
#include "polycor.h" // import relevant helper functions
#include "matrixalgebra.h" // same
#include "polycor_variance.h"
using namespace Rcpp;


// [[Rcpp::export]]
bool is_infinite_pos(double x)
{
  double Inf = std::numeric_limits<double>::infinity();
  bool out;
  if(x == Inf)
  {
    out = true;
  } else{
    out = false;
  }
  return out;
}


// [[Rcpp::export]]
bool is_infinite_neg(double x)
{
  double Inf = std::numeric_limits<double>::infinity();
  bool out;
  if(x == -Inf)
  {
    out = true;
  } else{
    out = false;
  }
  return out;
}


// [[Rcpp::export]]
bool is_infinite(double x)
{
  double Inf = std::numeric_limits<double>::infinity();
  bool out;
  if(x == Inf || x  == -Inf)
  {
    out = true;
  } else{
    out = false;
  }
  return out;
}


// stats::dnorm() at standard normal
// [[Rcpp::export]]
double dnorm(double x)
{
  Rcpp::Environment stats = Rcpp::Environment::namespace_env("stats");
  Rcpp::Function f = stats["dnorm"];
  NumericVector d = f(x);
  return d[0];
}


// stats::pnorm() at standard normal
// [[Rcpp::export]]
double pnorm(double x)
{
  Rcpp::Environment stats = Rcpp::Environment::namespace_env("stats");
  Rcpp::Function f = stats["pnorm"];
  NumericVector p = f(x);
  return p[0];
}


// bivariate standard gaussian density
// [[Rcpp::export]]
double dnorm2(double x, double y, double rho)
{
  double out;
  
  if(is_infinite(x) || is_infinite(y))
  {
    out = 0.0;
  } else
  {
    Rcpp::Environment mvtnorm = Rcpp::Environment::namespace_env("mvtnorm");
    Rcpp::Function f = mvtnorm["dmvnorm"];
    NumericMatrix cormat = make_cormat(rho);
    NumericVector mean(2);
    NumericVector eval = NumericVector::create(x, y);
    NumericVector d = f(eval, mean, cormat);
    out = d[0];
  }
  return out;
}


// derivative of F(theta) = \Phi_2(theta) w.r.t. rho; equals phi_2(x,y,;rho)
// [[Rcpp::export]]
double F_prime_rho(double x,
                   double y,
                   double rho)
{
  double out = dnorm2(x, y, rho);
  return out;
}


// derivative of F(theta) = \Phi_2(theta) w.r.t. a single threshold
// eq. (12) in Olsson (1979, Psychometrika)
// 'thres_variable' is the threshold that we  partially differentiate w.r.t.
// 'thres_fixed' is the threshold of the other dimension (constant in the partial differentiation)
// [[Rcpp::export]]
double F_prime_thres(double thres_variable,
                     double thres_fixed,
                     double rho)
{
  double out;
  if(is_infinite(thres_variable))
  {
    // thres_variable is infinite
    out = 0.0;
  } else if(is_infinite_neg(thres_fixed))
  {
    // thres_variable is finite, but fixed is -inf
    out = 0.0;
  } else if(is_infinite_pos(thres_fixed))
  {
    // thres_variable is finite, but fixed is +inf
    out = dnorm(thres_variable);
  } else
  {
    double phi    = dnorm(thres_variable);
    double rhofun = std::sqrt(1.0 - rho * rho);
    double arg    = thres_fixed - rho * thres_variable / rhofun;
    double p      = pnorm(arg);
    out           = phi * p;
  }
  
  return out;
}
  
  
// x and y are given responses in [Kx]x[Ky]
// the thres variables; first and last values are ∓∞
// eq. 10 in Olsson (1979, Psychometrika)
// returns gradient of pk(theta) at (x,y) if we differentiate with respect to all thresholds in thresX
// [[Rcpp::export]]
NumericVector pk_prime_thresX(
    int x, int y,
    NumericVector thresX,
    NumericVector thresY,
    double rho,
    int Kx)
{
  NumericVector out(Kx-1);
  int xminus1 = x - 1;
  int yminus1 = y - 1;
  
  // first and last threshold values are constant (∓∞), so skip them
  for(int k = 1; k < Kx; ++k)
  {
    if(k == x)
    {
      double deriv1 = F_prime_thres(thresX[k], thresY[y], rho);
      double deriv2 = F_prime_thres(thresX[k], thresY[yminus1], rho);
      out[k-1] = deriv1 - deriv2;
    } else if (k == xminus1)
    {
      double deriv1 = F_prime_thres(thresX[k], thresY[y], rho);
      double deriv2 = F_prime_thres(thresX[k], thresY[yminus1], rho);
      out[k-1] = deriv2 - deriv1;
    } else
    {
      out[k-1] = 0.0;
    }
  } // IF
  return out;
}


// x and y are given responses in [Kx]x[Ky]
// the thres variables; first and last values are ∓∞
// adapted from eq. 10 in Olsson (1979, Psychometrika)
// returns gradient of pk(theta) at (x,y) if we differentiate with respect to all finite thresholds in thresY
// [[Rcpp::export]]
NumericVector pk_prime_thresY(
    int x, int y,
    NumericVector thresX,
    NumericVector thresY,
    double rho,
    int Ky)
{
  NumericVector out(Ky-1);
  int xminus1 = x - 1;
  int yminus1 = y - 1;
  
  // first and last threshold values are constant (∓∞), so skip them
  for(int k = 1; k < Ky; ++k)
  {
    if(k == y)
    {
      double deriv1 = F_prime_thres(thresY[k], thresX[x], rho);
      double deriv2 = F_prime_thres(thresY[k], thresX[xminus1], rho);
      out[k-1] = deriv1 - deriv2;
    } else if (k == yminus1)
    {
      double deriv1 = F_prime_thres(thresY[k], thresX[x], rho);
      double deriv2 = F_prime_thres(thresY[k], thresX[xminus1], rho);
      out[k-1] = deriv2 - deriv1;
    } else
    {
      out[k-1] = 0.0;
    }
  } // IF
  return out;
}


// x and y are given responses in [Kx]x[Ky]
// the thres variables; first and last values are ∓∞
// eq. 8 in Olsson (1979, Psychometrika)
// returns gradient of pk(theta) at (x,y) if we differentiate with respect to rho
// [[Rcpp::export]]
double pk_prime_rho(
    int x, int y,
    NumericVector thresX,
    NumericVector thresY,
    double rho)
{
  double sum1 = F_prime_rho(thresX[x], thresY[y], rho);
  double sum2 = F_prime_rho(thresX[x-1], thresY[y], rho);
  double sum3 = F_prime_rho(thresX[x], thresY[y-1], rho);
  double sum4 = F_prime_rho(thresX[x-1], thresY[y-1], rho);
  double out = sum1 - sum2 - sum3 + sum4;
  return out;
}


// x and y are given responses in [Kx]x[Ky]
// the thres variables; first and last values are ∓∞
// returns gradient of pk(theta) at (x,y) if we differentiate with respect to theta
// [[Rcpp::export]]
NumericVector pk_prime_theta(
    int x, int y,
    NumericVector thresX,
    NumericVector thresY,
    double rho,
    int num_param, // = Kx + Ky - 1
    int Kx, int Ky)
{
  
  // compute gradients of pk(theta) w.r.t. to the variables in theta
  double grad_rho = pk_prime_rho(x, y, thresX, thresY, rho);
  NumericVector grad_thresX = pk_prime_thresX(x, y, thresX, thresY, rho, Kx);
  NumericVector grad_thresY = pk_prime_thresY(x, y, thresX, thresY, rho, Ky);
  
  // fill up the final gradient
  NumericVector out(num_param); 
  out[0] = grad_rho;
  
  int ct = 0;
  for(int k = 1; k < Kx; ++k)
  {
    out[k] = grad_thresX[ct];
    ct++;
  }
  
  ct = 0;
  for(int k = Kx; k < num_param; ++k)
  {
    out[k] = grad_thresY[ct];
    ct++;
  }
  return out;
}


// first derivative of univariate standard gaussian density 
// [[Rcpp::export]]
double dnorm1_prime(double x)
{
  double out;
  if(is_infinite(x))
  {
    out = 0.0;
  } else {
    const double pi = 3.14159265358979323846;
    double denom = std::sqrt(2.0 * pi);
    double a = -x / denom;
    double b = std::exp(-x*x / 2.0);
    out = a * b;
  }
  return out;
}


// first derivative of bivariate standard gaussian density w.r.t. rho
// [[Rcpp::export]]
double dnorm2_prime_rho(
    double x, double y, double rho)
{
  double out;
  
  if(is_infinite(x) || is_infinite(y))
  {
    out = 0.0;
  } else
  {
    double mrho2 = 1. - rho * rho;
    double sum1  = mrho2 * (rho + x * y); 
    double sum2 = rho * (x*x - 2.*rho*x*y + y*y);
    double phi = dnorm2(x, y, rho);
    out = phi * (sum1 - sum2) / (mrho2 * mrho2);
  }
  return out;
}


// 2nd derivative of bivariate standard gaussian CDF w.r.t. first argument
// \partial^2{\Phi_2(x,y)} / \partial{x^2}
// 'variable' is first argument that we differentiate w.r.t.
// 'fixed' is second argument (constant in differentiation)
// [[Rcpp::export]]
double pnorm2_prime2_x_x(double variable, double fixed, double rho)
{
  double out;
  
  if(is_infinite(variable))
  {
    out = 0.0;
  } else if(is_infinite_neg(fixed))
  {
    // if variable is finite, but fixed = -inf
    out = 0.0;
  } else if(is_infinite_pos(fixed))
  {
    // if variable is finite, but fixed = +inf
    out = dnorm1_prime(variable);
  } else
  {
    double sqrtrho = std::sqrt(1. - rho * rho);
    double arg = (fixed - rho * variable) / sqrtrho;
    double phiprime = dnorm1_prime(variable);
    double Phi = pnorm(arg);
    double phiarg = dnorm(arg);
    double phi = dnorm(variable);
    out = phiprime * Phi - rho * phi * phiarg / sqrtrho;
  }
  return out;
}


// 2nd derivative of bivariate standard gaussian CDF w.r.t. first argument and rho
// 'variable' is first argument that we differentiate w.r.t.
// 'fixed' is second argument (constant in differentiation)
// [[Rcpp::export]]
double pnorm2_prime2_x_rho(
    double variable, double fixed, double rho)
{
  double out;
  
  if(is_infinite(variable) || is_infinite(fixed))
  {
    out = 0.0;
  } else
  {
    double rhom = 1. - rho * rho;
    double sqrtrho = std::sqrt(rhom);
    double rhom15 = std::pow(rhom, 1.5);
    double phiarg = dnorm((fixed - rho * variable) / sqrtrho);
    double phi = dnorm(variable);
    out = phi* phiarg * (rho * fixed - variable) / rhom15; 
  }
  return out;
}


// 2nd derivative of bivariate standard gaussian CDF w.r.t. 1st argument and 2nd argument
// [[Rcpp::export]]
double pnorm2_prime2_x_y(
    double x, double y, double rho)
{
  double out;
  
  if(is_infinite(x) || is_infinite(y))
  {
    out = 0.0;
  } else
  {
    double sqrtrho = std::sqrt(1. - rho * rho);
    double phi = dnorm(x);
    double phiarg = dnorm((y - rho * x) / sqrtrho);
    out =  phi* phiarg / sqrtrho;
  }
  return out;
}


// x and y are given responses in [Kx]x[Ky]
// the thres variables; first and last values are ∓∞
// returns 2nd derivative of pk(theta) at (x,y) if we differentiate twice w.r.t. rho
// [[Rcpp::export]]
double pk_prime2_rho2(
    int x, int y,
    NumericVector thresX,
    NumericVector thresY,
    double rho)
{
  double sum1 = dnorm2_prime_rho(thresX[x], thresY[y], rho);
  double sum2 = dnorm2_prime_rho(thresX[x-1], thresY[y], rho);
  double sum3 = dnorm2_prime_rho(thresX[x], thresY[y-1], rho);
  double sum4 = dnorm2_prime_rho(thresX[x-1], thresY[y-1], rho);
  double out = sum1 - sum2 - sum3 + sum4;
  return out;
}


// x and y are given responses in [Kx]x[Ky]
// the thres variables; first and last values are ∓∞
// returns 2nd derivative of pk(theta) at (x,y) if we differentiate with respect to rho and then threshold a_k
// [[Rcpp::export]]
double pk_prime2_thresX_rho(
    int x, int y, int k,
    NumericVector thresX,
    NumericVector thresY,
    double rho)
{
  double out;
  
  if(k == x)
  {
    double deriv1 = pnorm2_prime2_x_rho(thresX[k], thresY[y], rho);
    double deriv2 = pnorm2_prime2_x_rho(thresX[k], thresY[y-1], rho);
    out = deriv1 - deriv2;
  } else if (k == x - 1)
  {
    double deriv1 = pnorm2_prime2_x_rho(thresX[k], thresY[y], rho);
    double deriv2 = pnorm2_prime2_x_rho(thresX[k], thresY[y-1], rho);
    out = deriv2 - deriv1;
  } else
  {
    out = 0.0;
  }
  return out;
}


// x and y are given responses in [Kx]x[Ky]
// the thres variables; first and last values are ∓∞
// returns 2nd derivative of pk(theta) at (x,y) if we differentiate with respect to rho and then threshold b_k
// [[Rcpp::export]]
double pk_prime2_thresY_rho(
    int x, int y, int k,
    NumericVector thresX,
    NumericVector thresY,
    double rho)
{
  double out;
  
  if(k == y)
  {
    double deriv1 = pnorm2_prime2_x_rho(thresY[k], thresX[x], rho);
    double deriv2 = pnorm2_prime2_x_rho(thresY[k], thresX[x-1], rho);
    out = deriv1 - deriv2;
  } else if (k == y - 1)
  {
    double deriv1 = pnorm2_prime2_x_rho(thresY[k], thresX[x], rho);
    double deriv2 = pnorm2_prime2_x_rho(thresY[k], thresX[x-1], rho);
    out = deriv2 - deriv1;
  } else
  {
    out = 0.0;
  }
  return out;
}


// x and y are given responses in [Kx]x[Ky]
// the thres variables; first and last values are ∓∞
// returns 2nd derivative of pk(theta) at (x,y) if we differentiate twice with respect to threshold a_k
// [[Rcpp::export]]
double pk_prime2_thresX2(
    int x, int y, int k,
    NumericVector thresX,
    NumericVector thresY,
    double rho)
{
  double out;
  
  if(k == x)
  {
    double deriv1 = pnorm2_prime2_x_x(thresX[k], thresY[y], rho);
    double deriv2 = pnorm2_prime2_x_x(thresX[k], thresY[y-1], rho);
    out = deriv1 - deriv2;
  } else if (k == x - 1)
  {
    double deriv1 = pnorm2_prime2_x_x(thresX[k], thresY[y], rho);
    double deriv2 = pnorm2_prime2_x_x(thresX[k], thresY[y-1], rho);
    out = deriv2 - deriv1;
  } else
  {
    out = 0.0;
  }
  return out;
}


// x and y are given responses in [Kx]x[Ky]
// the thres variables; first and last values are ∓∞
// returns 2nd derivative of pk(theta) at (x,y) if we differentiate twice with respect to threshold b_k
// [[Rcpp::export]]
double pk_prime2_thresY2(
    int x, int y, int k,
    NumericVector thresX,
    NumericVector thresY,
    double rho)
{
  double out;
  
  if(k == y)
  {
    double deriv1 = pnorm2_prime2_x_x(thresY[k], thresX[x], rho);
    double deriv2 = pnorm2_prime2_x_x(thresY[k], thresX[x-1], rho);
    out = deriv1 - deriv2;
  } else if (k == y - 1)
  {
    double deriv1 = pnorm2_prime2_x_x(thresY[k], thresX[x], rho);
    double deriv2 = pnorm2_prime2_x_x(thresY[k], thresX[x-1], rho);
    out = deriv2 - deriv1;
  } else
  {
    out = 0.0;
  }
  return out;
}


// x and y are given responses in [Kx]x[Ky]
// the thres variables; first and last values are ∓∞
// returns 2nd derivative of pk(theta) at (x,y) if we differentiate with respect to threshold b_k and b_l
// [[Rcpp::export]]
double pk_prime2_thresY_kl(
    int x, int y, int k, int l,
    NumericVector thresX,
    NumericVector thresY,
    double rho)
{
  double out;
  if(k == l)
  {
    out = pk_prime2_thresY2(x, y, k, thresX, thresY, rho);
  } else{
    out = 0.0;
  }
  return out; 
}


// x and y are given responses in [Kx]x[Ky]
// the thres variables; first and last values are ∓∞
// returns 2nd derivative of pk(theta) at (x,y) if we differentiate with respect to threshold a_k and a_l
// [[Rcpp::export]]
double pk_prime2_thresX_kl(
    int x, int y, int k, int l,
    NumericVector thresX,
    NumericVector thresY,
    double rho)
{
  double out;
  if(k == l)
  {
    out = pk_prime2_thresX2(x, y, k, thresX, thresY, rho);
  } else{
    out = 0.0;
  }
  return out; 
}


// x and y are given responses in [Kx]x[Ky]
// the thres variables; first and last values are ∓∞
// returns 2nd derivative of pk(theta) at (x,y) if we differentiate with respect to threshold a_k and b_l
// no reverse version needed thanks to Young's theorem
// [[Rcpp::export]]
double pk_prime2_thresXk_thresYl(
    int x, int y, int k, int l,
    NumericVector thresX,
    NumericVector thresY,
    double rho)
{
  int xm1 = x - 1;
  int ym1 = y - 1;
  double out;
  
  if((k == x && l == y) || (k == xm1 && l == ym1))
  {
    out = pnorm2_prime2_x_y(thresX[k], thresY[l], rho);
    
  } else if((k == xm1 && l == y) || (k == x && l == ym1))
  {
    out = -pnorm2_prime2_x_y(thresX[k], thresY[l], rho);
  } else{
    out = 0.0;
  }
  return out;
}


// index set (0-based) in theta that is associated with a_1,...,a_{Kx-1}
// [[Rcpp::export]]
std::vector<int> get_indices_thresX(int Kx)
{
  std::vector<int> out(Kx-1);
  for(int i = 0; i < Kx; i++)
  {
    out[i] = i + 1;
  }
  return out;
}


// index set (0-based) in theta that is associated with b_1,...,b_{Ky-1}
// [[Rcpp::export]]
std::vector<int> get_indices_thresY(int Kx, int Ky)
{
  std::vector<int> out(Ky-1);
  int ct = 0;
  for(int i = Kx; i < Kx+Ky-1; i++)
  {
    out[ct] = i;
    ct++;
  }
  return out;
}


// index set (0-based) in theta that is associated with rho
// [[Rcpp::export]]
std::vector<int> get_indices_rho()
{
  std::vector<int> out(1);
  out[0] = 0;
  return out;
}


// check if index i (0-based) is in index set v
// [[Rcpp::export]]
bool is_in(int i, std::vector<int> v)
{
  bool out = false;
  int n = v.size();
  
  for(int j = 0; j < n; j++)
  {
    if(i == v[j])
    {
      out = true;
      break;
    } else{
      continue;
    }
  }
  return out;
}


// these two functions are intended for the for loop in pk_prime2_theta2()
// i is an index in theta, and they transform i to the a_k or b_l space
// [[Rcpp::export]]
int to_a(int i)
{
  return i;
}

int to_b(int i, int Kx)
{
  return i - Kx + 1;
}


// x and y are given responses in [Kx]x[Ky]
// the thres variables; first and last values are ∓∞
// returns Hessian of pk(theta) at (x,y) if we differentiate twice with respect to theta
// [[Rcpp::export]]
NumericMatrix pk_prime2_theta2(
    int x, int y,
    NumericVector thresX,
    NumericVector thresY,
    double rho,
    int Kx, int Ky)
{
  // dimension of theta
  int d = Kx + Ky - 1;
  
  // index sets
  std::vector<int> Irho = get_indices_rho(); // rho (=0)
  std::vector<int> Ia = get_indices_thresX(Kx); // x-thresholds
  std::vector<int> Ib = get_indices_thresY(Kx, Ky); // y-thresholds
  
  // init
  NumericMatrix m(d);
  double z; // temp var
  int k; // counter for a_k or b_k (1-based)
  int l; // same
  
  for(int i = 0; i < d; i++)
  {
    for(int j = 0; j < d; j++)
    {
      if(i <= j)
      {
        if(is_in(i,Ia) && is_in(j, Ib))
        {
          k = to_a(i);     // i counts the a_k, k = 1,...,Kx-1
          l = to_b(j, Kx); // j counts the b_l, l = 1,...,Ky-1
          z = pk_prime2_thresXk_thresYl(x, y, k, l, thresX, thresY, rho);
          m(i,j) = z;
          m(j,i) = z;
        } else if(is_in(i,Ib) && is_in(j,Ia))
        {
          k = to_a(j);     // j counts the a_k, k = 1,...,Kx-1
          l = to_b(i, Kx); // i counts the b_l, l = 1,...,Ky-1
          z = pk_prime2_thresXk_thresYl(x, y, k, l, thresX, thresY, rho);
          m(i,j) = z;
          m(j,i) = z;
        } else if(is_in(i, Ia) && is_in(j,Ia))
        {
          k = to_a(i); // i,j count the a_k, k = 1,...,Kx-1
          l = to_a(j); 
          z = pk_prime2_thresX_kl(x, y, k, l, thresX, thresY, rho);
          m(i,j) = z;
          m(j,i) = z;
        } else if(is_in(i,Ib) && is_in(j,Ib))
        {
          k = to_b(i, Kx); // i,j count the b_k, k = 1,...,Ky-1
          l = to_b(j, Kx); 
          z = pk_prime2_thresY_kl(x, y, k, l, thresX, thresY, rho);
          m(i,j) = z;
          m(j,i) = z;
        } else if(is_in(i,Irho) && is_in(j, Ia))
        {
          k = to_a(j);
          z = pk_prime2_thresX_rho(x, y, k, thresX, thresY, rho);
          m(i,j) = z;
          m(j,i) = z;
        } else if(is_in(i,Ia) && is_in(j, Irho))
        {
          k = to_a(i);
          z = pk_prime2_thresX_rho(x, y, k, thresX, thresY, rho);
          m(i,j) = z;
          m(j,i) = z;
        } else if(is_in(i,Irho) && is_in(j,Ib))
        {
          k = to_b(j, Kx);
          z = pk_prime2_thresY_rho(x, y, k, thresX, thresY, rho);
          m(i,j) = z;
          m(j,i) = z;
        } else if(is_in(i,Ib) && is_in(j,Irho))
        {
          k = to_b(i, Kx);
          z = pk_prime2_thresY_rho(x, y, k, thresX, thresY, rho);
          m(i,j) = z;
          m(j,i) = z;
        } else if(is_in(i,Irho) && is_in(j,Irho))
        {
          z = pk_prime2_rho2(x, y, thresX, thresY, rho);
          m(i,i) = z; // this case only occurs when i=j=0
        } else{
          continue;
        } // if i <= j

      } else 
      {
        continue;
      }// if
    } // for
  } // for
  
  return m;
}


// x and y are given responses in [Kx]x[Ky]
// the thres variables; first and last values are ∓∞
// returns sk(theta) at (x,y)
// [[Rcpp::export]]
NumericVector sk_theta(
    int x, int y,
    NumericVector thresX,
    NumericVector thresY,
    double rho,
    int Kx, int Ky
)
{
  NumericMatrix cormat   = make_cormat(rho);
  NumericVector mean     = {0., 0.};
  int d                  = Kx + Ky - 1; // number of parameters
  double pk              = pk_theta(x, y, cormat, thresX, thresY, mean); 
  NumericVector pk_prime = pk_prime_theta(x, y, thresX, thresY, rho, d, Kx, Ky);
  return pk_prime / pk;
}


// [[Rcpp::export]]
double w_fun(double x, double c1, double c2)
{
  double out;
  if(x > c2)
  {
    out = c2 / x;
  } else if(x < c1)
  {
    out = c1 / x;
  } else
  {
    out = 1.0;
  }
  return out;
}


// [[Rcpp::export]]
double w_fun_prime(double x, double c1, double c2)
{
  double out;
  if(x > c2)
  {
    out = -c2 / (x*x);
  } else if(x < c1)
  {
    out = -c1 / (x*x);
  } else
  {
    out = 0.0;
  }
  return out;
}


// [[Rcpp::export]]
double in_interval(double x, double c1, double c2)
{
  double out;
  if(c1 <= x && x <= c2)
  {
    out = 1.0;
  } else
  {
    out = 0.0;
  }
  return out;
}


//[[Rcpp::export]]
List get_MW(
    double rho,
    NumericVector thresX,
    NumericVector thresY,
    NumericVector f, // density 
    double c1, double c2,
    int Kx, int Ky)
{
  // initialize
  int k = 0;
  int d = Kx + Ky - 1; // dim(theta)
  int K = Kx * Ky;
  double pearson, pk;
  NumericMatrix M(d);
  NumericMatrix W = no_init_matrix(d, K);
  
  // correlation matrix
  NumericMatrix cormat = make_cormat(rho);
  NumericVector mean = {0., 0.};
  
  // loop over all answer categories
  for(int x = 1; x <= Kx; x++)
  {
    for(int y = 1; y <= Ky; y++)
    { 
      // calculate sk(theta)
      pk = pk_theta(x, y, cormat, thresX, thresY, mean);
      NumericVector pk_prime = pk_prime_theta(x, y, thresX, thresY, rho, d, Kx, Ky);
      NumericVector sk = pk_prime / pk;
      
      // outer product of sk
      NumericMatrix ss = outer_vec(sk);
      
      // hessian of pk(theta)
      NumericMatrix hessian = pk_prime2_theta2(x, y, thresX, thresY, rho, Kx, Ky);
      
      // Qk(theta)
      NumericMatrix tmp = hessian / pk;
      NumericMatrix mss = ss * (-1.);
      NumericMatrix Qk = matplus(tmp, mss);
      
      // pearson residual
      pearson = f[k] / pk;
      
      // summands
      NumericMatrix sum1 = w_fun_prime(pearson, c1, c2) * pearson * ss;
      NumericMatrix sum2 = (-1.) * w_fun(pearson, c1, c2) * Qk;
      NumericMatrix sumk = f[k] * matplus(sum1, sum2);
      
      // increment to M
      M = matplus(M, sumk);
      
      // calculate W
      W(_, k) = sk * in_interval(pearson, c1, c2);
      
      // increment
      k++;
    }
  }
  
  return List::create(Rcpp::Named("M") = M, Rcpp::Named("W") = W);
}



// calculate fisher information matrix
//[[Rcpp::export]]
NumericMatrix get_fisher(
    double rho,
    NumericVector thresX,
    NumericVector thresY,
    int Kx, int Ky)
{
  // initialize
  int k = 0;
  int d = Kx + Ky - 1; // dim(theta)
  double pk;
  NumericMatrix J(d);

  // correlation matrix
  NumericMatrix cormat = make_cormat(rho);
  NumericVector mean = {0., 0.};
  
  // loop over all answer categories
  for(int x = 1; x <= Kx; x++)
  {
    for(int y = 1; y <= Ky; y++)
    { 
      // calculate sk(theta)
      pk = pk_theta(x, y, cormat, thresX, thresY, mean);
      NumericVector pk_prime = pk_prime_theta(x, y, thresX, thresY, rho, d, Kx, Ky);
      NumericVector sk = pk_prime / pk;
      
      // outer product of sk
      NumericMatrix ss = outer_vec(sk);
      
      // hessian of pk(theta)
      NumericMatrix hessian = pk_prime2_theta2(x, y, thresX, thresY, rho, Kx, Ky);
      
      // Qk(theta)
      NumericMatrix tmp = hessian / pk;
      NumericMatrix mss = ss * (-1.);
      NumericMatrix Qk = matplus(tmp, mss);
      
      // summand
      NumericMatrix sumk = Qk * pk;
      
      // increment matrix
      J = matplus(J, sumk);
      
      // increment
      k++;
    }
  }
  
  return J * (-1.);
}


////////////// GRADIENT COMPUTATION:

// RAF of Ruckstuhl + Welsh (2001, AoS)
// [[Rcpp::export]]
double raf_fast(double x, // PR (1-based)
                double c1, // c's are 1-based
                double c2)
{
  double out;
  
  if(x < c1)
  {
    out = c1;
  } else if(x > c2)
  {
    out = c2;
  } else
  {
    // occurs iff c1 <= x && x <= c2
    out = x;
  }
  
  return out;
} // FUN



// the thres variables; first and last values are ∓∞
// returns gradient of loss differentiated with respect to theta
// [[Rcpp::export]]
NumericVector gradient_loss_cpp(
    NumericVector f,      // density function
    double rho,
    NumericVector thresX, // X-thresholds; first and last values are ∓∞
    NumericVector thresY, // Y-thresholds; first and last values are ∓∞
    NumericVector mean, // = {0,0}
    int Kx, int Ky,
    int num_params, // = Kx + Ky -1
    double c1, // c's are 1-based
    double c2)
{
  NumericVector OUT(num_params);
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
      double pr = f[k] / prob_k;           // pearson residual (1-based)
      double raf_k = raf_fast(pr, c1, c2); // residual adjustment function
      
      // gradient of polychoric model w.r.t. theta at (x,y)
      NumericVector grad_pk = pk_prime_theta(x, y, thresX, thresY, rho, num_params, Kx, Ky);
      
      // update output
      OUT = OUT + grad_pk * raf_k;
      k++;
    } 
  }
  
  return OUT * (-1.);
} // FUN



// be careful when copy/pasting functions with variable/fixed inputs; these need to be flipped when differential changes
