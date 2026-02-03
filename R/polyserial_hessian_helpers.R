phi <- function(x) stats::dnorm(x = x, mean = 0, sd = 1)
phi_prime <- function(x) (-x) * phi(x)  # first derivative

get_taustar_y <- function(x, rho, mu, sigma2, tau_y)
{
  numer <- tau_y - rho * (x - mu) / sqrt(sigma2) 
  denom <- sqrt(1.0 - rho^2)
  taustar <- numer / denom 
  return(taustar)
}


## first order derivatives
tstar_d1_rho <- function(x, rho, mu, sigma2, tau_y)
{
  tstar <- get_taustar_y(
    x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = tau_y)
  rootscale <- 1.0 / sqrt(1 - rho^2)
  z <- (x - mu) / sqrt(sigma2)
  rootscale * (rho * tstar * rootscale - z) 
}


tstar_d1_mu <- function(x, rho, mu, sigma2)
{
  denom <- sqrt(sigma2 * (1 - rho^2))
  rho / denom
}


tstar_d1_sigma2 <- function(x, rho, mu, sigma2)
{
  numer <- rho * (x - mu)
  denom <- 2.0 * sigma2^(1.5) * sqrt(1.0 - rho^2)
  numer / denom
}

# k is k-th threshold
tstar_d1_tauk <- function(x, y, k, rho, mu, sigma2)
{
  if(y == k)
  {
    out <- 1.0 / sqrt(1.0 - rho^2)
  } else
  {
    out <- 0.0
  }
  out
}
