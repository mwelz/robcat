## hessian matrix of p_X

px_d2_mu2 <- function(x, mu, sigma2)
{
  pxtheta <- px(x = x, mu = mu, sigma = sqrt(sigma2))
  pxprime_mu <- px_d1_mu(x = x, mu = mu, sigma2 = sigma2)
  out <- ((x - mu) * pxprime_mu - pxtheta) / sigma2
  if(is.nan(out)) out <- 0.0 ## account for 0 * Inf
  return(out)
}


px_d2_sigma22 <- function(x, mu, sigma2)
{
  pxtheta_scaled <- px(x = x, mu = mu, sigma = sqrt(sigma2)) / sigma2
  pxprime_sigma2 <- px_d1_sigma2(x = x, mu = mu, sigma2 = sigma2)
  z2 <- ((x - mu)^2) / sigma2
  out <- (pxprime_sigma2 * (z2 - 1) + pxtheta_scaled * (1 - 2*z2)) / (2*sigma2)
  if(is.nan(out)) out <- 0.0 ## account for 0 * Inf
  return(out)
}


# cross-derivative
px_d2_musigma2 <- function(x, mu, sigma2)
{
  pxtheta <- px(x = x, mu = mu, sigma = sqrt(sigma2))
  pxprime_sigma2 <- px_d1_sigma2(x = x, mu = mu, sigma2 = sigma2)
  scaling <- (x - mu) / sigma2
  out <- scaling * (pxprime_sigma2 - pxtheta / sigma2)
  if(is.nan(out)) out <- 0.0 ## account for 0 * Inf
  return(out)
}


## hessian matrix of marginal density wrt theta
px_hessian <- function(x, mu, sigma2, num_y)
{
  ## initialize
  d <- num_y + 2L
  hess <- matrix(0.0, nrow = d, ncol = d)
  
  ## fill up nonzero entries; rest are structural zeros
  hess[2L,2L] <-      px_d2_mu2(x = x, mu = mu, sigma2 = sigma2)
  hess[3L,3L] <-  px_d2_sigma22(x = x, mu = mu, sigma2 = sigma2)
  cross       <- px_d2_musigma2(x = x, mu = mu, sigma2 = sigma2)
  hess[2L,3L] <- hess[3L,2L] <- cross
  return(hess)
}


#### sanity check
# mu <- 0
# sigma2 <- 1.2
# rho <- 0.5
# thres <- c(-1.5, -1)
# y_set <- c(1L, 2L, 3L)
# num_y <- length(y_set)
# theta <- c(rho, mu, sigma2, thres)
# d <- length(theta)
# set.seed(1)
# x <- rnorm(n = 1, mean = mu, sd = sqrt(sigma2))
# tol <- 0.0001 
# 
# f <- function(theta) px(x = x, mu = theta[2], sigma = sqrt(theta[3]))
# apprx <- numDeriv::hessian(f, theta)
# exact <- px_hessian(x = x, mu = mu, sigma2 = sigma2, num_y = num_y)
# all.equal(exact, apprx, tolerance = tol) # TRUE
