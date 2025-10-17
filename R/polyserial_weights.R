weights_mainfun <- function(x, y, alpha, rho, mu, sigma, tau, num_y)
{
  marg <- px(x = x, mu = mu, sigma = sigma)
  cond <- py_x_vectorized(y = y, x = x, rho = rho,
                          mu = mu, sigma = sigma, thres = tau,num_y = num_y)
  pxy <- unname(marg * cond)
  wgt <- pxy^alpha
  return(wgt)
}


## maximum theoretically possible weight
maxweight <- function(theta, alpha)
{
  d     <- length(theta)
  rho   <- theta[1L]
  mu    <- theta[2L]
  sigma <- sqrt(theta[3L])
  tau   <- theta[4L:d]
  num_y <- length(tau) + 1L
  maxvals <- rep(NA_real_, num_y)
  
  for(y in seq_len(num_y))
  {
    fn <- function(x)
    {
      wgt <- weights_mainfun(
        x = x, y = y, alpha = alpha, rho = rho, 
        mu = mu, sigma = sigma, tau = tau, num_y = num_y)
      return(-wgt) # we maximize, so minimize negative objective
    }
    
    opt_y <- stats::optim(mu, fn = fn, method = "BFGS") # max should lie close to location param
    maxval_y <- (-1) * opt_y$value # recast to maximization problem
    maxvals[y] <- maxval_y
  }
  return(max(maxvals))
}

## compute rescaled weights
weights <- function(x, y, theta, alpha)
{
  d     <- length(theta)
  rho   <- theta[1L]
  mu    <- theta[2L]
  sigma <- sqrt(theta[3L])
  tau   <- theta[4L:d]
  num_y <- length(tau) + 1L
  
  ## compute unscaled weights
  wgt_unscaled <- weights_mainfun(
    x = x, y = y, alpha = alpha, rho = rho, 
    mu = mu, sigma = sigma, tau = tau, num_y = num_y)
  
  ## maximum weight
  maxwgt <- maxweight(theta = theta, alpha = alpha)
  
  ## rescaled weights
  # due to numerical approximation error, it may happen that
  # the computed maxweight is by a tiny bit smaller than the
  # largest empirical weight, so we truncate
  wgt_scaled <- wgt_unscaled / maxwgt
  wgt_scaled[wgt_scaled > 1] <- 1.0
  
  out <- list(weights_rescaled = wgt_scaled, weights_raw = wgt_unscaled, sup_weight = maxwgt)
  return(out)
}
