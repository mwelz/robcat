#### first order derivatives: all scalar-input, vectorize later in componentwise integral! ### ----

## gradient of p(x; theta)
px_grad <- function(x, mu, sigma, num_y)
{
  pxtheta <- px(x = x, mu = mu, sigma = sigma)
  sigma2  <- sigma^2
  xmmu    <- x - mu
  
  ## derivative wrt mu
  prime_mu <- pxtheta * xmmu / sigma2
  
  ## derivative wrt sigma2
  prime_sigma2 <- pxtheta * (xmmu^2 / sigma2 - 1.0) / (2.0 * sigma2)
  
  ## the remaining derivatives are all zero
  out <- c(0.0, prime_mu, prime_sigma2, rep(0.0, num_y - 1L))
  return(out)
}


## gradient of p(y | x; theta)
py_x_grad <- function(y, x, rho, mu, sigma, thres, num_y)
{
  ## FIXME: this is inefficient; we don't need the whole PMF or all associated taustars
  taustar_all  <- get_taustar(x = x, rho = rho, mu = mu, sigma = sigma, thres = thres) # len = r-1
  
  ## standardized x
  x_st <- (x - mu) / sigma
  
  ## denominators
  onem_rho2      <- 1.0 - rho^2
  sqrt_onem_rho2 <- sqrt(onem_rho2)
  
  ## initialize derivatives wrt tau
  grad_tau <- rep(0.0, num_y - 1L)
  
  
  ## evaluate relevant taustars at N(0,1) density
  if(y == 1L)
  {
    phi_y   <- stats::dnorm(taustar_all[y])
    phi_ym1 <- 0.0
    
    phi_y_tau   <- phi_y * (thres[y] * rho - x_st)
    phi_ym1_tau <- 0.0
    
    ## derivatives wrt tau
    grad_tau[y] <- phi_y / sqrt_onem_rho2
    
  } else if(y == num_y)
  {
    phi_y   <- 0.0
    phi_ym1 <- stats::dnorm(taustar_all[y - 1L])
    
    phi_y_tau <- 0.0
    phi_ym1_tau <- phi_ym1 * (thres[y-1L] * rho - x_st)
    
    ## derivatives wrt tau
    grad_tau[y-1L] <- (-1) * phi_ym1 / sqrt_onem_rho2
    
  } else
  {
    phi_y   <- stats::dnorm(taustar_all[y])
    phi_ym1 <- stats::dnorm(taustar_all[y - 1L])
    
    phi_y_tau   <- phi_y * (thres[y] * rho - x_st)
    phi_ym1_tau <- phi_ym1 * (thres[y-1L] * rho - x_st)
    
    ## derivatives wrt tau
    grad_tau[y]    <- phi_y / sqrt_onem_rho2
    grad_tau[y-1L] <- (-1) * phi_ym1 / sqrt_onem_rho2
  } # IF
  
  
  ## derivative wrt rho
  grad_rho <- onem_rho2^(-1.5) * (phi_y_tau - phi_ym1_tau)
  
  ## derivative wrt mu
  phi_diff <- phi_y - phi_ym1
  grad_mu <- rho * phi_diff / (sigma * sqrt_onem_rho2)
  
  ## derivative wrt sigma2
  grad_sigma2 <- rho * x_st * phi_diff / (2.0 * sigma^2 * sqrt_onem_rho2)
  
  ## return gradient wrt theta
  return(c(grad_rho, grad_mu, grad_sigma2, grad_tau))
}


## first-order derivatives of p_{X}
px_d1_mu <- function(x, mu, sigma2)
{
  if(is.infinite(x))
  {
    out <- 0.0
  } else
  {
    pxtheta <- px(x = x, mu = mu, sigma = sqrt(sigma2))
    xmmu    <- x - mu
    out <- pxtheta * xmmu / sigma2
  }
  return(out)
}

px_d1_sigma2 <- function(x, mu, sigma2)
{
  if(is.infinite(x))
  {
    out <- 0.0
  } else
  {
    pxtheta <- px(x = x, mu = mu, sigma = sqrt(sigma2))
    xmmu    <- x - mu
    out <- pxtheta * (xmmu^2 / sigma2 - 1.0) / (2.0 * sigma2)
  }
  return(out)
}



## first-order derivatives of p_{Y|X}
py_x_d1_rho <- function(y, x, rho, mu, sigma2, thres, num_y)
{
  scaling <- (1.0 - rho^2)^(-1.5)
  z <- (x - mu) / sqrt(sigma2)
  ym1 <- y - 1L
  
  if(y == 1L)
  {
    p1_1 <- phi(get_taustar_y(x = x, rho = rho, mu = mu, sigma2 = sigma2, 
                              tau_y = thres[y]))
    p1_2 <- rho * thres[y] - z
    p2_1 <- p2_2 <- 0.0
  } else if(y == num_y)
  {
    p1_1 <- p1_2 <- 0.0
    p2_1 <- phi(get_taustar_y(x = x, rho = rho, mu = mu, sigma2 = sigma2, 
                              tau_y = thres[ym1]))
    p2_2 <- rho * thres[ym1] - z
  } else
  {
    p1_1 <- phi(get_taustar_y(x = x, rho = rho, mu = mu, sigma2 = sigma2, 
                              tau_y = thres[y]))
    p1_2 <- rho * thres[y] - z
    
    p2_1 <- phi(get_taustar_y(x = x, rho = rho, mu = mu, sigma2 = sigma2, 
                              tau_y = thres[ym1]))
    p2_2 <- rho * thres[ym1] - z
  }
  out <- scaling * (p1_1 * p1_2 - p2_1 * p2_2)
  if(is.nan(out)) out <- 0.0 ## account for 0 * Inf
  return(out)
}


py_x_d1_mu <- function(y, x, rho, mu, sigma2, thres, num_y)
{
  scaling <- rho / sqrt(sigma2 * (1.0 - rho^2))
  
  if(y == 1L)
  {
    p1 <- phi(get_taustar_y(x = x, rho = rho, mu = mu, sigma2 = sigma2, 
                            tau_y = thres[y]))
    p2 <- 0.0
  } else if(y == num_y)
  {
    p1 <- 0.0
    p2 <- phi(get_taustar_y(x = x, rho = rho, mu = mu, sigma2 = sigma2, 
                            tau_y = thres[y-1L]))
  } else
  {
    p1 <- phi(get_taustar_y(x = x, rho = rho, mu = mu, sigma2 = sigma2, 
                            tau_y = thres[y]))
    p2 <- phi(get_taustar_y(x = x, rho = rho, mu = mu, sigma2 = sigma2, 
                            tau_y = thres[y-1L]))
  }
  return(scaling * (p1 - p2))
}


py_x_d1_sigma2 <- function(y, x, rho, mu, sigma2, thres, num_y)
{
  scaling <- rho * (x-mu) / ( 2.0 * sigma2^1.5 * sqrt(1.0 - rho^2))
  
  if(y == 1L)
  {
    p1 <- phi(get_taustar_y(x = x, rho = rho, mu = mu, sigma2 = sigma2, 
                            tau_y = thres[y]))
    p2 <- 0.0
  } else if(y == num_y)
  {
    p1 <- 0.0
    p2 <- phi(get_taustar_y(x = x, rho = rho, mu = mu, sigma2 = sigma2, 
                            tau_y = thres[y-1L]))
  } else
  {
    p1 <- phi(get_taustar_y(x = x, rho = rho, mu = mu, sigma2 = sigma2, 
                            tau_y = thres[y]))
    p2 <- phi(get_taustar_y(x = x, rho = rho, mu = mu, sigma2 = sigma2, 
                            tau_y = thres[y-1L]))
  }
  
  out <- scaling * (p1 - p2)
  if(is.nan(out)) out <- 0.0 ## account for 0 * Inf
  return(out)
}

# defined for k \in [num_y-1]
py_x_d1_tauk <- function(y, x, k, rho, mu, sigma2, thres, num_y)
{
  stopifnot(k < num_y)
  ym1     <- y - 1L
  scaling <- 1.0 / sqrt(1.0 - rho^2)
  
  if(k == y)
  {
    tmp <- phi(get_taustar_y(x = x, rho = rho, mu = mu, sigma2 = sigma2, 
                             tau_y = thres[y]))
    out <- scaling * tmp
  } else if(k == ym1)
  {
    tmp <- phi(get_taustar_y(x = x, rho = rho, mu = mu, sigma2 = sigma2, 
                             tau_y = thres[ym1]))
    out <- (-1.0) * scaling * tmp
  } else
  {
    out <- 0.0
  }
  return(out)
}


## sanity check
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
# y <- 3L
# tol <- 0.0001
# target <- py_x_grad(y = y, x = x, rho = rho, mu = mu, 
#                     sigma = sqrt(sigma2), thres = thres, num_y = num_y)
# args <- list(y=y, x=x, rho=rho, mu=mu, sigma2=sigma2, thres=thres, num_y=num_y)
# grad <- rep(NA, d)
# grad[1] <- do.call(py_x_d1_rho, args)
# grad[2] <- do.call(py_x_d1_mu, args)
# grad[3] <- do.call(py_x_d1_sigma2, args)
# grad[4:d] <- sapply(seq_along(thres), function(k){
#   py_x_d1_tauk(y,x,k,rho,mu,sigma2,thres,num_y)
# })
# all.equal(target, grad, tolerance = tol) # TRUE
