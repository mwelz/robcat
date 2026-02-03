## assemble hessian matrix of p_{Y|X}
py_x_hessian <- function(y, x, rho, mu, sigma2, thres, num_y)
{
  ## initialize
  d <- num_y + 2L
  hess <- matrix(0.0, nrow = d, ncol = d)
  args <- list(y = y, x = x, rho = rho, mu = mu, 
               sigma2 = sigma2, thres = thres, num_y = num_y)
  
  ## submatrix w/o threshold entries
  hess[1L,1L] <- do.call(py_x_d2_rho2,    args)
  hess[2L,2L] <- do.call(py_x_d2_mu2,     args)
  hess[3L,3L] <- do.call(py_x_d2_sigma22, args)
  hess[1L,2L] <- hess[2L,1L] <- do.call(py_x_d2_rhomu,     args)
  hess[1L,3L] <- hess[3L,1L] <- do.call(py_x_d2_rhosigma2, args)
  hess[2L,3L] <- hess[3L,2L] <- do.call(py_x_d2_musigma2,  args)
  
  
  ## now fill up remaining derivatives involving the thresholds
  for(k in seq_along(thres))
  {
    idx <- k + 3L # index in theta vector
    argsk <- list(y = y, x = x, k = k, rho = rho, mu = mu, 
                  sigma2 = sigma2, thres = thres, num_y = num_y)
    
    hess[idx, idx] <- do.call(py_x_d2_tauk2, argsk)
    hess[1L,idx]   <- hess[idx,1L] <- do.call(py_x_d2_rhotauk, argsk)
    hess[2L,idx]   <- hess[idx,2L] <- do.call(py_x_d2_mutauk, argsk)
    hess[3L,idx]   <- hess[idx,3L] <- do.call(py_x_d2_sigma2tauk, argsk)
    # NB: the cross-derivatives of distinct tau parameters are structural 0s,
    # so no need to explicitly code them here
  } # FOR
  
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
# y <- 2L
# 
# f <- function(theta) py_x_scalars(y = y, x = x, 
#                                   rho = theta[1], mu = theta[2],
#                                   sigma = sqrt(theta[3]), thres = theta[4:d],
#                                   num_y = num_y)
# apprx <- numDeriv::hessian(f, theta)
# exact <- py_x_hessian(y = y, x = x, rho = rho, mu = mu, sigma2 = sigma2,
#                       thres = thres, num_y = num_y)
# all.equal(exact, apprx, tolerance = tol) # TRUE
