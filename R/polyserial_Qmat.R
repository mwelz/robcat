## hessian of log density, times (-1)
polyserial_Qmat <- function(x, y, rho, mu, sigma2, thres, num_y)
{
  sigma <- sqrt(sigma2)
  
  ## marginal density
  marg_dens <- px(x = x, mu = mu, sigma = sigma)
  marg_grad <- px_grad(x = x, mu = mu, sigma = sigma, num_y = num_y)
  marg_hess <- px_hessian(x = x, mu = mu, sigma2 = sigma2, num_y = num_y)
  
  ## conditional density
  cond_dens <- py_x_scalars(y = y, x = x, rho = rho, mu = mu, 
                            sigma = sigma, thres = thres, num_y = num_y)
  cond_grad <- py_x_grad(y = y, x = x, rho = rho, mu = mu, 
                         sigma = sigma, thres = thres, num_y = num_y)
  cond_hess <- py_x_hessian(y = y, x = x, rho = rho, mu = mu, 
                            sigma2 = sigma2, thres = thres, num_y = num_y)
  
  ## if a density is exactly zero, cast its reciprocal as 0
  # this has no effect because of our convention that Inf * 0 = 0 :
  # in the case of a zero density, density's gradient and hessian
  # will be zero as well
  inv_px   <- ifelse(marg_dens > 0.0, 1.0 / marg_dens, 0.0)
  inv_py_x <- ifelse(cond_dens > 0.0, 1.0 / cond_dens, 0.0)
  
  ## assemble final term (matrix Q)
  diff_cond <- 
    inv_py_x * cond_hess - inv_py_x^2 * outer(cond_grad, cond_grad)
  diff_marg <- 
    inv_px * marg_hess - inv_px^2 * outer(marg_grad, marg_grad)
  hess_logdens <- diff_cond + diff_marg
  return((-1) * hess_logdens)
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
# apprx <- Qz(x = x, y = y, rho = rho, mu = mu, sigma = sqrt(sigma2), 
#             thres = thres, num_y = num_y)
# exact <- polyserial_Qmat(x = x, y = y, rho = rho, mu = mu, sigma2 = sigma2,
#                          thres = thres, num_y = num_y)
# all.equal(exact, apprx, tolerance = tol) # TRUE
