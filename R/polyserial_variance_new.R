get_matrix_A <- function(rho, mu, sigma2, thres, num_y, alphaplus1)
{
  ## initialize
  num_params <- num_y + 2L
  mat_out <- matrix(NA_real_, num_params, num_params)


  for(i in seq_len(num_params))
  {
    ## get functions
    px_d1_i   <- get_px_d1_i(i)
    py_x_d1_i <- get_py_x_d1_i(i)
    
    ## j: j >= i
    for(j in which(1:num_params >= i))
    {
      ## get functions
      px_d1_j    <- get_px_d1_i(j)
      py_x_d1_j  <- get_py_x_d1_i(j)
      px_d2_ij   <- get_px_d2_ij(i,j)
      py_x_d2_ij <- get_py_x_d2_ij(i,j)
      
      ## function to solve integral over for dimension (i,j)
      fn <- Vectorize(function(x) {
        A_ij_integrand_x(x = x, rho = rho, mu = mu, 
                         sigma2 = sigma2, thres = thres, num_y = num_y, 
                         alphaplus1 = alphaplus1, 
                         px_d1_i = px_d1_i, px_d1_j = px_d1_j, 
                         py_x_d1_i = py_x_d1_i, py_x_d1_j = py_x_d1_j, 
                         px_d2_ij = px_d2_ij, py_x_d2_ij = py_x_d2_ij)
      }) 
      
      integrated <- stats::integrate(fn, lower = -Inf, upper = Inf)$value
      mat_out[i,j] <- mat_out[j,i] <- integrated
    } # FOR i
  } # FOR j
  
  return(mat_out)
}


A_ij_integrand_x <- function(x, rho, mu, sigma2, thres, num_y, 
                             alphaplus1,
                             px_d1_i, 
                             px_d1_j,
                             py_x_d1_i,
                             py_x_d1_j,
                             px_d2_ij,
                             py_x_d2_ij)
{
  ## marginal information
  sigma <- sqrt(sigma2)
  marg_p <- px(x = x, mu = mu, sigma = sigma)
  marg_p_ap1 <- marg_p^alphaplus1
  marg_p_inv <- ifelse(marg_p > 0.0, 1.0 / marg_p, 0.0) # avoid divide by 0
  marg_grad_i <- px_d1_i(x = x, mu = mu, sigma2 = sigma2)
  marg_grad_j <- px_d1_j(x = x, mu = mu, sigma2 = sigma2)
  marg_hess_ij <- px_d2_ij(x = x, mu = mu, sigma2 = sigma2)
  marg_score_i <- marg_p_inv * marg_grad_i
  marg_score_j <- marg_p_inv * marg_grad_j
  marg_score_prod <- marg_score_i * marg_score_j
  marg_Qij <- marg_score_prod - marg_p_inv * marg_hess_ij
  
  ## init
  out <- 0.0
  
  for(y in seq_len(num_y))
  {
    cond_p <- py_x_scalars(y = y, x = x, rho = rho, mu = mu, 
                           sigma = sigma, thres = thres, num_y = num_y)
    cond_p_ap1 <- cond_p^alphaplus1
    cond_p_inv <- ifelse(cond_p > 0.0, 1.0 / cond_p, 0.0) # avoid divide by 0
    cond_grad_i <- py_x_d1_i(y = y, x = x, rho = rho, mu = mu, 
                             sigma2 = sigma2, thres = thres, num_y = num_y)
    cond_grad_j <- py_x_d1_j(y = y, x = x, rho = rho, mu = mu, 
                             sigma2 = sigma2, thres = thres, num_y = num_y)
    cond_hess_ij <- py_x_d2_ij(y = y, x = x, rho = rho, mu = mu, 
                               sigma2 = sigma2, thres = thres, num_y = num_y)
    cond_score_i <- cond_p_inv * cond_grad_i
    cond_score_j <- cond_p_inv * cond_grad_j
    cond_score_prod <- cond_score_i * cond_score_j
    cond_Qij <- cond_score_prod - cond_p_inv * cond_hess_ij
    
    ## assemble final summand
    Q_ij <- cond_Qij + marg_Qij
    uiuj <- (marg_score_i + cond_score_i) * (marg_score_j + cond_score_j)
    increment <- cond_p_ap1 * (alphaplus1 * uiuj - Q_ij)
    
    ## update
    out <- out + increment
  }
  return(marg_p_ap1 * out)
}





# mu <- 0
# sigma2 <- 1.2
# rho <- 0.5
# thres <- c(-1, -0)
# y_set <- c(1L, 2L, 3L)
# num_y <- length(y_set)
# theta <- c(rho, mu, sigma2, thres)
# alphaplus1 = 1.1
# i = 1L
# j = 2L
# d <- length(theta)
# set.seed(1)
# x <- rnorm(n = 1, mean = mu, sd = sqrt(sigma2))
# f = get_px_d1_i(4L)
# f(x,mu,sigma2)
