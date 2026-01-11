# exception: here x and y are vectors of observations, no longer scalars
get_matrix_B <- function(x, y, rho, mu, sigma2, thres, num_y, alpha)
{
  ## calculate densities at empirical observations
  sigma <- sqrt(sigma2)
  p_marg_alpha <- px(x = x, mu = mu, sigma = sigma)^alpha
  p_cond_alpha <- py_x_vectorized(
    y = y, x = x, rho = rho, mu = mu, sigma = sigma, thres = thres, num_y = num_y)^alpha
  
  ## initialize
  num_params <- num_y + 2L
  out_mat <- matrix(0.0, nrow = num_params, ncol = num_params)
  N <- length(x)
  
  for(i in seq_len(N))
  {
    Q <- polyserial_Qmat(x = x[i], y = y[i], rho = rho, mu = mu, 
                         sigma2 = sigma2, thres = thres, num_y = num_y)
    u <- score_z(x = x[i], y = y[i], rho = rho, mu = mu, sigma = sigma, 
                 thres = thres, num_y = num_y)
    
    increment <- p_marg_alpha[i] * p_cond_alpha[i] * (Q - alpha * u %*% t(u))
    out_mat <- out_mat + increment
  }
  return(out_mat / N)
}


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
      
      ## solve integral
      # in rare cases, stats::integrate throws a cryptic error. Using 
      # a different solver seems to fix the issue
      # FIXME: it seems that fn() sometimes returns an NaN, which triggers
      # the issue. I suspect this is caused by a 0 * Inf operation somewhere in 
      # the finite grid we evaluate over; investigate!
      integrated <- try(
        stats::integrate(fn, -Inf, Inf)$value, 
        silent = TRUE
      )
      if (inherits(integrated, "try-error"))
      {
        invisible(utils::capture.output(
          integrated <- pracma::integral(fn, -Inf, Inf)))
        
      }

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


# exception: here x and y are vectors of observations, no longer scalars
get_Jhat <- function(x, y, rho, mu, sigma2, thres, num_y, alpha)
{
  A <- get_matrix_A(
    rho = rho, mu = mu, sigma2 = sigma2, thres = thres, 
    num_y = num_y, alphaplus1 = alpha + 1.0)
  B <- get_matrix_B(
    x = x, y = y, rho = rho, mu = mu, sigma2 = sigma2, thres = thres,
    num_y = num_y, alpha = alpha)
  return(A + B)
}


# exception: here x and y are vectors of observations, no longer scalars
get_K_xihat <- function(x, y, rho, mu, sigma2, thres, num_y, alpha)
{
  
  ## densities
  sigma  <- sqrt(sigma2)
  p_marg <- px(x = x, mu = mu, sigma = sigma)
  p_cond <- py_x_vectorized(y = y, x = x, rho = rho, mu = mu, sigma = sigma,
                            thres = thres, num_y = num_y)
  prob <- p_marg * p_cond
  prob_alpha <- prob^alpha
  prob_alpha2 <- prob^(2*alpha)
  N <- length(x)
  
  ## score vectors
  score_vecs <- sapply(seq_len(N), function(i){
    score_z(x = x[i], y = y[i], rho = rho, mu = mu, sigma = sigma,
            thres = thres, num_y = num_y)
  })
  
  ## xi hat
  summands_xi <- sapply(seq_len(N), function(i){
    prob_alpha[i] * score_vecs[,i]
  })
  xihat <- rowMeans(summands_xi)
  
  ## K hat: initialize
  d <- length(xihat)
  summands_K <- matrix(0.0, d, d)
  
  for(i in seq_len(N))
  {
    u <- score_vecs[,i]
    summand <- prob_alpha2[i] * u %*% t(u)
    summands_K <- summands_K + summand
  } # FOR
  
  Khat <- summands_K / N - xihat %*% t(xihat)
  
  return(list(Khat = Khat, xihat = xihat))
} # FUN


# exception: here x and y are vectors of observations, no longer scalars
get_Sigmahat_polyserial <- function(x, y, rho, mu, sigma2, thres, num_y, alpha)
{
  J <- get_Jhat(x = x, y = y, rho = rho, mu = mu, sigma2 = sigma2, 
                thres = thres, num_y = num_y, alpha = alpha)
  K <- get_K_xihat(x = x, y = y, rho = rho, mu = mu, sigma2 = sigma2,
                   thres = thres, num_y = num_y, alpha = alpha)$Khat
  
  Jinv <- solve(J)
  return(Jinv %*% K %*% Jinv)
}


# exception: here x and y are vectors of observations, no longer scalars
variance_polyserial <- function(theta, x, y, alpha)
{
  d <- length(theta)
  num_y <- d - 2L
  N <- length(x)
  
  rho <- theta[1L]
  mu <- theta[2L]
  sigma2 <- theta[3L]
  tau <- theta[4L:d]
  
  # ASV of sqrt(N) * thetahat
  asv <- get_Sigmahat_polyserial(
    x = x, y = y, rho = rho, mu = mu, sigma2 = sigma2,
    thres = tau, num_y = num_y, alpha = alpha)
  
  # ASV of thetahat
  vcov <- asv / N
  colnames(vcov) <- rownames(vcov) <- polyserial_thetanames(num_y)
  
  se <- sqrt(diag(vcov))
  
  return(list(se = se, vcov = vcov))
}



variance_polyserial_mle <- function(theta, x, y)
{
  d      <- length(theta)
  N      <- length(x)
  num_y  <- d - 2L
  rho    <- theta[1L]
  mu     <- theta[2L]
  sigma  <- sqrt(theta[3L])
  thres  <- theta[4L:d]
  
  ## score vectors
  score_vecs <- sapply(seq_len(N), function(i){
    score_z(x = x[i], y = y[i], rho = rho, mu = mu, sigma = sigma,
            thres = thres, num_y = num_y)
  })
  
  hess <- matrix(0.0, d, d)
  for(i in seq_len(N))
  {
    u <- score_vecs[,i]
    summand <- u %*% t(u)
    hess <- hess + summand
  } # FOR
  
  # ASV of thetahat (not sqrt(N) * thetahat b/c of the summation in objective)
  vcov <- solve(hess)
  colnames(vcov) <- rownames(vcov) <- names(theta)
  se <- sqrt(diag(vcov))
  return(list(se = se, vcov = vcov))
}


# mu <- 0
# sigma2 <- 1.2
# rho <- 0.5
# thres <- c(-1, 0, 1,1.5)
# y_set <- c(1L, 2L, 3L, 4L, 5L)
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