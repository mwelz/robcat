## the functions here are intended for efficiency comparisons. They assume that the polyserial model is correctly specified and that the supplied theta parameter is the true parameter value


## integrates p(z)^(1+alpha)*Q(z) over z
integral_over_Q <- function(theta, alpha)
{
  num_params <- length(theta)
  num_y      <- num_params - 2L
  alphaplus1 <- alpha + 1.0
  mat_out    <- matrix(NA_real_, num_params, num_params)
  
  for(i in seq_len(num_params))
  {
    ## j: j >= i
    for(j in which(1:num_params >= i))
    {
      
      ## function to solve integral over for dimension (i,j)
      fn <- Vectorize(function(x) {
        integrand_over_Qij(theta = theta, i = i, j = j, 
                           x = x, num_params = num_params, num_y = num_y, 
                           alphaplus1 = alphaplus1)
      })
      
      integrated <- stats::integrate(fn, lower = -Inf, upper = Inf)$value
      mat_out[i,j] <- mat_out[j,i] <- integrated
      
    } # FOR i
  } # FOR j
  return(mat_out)
}



## take integral in dimension (i,j) of p(z)^(1+alpha)*Q(z) over z
integrand_over_Qij <- function(theta, i, j, x, num_params, num_y, alphaplus1)
{
  rho <- theta[1L]
  mu <- theta[2L]
  sigma2 <- theta[3L]
  thres <- theta[4L:num_params]
  sigma <- sqrt(sigma2)
  p_marg_alphap1 <- px(x = x, mu = mu, sigma = sigma)^alphaplus1
  out <- 0.0
  
  for(y in seq_len(num_y))
  {
    p_cond_alphap1 <- py_x_scalars(y = y, x = x, rho = rho, mu = mu, 
                                   sigma = sigma, thres = thres, num_y = num_y)^alphaplus1
    
    Qij. <- Qij(theta = theta, i = i, j = j, x = x, y = y, num_params = num_params, num_y = num_y)
    summand <- p_cond_alphap1 * Qij.
    out <- out + summand
  } # FOR
  return(p_marg_alphap1 * out)
}



## integrates p(z)^(1+alpha)*u(z)*u(z)' over z
integral_over_uu <- function(theta, alpha)
{
  num_params <- length(theta)
  num_y      <- num_params - 2L
  alphaplus1 <- alpha + 1.0
  mat_out    <- matrix(NA_real_, num_params, num_params)
  
  for(i in seq_len(num_params))
  {
    ## j: j >= i
    for(j in which(1:num_params >= i))
    {
      
      ## function to solve integral over for dimension (i,j)
      fn <- Vectorize(function(x) {
        integrand_over_uij(theta = theta, i = i, j = j, 
                           x = x, num_params = num_params, num_y = num_y, 
                           alphaplus1 = alphaplus1)
      })
      
      integrated <- stats::integrate(fn, lower = -Inf, upper = Inf)$value
      mat_out[i,j] <- mat_out[j,i] <- integrated
      
    } # FOR i
  } # FOR j
  return(mat_out)
}



## take integral in dimension (i,j) of p(z)^(1+alpha)*u(z)*u(z)' over z
integrand_over_uij <- function(theta, i, j, x, num_params, num_y, alphaplus1)
{
  rho <- theta[1L]
  mu <- theta[2L]
  sigma2 <- theta[3L]
  thres <- theta[4L:num_params]
  sigma <- sqrt(sigma2)
  p_marg_alphap1 <- px(x = x, mu = mu, sigma = sigma)^alphaplus1
  out <- 0.0
  
  for(y in seq_len(num_y))
  {
    p_cond_alphap1 <- py_x_scalars(y = y, x = x, rho = rho, mu = mu, 
                                   sigma = sigma, thres = thres, num_y = num_y)^alphaplus1
    
    uij <- crossscore(theta = theta, i = i, j = j, x = x, y = y, num_params = num_params, num_y = num_y)
    summand <- p_cond_alphap1 * uij
    out <- out + summand
  } # FOR
  return(p_marg_alphap1 * out)
}


## the xi vector: take integral of p(z)^(1+alpha)*u(z) over z
xifun <- function(theta, alpha)
{
  num_params <- length(theta)
  num_y      <- num_params - 2L
  alphaplus1 <- alpha + 1.0
  vec_out    <- rep(NA_real_, num_params)
  
  for(i in seq_len(num_params))
  {
    ## function to solve integral over for dimension i
    fn <- Vectorize(function(x) {
      integrand_for_xi_i(theta = theta, i = i, 
                         x = x, num_params = num_params, num_y = num_y, 
                         alphaplus1 = alphaplus1)
    })
    
    integrated <- stats::integrate(fn, lower = -Inf, upper = Inf)$value
    vec_out[i] <- integrated
  }
  return(vec_out)
}



## take integral in dimension i of p(z)^(1+alpha)*u(z) over z
integrand_for_xi_i <- function(theta, i, x, num_params, num_y, alphaplus1)
{
  rho <- theta[1L]
  mu <- theta[2L]
  sigma2 <- theta[3L]
  thres <- theta[4L:num_params]
  sigma <- sqrt(sigma2)
  p_marg_alphap1 <- px(x = x, mu = mu, sigma = sigma)^alphaplus1
  out <- 0.0
  
  for(y in seq_len(num_y))
  {
    p_cond_alphap1 <- py_x_scalars(y = y, x = x, rho = rho, mu = mu, 
                                   sigma = sigma, thres = thres, num_y = num_y)^alphaplus1
    
    ui <- score_i(theta = theta, i = i, x = x, y = y, num_params = num_params, num_y = num_y)
    summand <- p_cond_alphap1 * ui
    out <- out + summand
  } # FOR
  return(p_marg_alphap1 * out)
}


population_Jmat_trumemodel <- function(theta, alpha)
{
  int_uu <- integral_over_uu(theta = theta, alpha = alpha)
  return(int_uu)
}

population_Kmat_trumemodel <- function(theta, alpha)
{
  int_uu <- integral_over_uu(theta = theta, alpha = 2*alpha)
  xi     <- xifun(theta = theta, alpha = alpha)
  return(int_uu - outer(xi, xi))
}


## calculate robust variance at correct model and true parameter
efficiency_robust <- function(theta, alpha)
{
  J <- population_Jmat_trumemodel(theta = theta, alpha = alpha)
  K <- population_Kmat_trumemodel(theta = theta, alpha = alpha)
  Jinv <- solve(J)
  return(Jinv %*% K %*% Jinv)
}


## calculate MLE variance at correct model and true parameter
efficiency_mle <- function(theta)
{
  fisher <- integral_over_Q(theta = theta, alpha = 0)
  return(solve(fisher))
}


#' Efficiency of minimum DPD estimators of polyserial model
#' 
#' Calculate population asymptotic variance-covariance matrix associated with a parameter vector \code{theta}, assuming that the polyserial model is correctly specified and that \code{theta} is the true model parameter. May take a few moments to compute because a relatively large number of integrals need to be numerically solved.
#' 
#' @param theta Parameter vector of polyserial model; assumed to be the true one. First element is polyserial correlatio coefficient, second is the population mean of X, third is population variance of X, and the remaining elements are the thresholds associated with the ordnal Y (must be in increasing order)
#' @param alpha Tuning constant governing robustness-efficiency tradeoff. Set to 0 for maximum likelihood.
#' 
#' @return 
#' A numeric matrix that is the population asymptotic variance-covariance matrix associated with a parameter vector \code{theta} and tuning constant \code{alpha}.
#' 
#' @examples
#' theta <- c(rho = 0, mu = 0, sigma2 = 1, tau1 = 0) # true parameter vector
#' polyserial_efficiency(theta, alpha = 0.5)
#' 
#' @export
polyserial_efficiency <- function(theta, alpha)
{
  stopifnot(alpha >= 0)
  if(alpha > 0)
  {
    asv <- efficiency_robust(theta = theta, alpha = alpha)
  } else
  {
    asv <- efficiency_mle(theta = theta)
  }
  nam <- polyserial_thetanames(length(theta) - 2L)
  rownames(asv) <- colnames(asv) <- nam
  return(asv)
}
