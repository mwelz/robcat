## functions for the ASV

#### first order derivatives: all scalar-input, vectorize later in componentwise integral! ### ----


## u_theta(z): derivative of log p_theta(z) wrt theta
score_z <- function(x, y, rho, mu, sigma, thres, num_y)
{
  ## We run into a divide-by-zero issue for extremely large values of x because
  ## densities are (computationally) zero at such values, regardless of parameter
  ## values. To avoid this issue, return 0 if standardized x exceeds value 15,
  ## which has virtually zero probability mass at the assumed normal model
  if(abs((x - mu) / sigma) > 10) # TODO: maybe more general is better: as soon as p = 0
  {
    out <- rep(0.0, num_y + 2L)
  } else
  {
    ## densities
    p_cond <- py_x_scalars(y = y, x = x, rho = rho, mu = mu, sigma = sigma, thres = thres, num_y = num_y)
    p_marg <- px(x = x, mu = mu, sigma = sigma)
    
    ## gradients
    p_cond_grad <- py_x_grad(y = y, x = x, rho = rho, mu = mu, sigma = sigma, thres = thres, num_y = num_y)
    p_marg_grad <- px_grad(x = x, mu = mu, sigma = sigma, num_y = num_y)
    
    out <- p_cond_grad / p_cond + p_marg_grad / p_marg
  } # IF
  
  return(out)
} # FUN


## derivative of score wrt theta
# approximate numerically for now
# NB: behavior of score_z for large x could cause issues again
Qz <- function(x, y, rho, mu, sigma, thres, num_y)
{
  theta <- c(rho, mu, sigma^2, thres)
  fn <- function(theta)
  {
    score_z(x = x, y = y, 
            rho = theta[1L],
            mu = theta[2L], 
            sigma = sqrt(theta[3L]), 
            thres = theta[4L:(num_y+2L)],
            num_y = num_y)
  }
  
  jac <- numDeriv::jacobian(fn, theta)
  return(jac * (-1.0))
}


# ## for the integral in J-hat
# J_integrand <- function(x, rho, mu, sigma, thres, num_y, alphaplus1)
# {
#   p_marg_alphap1 <- px(x = x, mu = mu, sigma = sigma)^alphaplus1
#   num_params <- num_y + 2L
#   matsum <- matrix(0.0, nrow = num_params, ncol = num_params)
#   
#   for(y in seq_len(num_y))
#   {
#     p_cond_alphap1 <- py_x_scalars(y = y, x = x, rho = rho, mu = mu, 
#                                    sigma = sigma, thres = thres, num_y = num_y)
#     
#     u <- score_z(x = x, y = y, rho = rho, mu = mu, sigma = sigma, thres = thres, num_y = num_y)
#     Q <- Qz(x = x, y = y, rho = rho, mu = mu, sigma = sigma, thres = thres, num_y = num_y)
#     
#     summand <- p_cond_alphap1 * (alphaplus1 * u %*% t(u) - Q)
#     matsum <- matsum + summand
#   } # FOR
#   
#   return(p_marg_alphap1 * matsum)
# }
# 
# 
# ##  solve the integral in Jhat component-wisely
# # TODO: this is very slow, can be improved a lot!
# Jhat_integral <- function(rho, mu, sigma, thres, num_y, alphaplus1)
# {
#   ## initialize
#   num_params <- num_y + 2L
#   mat_out <- matrix(NA_real_, num_params, num_params)
#   
#   for(i in seq_len(num_params))
#   {
#     ## j: j >= i
#     for(j in which(1:num_params >= i))
#     {
#       ## function that picks out (i,j)-th coordinate of the integrand matrix
#       fn <- function(x) {
#         J_integrand(x = x, rho = rho, mu = mu, sigma = sigma, thres = thres, num_y = num_y, alphaplus1 = alphaplus1)[i,j]
#       }
#       
#       ## solve integral at coordinate (i,j)
#       # FIXME: handling of large x values is still a bit shakey, see score fun!
#       # Also, this is shittily programmed: very wasteful since in each iteration we compute the entire matrix just to drop everything except component (i,j). can be sped up a lot!
#       integrated <- stats::integrate(Vectorize(fn), lower = -Inf, upper = Inf)$value
#       mat_out[i,j] <- mat_out[j,i] <- integrated
#       
#     } # FOR i
#   } # FOR j
#   
#   return(mat_out)
# } # FUN


# order: rho, mu, sigma2, thres
logl_scalar <- function(theta, x, y, num_params, num_y)
{
  rho <- theta[1L]
  mu <- theta[2L]
  sigma2 <- theta[3L]
  thres <- theta[4L:num_params]
  sigma <- sqrt(sigma2)
  
  marg <- px(x = x, mu = mu, sigma = sigma)
  cond <- py_x_scalars(y = y, x = x, rho = rho, mu = mu, sigma = sigma, thres = thres, num_y = num_y)
  prob <- marg * cond
  out <- ifelse(prob > 0.0, log(prob), 0.0) # pathological cases with 0 prob get 0 likelihood (e.g. Inf; cancels anyway)
  return(out)
}


## first derivative of logl wrt i-th parameter
score_i <- function(theta, i, x, y, num_params, num_y)
{
  theta_before_i <- theta[seq_len(i-1L)]
  theta_after_i <- theta[i + seq_len(num_params - i)]
  
  fn <- function(thetai) logl_scalar(theta = c(theta_before_i, thetai, theta_after_i),
                                     x = x, y = y, num_params = num_params, num_y = num_y)
  
  return(numDeriv::grad(fn, x = theta[i]))
}

## multiply i-th score by j-th score
crossscore <- function(theta, i, j, x, y, num_params, num_y)
{
  ui <- score_i(theta = theta, i = i, x = x, y = y, num_params = num_params, num_y = num_y)
  
  if(i == j)
  {
    uj <- ui
  } else
  {
    uj <- score_i(theta = theta, i = j, x = x, y = y, num_params = num_params, num_y = num_y)
  }
  
  return(ui * uj)
}


## negative second derivative of logl wrt i-th parameter
Qii <- function(theta, i, x, y, num_params, num_y)
{
  theta_before_i <- theta[seq_len(i-1L)]
  theta_after_i <- theta[i + seq_len(num_params - i)]
  
  fn <- function(thetai) logl_scalar(theta = c(theta_before_i, thetai, theta_after_i),
                                     x = x, y = y, num_params = num_params, num_y = num_y)
  
  return(-as.numeric(numDeriv::hessian(fn, x = theta[i])))
}


## negative cross-derivatives of logl
# i <= j!!!
Qij <- function(theta, i, j, x, y, num_params, num_y)
{
  if(i == j)
  {
    out <- Qii(theta = theta, i = i, x = x, y = y, num_params = num_params, num_y = num_y)
  } else if (i < j)
  {
    theta_before_i <- theta[seq_len(i-1L)]
    theta_btw_i_and_j <- theta[i + seq_len(j - i - 1L)]
    theta_after_j <- theta[j + seq_len(num_params - j)]
    
    fn <- function(thetaij) logl_scalar(theta = c(theta_before_i, thetaij[1], theta_btw_i_and_j, thetaij[2], theta_after_j),
                                        x = x, y = y, num_params = num_params, num_y = num_y)
    out <- -numDeriv::hessian(fn, x = c(theta[i], theta[j]))[1,2]
  } else stop(" i <= j violated")
  
  return(out)
}



#score_z(x = x, y = y, rho = rho, mu = mu, sigma = sigma, thres = thres, num_y = num_y)
#sapply(1:num_params, function(i) score_i(theta = theta, i = i, x = x, y = y, num_y = num_y, num_params = num_params)) # yay



## for the integral in J-hat: x is a scalar
J_integrand_ij <- function(theta, i, j, x, num_params, num_y, alphaplus1)
{
  rho <- theta[1L]
  mu <- theta[2L]
  sigma2 <- theta[3L]
  thres <- theta[4L:num_params]
  sigma <- sqrt(sigma2)
  p_marg_alphap1 <- px(x = x, mu = mu, sigma = sigma)^alphaplus1
  
  ## init
  out <- 0.0
  
  for(y in seq_len(num_y))
  {
    p_cond_alphap1 <- py_x_scalars(y = y, x = x, rho = rho, mu = mu, 
                                   sigma = sigma, thres = thres, num_y = num_y)^alphaplus1
    
    uiuj. <- crossscore(theta = theta, i = i, j = j, x = x, y = y, num_params = num_params, num_y = num_y)
    Qij. <- Qij(theta = theta, i = i, j = j, x = x, y = y, num_params = num_params, num_y = num_y)
    
    summand <- p_cond_alphap1 * (alphaplus1 * uiuj. - Qij.)
    out <- out + summand
  } # FOR
  
  return(p_marg_alphap1 * out)
}


## the 1st part in Jhat (the integral one): matrix A
Jhat_integral <- function(rho, mu, sigma, thres, num_y, alphaplus1)
{
  ## initialize
  num_params <- num_y + 2L
  mat_out <- matrix(NA_real_, num_params, num_params)
  theta <- c(rho, mu, sigma^2, thres)
  
  
  for(i in seq_len(num_params))
  {
    ## j: j >= i
    for(j in which(1:num_params >= i))
    {
      
      ## function to solve integral over for dimension (i,j)
      fn <- Vectorize(function(x) {
        J_integrand_ij(theta = theta, i = i, j = j, x = x, num_params = num_params, num_y = num_y, alphaplus1 = alphaplus1)
      })
      
      integrated <- stats::integrate(fn, lower = -Inf, upper = Inf)$value
      mat_out[i,j] <- mat_out[j,i] <- integrated
      
    } # FOR i
  } # FOR j
  
  return(mat_out)
} # FUN


## the 2nd part in Jhat (the empirical one): matrix Bhat
# exception: here x and y are vectors of observations, no longer scalars
Jhat_empirical <- function(x, y, rho, mu, sigma, thres, num_y, alpha)
{
  
  ## calculate densities at empirical observations
  p_marg_alpha <- px(x = x, mu = mu, sigma = sigma)^alpha
  p_cond_alpha <- py_x_vectorized(
    y = y, x = x, rho = rho, mu = mu, sigma = sigma, thres = thres, num_y = num_y)^alpha
  
  ## initialize
  num_params <- num_y + 2L
  out_mat <- matrix(0.0, nrow = num_params, ncol = num_params)
  N <- length(x)
  
  for(i in seq_len(N))
  {
    Q <- Qz(x = x[i], y = y[i], rho = rho, mu = mu, sigma = sigma, thres = thres, num_y = num_y)
    u <- score_z(x = x[i], y = y[i], rho = rho, mu = mu, sigma = sigma, thres = thres, num_y = num_y)
    
    increment <- p_marg_alpha[i] * p_cond_alpha[i] * (Q - alpha * u %*% t(u))
    out_mat <- out_mat + increment
  }
  return(out_mat / N)
} # FUN


# exception: here x and y are vectors of observations, no longer scalars
get_Jhat <- function(x, y, rho, mu, sigma, thres, num_y, alpha)
{
  A <- Jhat_integral(
    rho = rho, mu = mu, sigma = sigma, thres = thres, num_y = num_y, alphaplus1 = alpha + 1.0)
  B <- Jhat_empirical(
    x = x, y = y, rho = rho, mu = mu, sigma = sigma, thres = thres, num_y = num_y, alpha = alpha)
  return(A + B)
}


# exception: here x and y are vectors of observations, no longer scalars
get_K_xihat <- function(x, y, rho, mu, sigma, thres, num_y, alpha)
{
  
  ## densities
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
get_Sigmahat_polyserial <- function(x, y, rho, mu, sigma, thres, num_y, alpha)
{
  J <- get_Jhat(x = x, y = y, rho = rho, mu = mu, sigma = sigma, 
                thres = thres, num_y = num_y, alpha = alpha)
  K <- get_K_xihat(x = x, y = y, rho = rho, mu = mu, sigma = sigma,
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
  asv <- get_Sigmahat_polyserial(x = x, y = y, rho = rho, mu = mu, sigma = sqrt(sigma2),
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

# keep y fixed and have x vectorized
#py_x_vectorized_x <- Vectorize(py_x_scalars, vectorize.args = "x")

# x = 1.548
# y <- 3L
# rho <- 0.5
# sigma <- 1
# sigma2 <- sigma^2
# mu <- 0
# thres <- c(-1.5, -0.5, 0.5, 1.5)
# theta <- c(rho, mu, sigma2, thres)
# num_y <- 5L
# num_params <- num_y + 2L
# i <- 5L
# j <- 7L
# 
# Qz(x=x, y=y, rho = rho, mu = mu, sigma = sigma, thres = thres, num_y = num_y)[i,j]
# Qij(theta = theta, i = i, j = j, x = x, y = y, num_params = num_params, num_y = num_y)




## TODO: double check Qz, always returns 0 in certain cells (probably taus that are uncorrelated)

## TODO: test score function numerically for correctness!
# once you have that, test the elementwise integral strategy by solving the expectation integral over the score at true theta, should be zero!

## tester for correctness!
# TODO: funs depend on sigma, bbut the param is theta^2! needs to be handled as below:
# THIS IS IMPORTANT FOR NUMERIC APPROX OF FISHER INFO!
# f = function(theta) py_x_scalars(y,x,theta[1],theta[2],sqrt(theta[3]),theta[4:7], num_y)
# numDeriv::grad(f, theta_true)
# py_x_grad(y,x,rho,mu,sigma,thres,num_y)
