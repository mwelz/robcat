
## x here is a vector of length N
px <- function(x, mu, sigma) stats::dnorm(x = x, mean = mu, sd = sigma)


## calculate taustar for entire support of y (except last value) for a given scalar x
# thres only has finite values
# x must be a scalar here!!
# -> this will be a bottleneck, candidate for Rcpp!
get_taustar <- function(x, rho, mu, sigma, thres)
{
  numer <- thres - rho * (x - mu) / sigma 
  denom <- sqrt(1.0 - rho^2)
  taustar <- numer / denom # length equals r - 1
  return(taustar)
}


# get PMF of entire support of (Y | X) for a given scalar value of x
# unlike robcat, thres here only has finite values
py_x <- function(x, rho, mu, sigma, thres)
{
  taustar <- get_taustar(
    x = x, rho = rho, mu = mu, sigma = sigma, 
    thres = thres)
  
  Phi <- stats::pnorm(q = taustar, mean = 0, sd = 1)
  return(c(Phi, 1.0) - c(0.0, Phi))
}

## sum over PMF^alphaplus1 of entire support of (Y | X) for a given vector of x
sum_py_x_alpha <- Vectorize(
  function(x, rho, mu, sigma, thres, alphaplus1)
  {
    sum(py_x(x = x, rho = rho, mu = mu, sigma = sigma, thres = thres)^alphaplus1)
  }, vectorize.args = "x")


## the function we integrate over in int_ptheta()
pz_int <- function(x, rho, mu, sigma, thres, alphaplus1)
{
  part1 <- px(x = x, mu = mu, sigma = sigma)^alphaplus1
  part2 <- sum_py_x_alpha(x = x, rho = rho, mu = mu, sigma = sigma, thres = thres, alphaplus1)
  
  return(part1 * part2)
} # FUN


## integrate over p_theta(z)^(alpha+1)
integrate_ptheta <- function(rho, mu, sigma, thres, alphaplus1)
{
  fun <- function(x){
    pz_int(x = x, rho = rho, mu = mu, sigma = sigma, thres = thres, alphaplus1 = alphaplus1)
  }
  
  ## integrate
  integrated <- stats::integrate(f = fun, lower = -Inf, upper = Inf)
  
  return(integrated$value)
} # FUN



## evaluate density of (Y|X) when BOTH x and y are scalars
# FIXME: rewrite in C++!
py_x_scalars <- function(y, x, rho, mu, sigma, thres, num_y = length(thres) + 1L)
{
  if(abs(rho) >= 1){return(0)} ## pathological case where rho violates a boundary condition
  
  ## calculate tau_star
  summand <- rho * (x - mu) / sigma
  denom   <- sqrt(1.0 - rho^2)
  
  if(y == 1L)
  {
    taustar_y <- (thres[y] - summand) / denom
    out <- stats::pnorm(taustar_y)
  } else if(y == num_y)
  {
    taustar_ym1 <- (thres[y-1L] - summand) / denom
    out <- 1.0 - stats::pnorm(taustar_ym1)
  } else
  {
    taustar_y   <- (thres[y] - summand) / denom
    taustar_ym1 <- (thres[y-1L] - summand) / denom
    out <- stats::pnorm(taustar_y) - stats::pnorm(taustar_ym1)
  } # IF
  
  return(out)
} # FUN


## the above function where both y and x can be vectors
py_x_vectorized <- Vectorize(py_x_scalars, vectorize.args = c("y", "x"))


## marginal probabilities of Y
py <- function(thres)
{
  p <- stats::pnorm(q = thres, mean = 0, sd = 1)
  return(c(p, 1) - c(0, p))
}

## \sum_{i=1} (p_theta(X_i, Y_i))^alpha
# both x and y are vectors here
# FIXME: bottleneck; rewrite in C++ (may require scalar version of px())
sum_pz_alpha <- function(x, y, rho, mu, sigma, thres, num_y, alpha)
{
  ## p(X_i) for all X_i
  px_vec <- px(x = x, mu = mu, sigma = sigma)
  
  ## p(Y_i | X_i) for all i 
  py_x_vec <- py_x_vectorized(
    y = y, x = x, rho = rho, mu = mu, sigma = sigma, thres = thres, num_y = num_y)
  
  # p(Z_i) for all i
  pz <- px_vec * py_x_vec
  
  ## sum over all p(Z_i) to the alpha
  return(sum(pz^alpha))
} # FUN


## density power divergence btw fhat and ptheta
# const = (1 + alpha_inv) / N
DPD <- function(x, y, 
                rho, mu, sigma, thres,
                alpha, 
                alphaplus1 = alpha + 1,
                alpha_inv = 1 / alpha,
                N = length(x),
                const = (1 + alpha_inv) / N,
                num_y = length(thres) + 1L)
{
  
  ## summand 1: integral over p(x)^alphap1
  sum1 <- integrate_ptheta(
    rho = rho, mu = mu, sigma = sigma, thres = thres, alphaplus1 = alphaplus1)
  
  ## summand 2: sum over the joint densities at observations
  sums <- 
    sum_pz_alpha(x = x, y = y, rho = rho, mu = mu, sigma = sigma, 
                 thres = thres, num_y = num_y, alpha = alpha)
  sum2 <- const * sums
  
  ## summand 3: 1 / alpha
  sum3 <- alpha_inv
  
  return(sum1 - sum2 + sum3)
} # FUN


polyserial_thetanames <- function(num_y)
{
  c("rho", "mu", "sigma2", paste0("tau", seq_len(num_y - 1L)))
}



#' Robust estimation of polyserial correlation 
#' 
#' Implements the robust estimator of  Welz (2025)  for the polyserial correlation model.
#' 
#' @param x Vector of numeric values.
#' @param y Vector of integer-valued ordinal values.
#' @param alpha Tuning constant that governs robustness-efficiency tradeoff; must be in \code{[0, Inf]}. Defaults to 0.5.
#' @param num_y Number of response categories in y; defaults to max(y)
#' @param variance Shall an estimated asymptotic covariance matrix be returned? Default is \code{TRUE}.
#' @param method Numerical optimization method, see \code{\link[stats]{optim}()} and \code{\link[stats]{constrOptim}()}. Default is to use \code{"BFGS"} in case of unconstrained optimization and \code{"Nelder-Mead"} in case of constrained optimization.
#' @param constrained Shall parameter restructions be enforced by linear constraints? This can be a logical (\code{TRUE} or \code{FALSE}), or \code{"ifneeded"} to first try unconstrained optimization and in case of an error perform constrained optimization. Default is \code{"ifneeded"}.
#' @param maxcor Maximum absolute correlation (to ensure numerical stability). Default is 0.999.
#' @param tol_thresholds Minimum distance between consecutive thresholds (to enforce strict monotonicity); only relevant in case of constrained optimization. Default is 0.01.
#' @param tol_sigma2 Minimum value of sigma2 parameter (population variance of X); only relevant in case of constrained optimization. Default is 0.01.
#' @param init Initialization of numerical optimization. Default is neutral.
#'
#' @return 
#' An object of class \code{"polyserial"}, which is a list with the following components. 
#' \describe{
#'   \item{\code{thetahat}}{A vector of estimates for the polyserial correlation coefficient (\code{rho}), population mean of X (\code{mu}), population variance of Y (\code{sigma2}), as well as thresholds for \code{y} (named \code{tau1,tau2,...,tau_{r-1}}).}
#'   \item{\code{stderr}}{A vector of standard errors for each estimate in \code{thetahat}.}
#'   \item{\code{vcov}}{Estimated asymptotic covariance matrix of \code{thetahat}. The matrix \eqn{\Sigma} in the paper (asymptotic covariance matrix of \eqn{\sqrt{N} \hat{\theta}}) can be obtained via \code{vcov * N}, where \code{N} is the sample size.}
#'   \item{\code{pointpolyserial}}{Estimated polyserial correlation coefficient, calculated with provided scoring of Y}
#'   \item{\code{weights}}{List of rescaled and raw outlyingness weights for each observation as well as maximum possible raw weight that was used for rescaling (\code{sup}).}
#'   \item{\code{objective}}{Value of minimized loss function.}
#'   \item{\code{optim}}{Object of class \code{optim}.}
#'   \item{\code{inputs}}{List of provided inputs.}
#' }
#' 
#' @examples
#' ## example data
#' set.seed(123)
#' x <- rnorm(n = 100)
#' y <- sample(c(1,2), size = 100, replace = TRUE)
#' 
#' polyserial(x,y)
#' 
#' @export
polyserial <- function(x, y, 
                       alpha = 0.5, 
                       num_y = max(y),
                       constrained = "ifneeded",
                       method = NULL,
                       variance = TRUE,
                       init = polyserial_initialize_param(x = x, num_y = num_y, robust = TRUE),
                       maxcor = 0.999,
                       tol_thresholds = 0.01,
                       tol_sigma2 = 0.01)
{
  ## important constants
  N          <- length(x)
  alphaplus1 <- 1.0 + alpha
  alpha_inv  <- 1.0 / alpha
  const      <- (1.0 + alpha_inv) / N
  stopifnot(N == length(y))
  stopifnot(alpha > 0.0)
  check_NA(x) ; check_NA(y)
  is_numeric(x, checkinteger = FALSE)
  is_numeric(y, checkinteger = TRUE)
  
  ## initialize params
  init <- polyserial_initialize_param(x = x, num_y = num_y, robust = TRUE)
  num_params <- num_y + 2L
  is_null_method <- is.null(method)
  
  ## prepare function to optimize over
  fn <- function(theta)
  {
    
    ## extract the parameter values
    rho    <- theta[1L]
    mu     <- theta[2L]
    sigma2 <- theta[3L]
    thres  <- theta[4L:num_params]
    
    DPD(x = x, y = y, rho = rho, mu = mu, sigma = sqrt(sigma2), thres = thres,
        alpha = alpha, alphaplus1 = alphaplus1, alpha_inv = alpha_inv,
        N = N, const = const, num_y = num_y)
  } 
  
  
  if (constrained == "ifneeded") {
    ## try unconstrained optimization first and use constrained optimization if 
    ## unconstrained fails (quick & dirty, mostly copy & paste from below)
    if (is_null_method) method <- "BFGS"
    opt <- try(
      stats::optim(par = init, fn = fn, gr = NULL, 
                   method = method, hessian = FALSE), 
      silent = TRUE
    )
    # TODO: Should we also perform constrained optimization if the solution is 
    #       identical to starting values? I'm not sure since the issue may also 
    #       be solved with a different (unconstrained) optimization method, and 
    #       we give a corresponding warning below.
    # if (inherits(opt, "try-error") || identical(init, opt$par)) {
    if (inherits(opt, "try-error")) {
      if (is_null_method) method <- "Nelder-Mead"
      lincon <- polyserial_constrOptim_constraints(
        num_cat = num_y, maxcor = maxcor,
        tol_thresholds = tol_thresholds,
        tol_sigma2 = tol_sigma2)
      opt <- 
        stats::constrOptim(theta = init,
                           f = fn,
                           grad = NULL,
                           ui = lincon$ui,
                           ci = lincon$ci,
                           method = method,
                           hessian = FALSE)
    }
  } else if (constrained) {
    ## set default optimization method if necessary
    if (is_null_method) method <- "Nelder-Mead"
    ## construct the linear constraints
    lincon <- polyserial_constrOptim_constraints(
      num_cat = num_y, maxcor = maxcor,
      tol_thresholds = tol_thresholds,
      tol_sigma2 = tol_sigma2)
    ## linearly constrained optimization
    opt <-
      stats::constrOptim(theta = init,
                         f = fn,
                         grad = NULL,
                         ui = lincon$ui,
                         ci = lincon$ci,
                         method = method,
                         hessian = FALSE)
  } else {
    ## set default optimization method if necessary
    if (is_null_method) method <- "BFGS"
    ## without monotonicity or rho-boundary constraints: 
    ## runs faster, but less accurate
    opt <- 
      stats::optim(par = init, fn = fn, gr = NULL, 
                   method = method, hessian = FALSE)
  }
  
  ## extract and name estimated parameters
  thetahat <- opt$par
  if (identical(init, thetahat)) {
    warning("estimates are identical to starting values; ", 
            "try a different optimization method or different starting values")
  }
  names(thetahat) <- polyserial_thetanames(num_y)
  
  
  ## calculate weights
  wgt_ls <- weights(x = x, y = y, theta = thetahat, alpha = alpha)
  
  ## point polyserial correlation
  pp <- pointpolyserial(theta = thetahat, ycat = sort(unique(y), decreasing = FALSE))
  
  ## if requested, calculate asymptotic covariance matrix
  if(variance)
  {
    asv <- variance_polyserial(theta = thetahat, x = x, y = y, alpha = alpha)
  } else
  {
    asv <- NULL
  } # IF
  
  return(structure(
    list(thetahat = thetahat, 
         stderr = asv$se,
         vcov = asv$vcov,
         pointpolyserial = pp,
         weights = list(
           rescaled = wgt_ls$weights_rescaled, 
           raw = wgt_ls$weights_raw,
           sup = wgt_ls$sup_weight),
         objective = opt$value,
         optim = opt, 
         inputs = list(
           x = x, y = y, num_y = num_y, N = N, alpha = alpha, method = method)), 
    class = c("robpolyserial", "polyserial")))
} # FUN


### functions for MLE
# log-likelihood
logl <- function(x, y, rho, mu, sigma, thres, num_y)
{
  marg <- px(x = x, mu = mu, sigma = sigma)
  cond <- py_x_vectorized(y = y, x = x, rho = rho, mu = mu, sigma = sigma, thres = thres, num_y = num_y)
  logprob <- log(marg * cond)
  
  return(sum(logprob))
}


polyserial_mle <- function(x, y, num_y = max(y),
                           constrained = "ifneeded",
                           method = NULL,
                           maxcor = 0.999,
                           tol_thresholds = 0.01,
                           tol_sigma2 = 0.01)
{
  ## initialize params
  # the problem is convex, so initialization shouldn't matter. However, if x has at least one extreme 
  # value, it may happen that px(x) is exactly zero somewhere so that log(x) is -Inf somewhere
  # Capping doesn't work here because the loss becomes constant in the the infinite dimension,
  # so it won't matter for optimization, thereby negating the adverse effect of the outlier
  # Thus, we need to avoid that loss becomes infinite anywhere (esp. starting value), which
  # we can achieve by nonrobust initialization. Again, this is a computational problem
  # because in theory, initialization shouldn't matter due to convexity. 
  # Also, nonrobust initialization is equivalent to what the Fox package does (they even fix the moments like that)
  init <- polyserial_initialize_param(x = x, num_y = num_y, robust = FALSE)
  num_params <- num_y + 2L
  N <- length(x)
  is_null_method <- is.null(method)
  
  ## prepare function to minimize over
  # the negative log-likelihood here
  fn <- function(theta)
  {
    
    ## extract the parameter values
    rho    <- theta[1L]
    mu     <- theta[2L]
    sigma2 <- theta[3L]
    thres  <- theta[4L:num_params]
    
    -logl(x = x, y = y, rho = rho, mu = mu, sigma = sqrt(sigma2), thres = thres, num_y = num_y)
  } 
  
  
  if (constrained == "ifneeded") {
    ## try unconstrained optimization first and use constrained optimization if 
    ## unconstrained fails (quick & dirty, mostly copy & paste from below)
    if (is_null_method) method <- "BFGS"
    opt <- try(
      stats::optim(par = init, fn = fn, gr = NULL, 
                   method = method, hessian = FALSE), 
      silent = TRUE
    )
    # TODO: Should we also perform constrained optimization if the solution is 
    #       identical to starting values? I'm not sure since the issue may also 
    #       be solved with a different (unconstrained) optimization method, and 
    #       we give a corresponding warning below.
    # if (inherits(opt, "try-error") || identical(init, opt$par)) {
    if (inherits(opt, "try-error")) {
      if (is_null_method) method <- "Nelder-Mead"
      lincon <- polyserial_constrOptim_constraints(
        num_cat = num_y, maxcor = maxcor,
        tol_thresholds = tol_thresholds,
        tol_sigma2 = tol_sigma2)
      opt <- 
        stats::constrOptim(theta = init,
                           f = fn,
                           grad = NULL,
                           ui = lincon$ui,
                           ci = lincon$ci,
                           method = method,
                           hessian = FALSE)
    }
  } else if (constrained) {
    ## set default optimization method if necessary
    if (is_null_method) method <- "Nelder-Mead"
    ## construct the linear constraints
    lincon <- polyserial_constrOptim_constraints(
      num_cat = num_y, maxcor = maxcor,
      tol_thresholds = tol_thresholds,
      tol_sigma2 = tol_sigma2)
    ## linearly constrained optimization
    opt <-
      stats::constrOptim(theta = init,
                         f = fn,
                         grad = NULL,
                         ui = lincon$ui,
                         ci = lincon$ci,
                         method = method,
                         hessian = FALSE)
  } else {
    ## set default optimization method if necessary
    if (is_null_method) method <- "BFGS"
    ## without monotonicity or rho-boundary constraints: 
    ## runs faster, but less accurate
    opt <- 
      stats::optim(par = init, fn = fn, gr = NULL, 
                   method = method, hessian = FALSE)
  }
  
  ## extract and name estimated parameters
  thetahat <- opt$par
  if (identical(init, thetahat)) {
    warning("estimates are identical to starting values; ", 
            "try a different optimization method or different starting values")
  }
  names(thetahat) <- polyserial_thetanames(num_y)
  
  ## value of objective function
  # negative logL divided by N. Not directly comparable
  # to DPD loss, but useful to spot identification issues
  objective <- opt$value / N
  
  return(list(thetahat = thetahat, objective = objective, optim = opt))
}


#' Neutral initialization of starting values for polyserial correlation
#' 
#' Initializes starting values for numerical optimization in a neutral way.
#' 
#' @param x Vector of numeric values
#' @param num_y Number of response options of ordinal variable
#' @param robust Should values of \code{mu} and \code{sigma2} be initialized in robust way (that is, median and squared MAD)? If \code{FALSE} (default), then nonrobust sample mean and variance are used.
#' @return A vector of initial values for the polyserial correlation coefficient, mu, sigma2, and Y-threshold parameters
#' @examples
#' ## example data
#' set.seed(123)
#' x <- rnorm(100)
#' polyserial_initialize_param(x = x, num_y = 3)
#' 
#' @export
polyserial_initialize_param <- function(x, num_y, robust = FALSE)
{
  # order: rho, mu, sigma2, thres
  init_thres <- 1L:(num_y - 1L) - stats::median(1:(num_y - 1))
  
  if(robust)
  {
    mu <- stats::median(x)
    sigma2 <- stats::mad(x)^2
  } else
  {
    mu <- mean(x)
    sigma2 <- stats::var(x)
  }
  
  init <- c(0.0, mu, sigma2, init_thres)
  return(init)
}


## for constrained optimization
polyserial_constrOptim_constraints <- function(num_cat, maxcor = 0.999, tol_thresholds = 0.001, tol_sigma2 = 0.001)
{
  num_param <- num_cat + 2L
  num_constr <- num_cat + 1L

  # 2 constraints on rho (in (-1,1)), 0 on mu, 1 on sigma2 (>0), r-2 on thresholds (monotonicity)
  U_rho <- matrix(c(1L, -1L), nrow = 2L, ncol = 1L)
  U_sigma2 <- matrix(1L, nrow = 1L, ncol = 1L)
  U_rhomu <- cbind(U_rho, 0L)
  
  if(num_cat > 2L)
  {
    U_tau <- constrOptim_submatrix(K = num_cat)
    U <- as.matrix(Matrix::bdiag(U_rhomu, U_sigma2, U_tau)) # no constraints on mu!
    C <- c(-maxcor, -maxcor, tol_sigma2, rep(tol_thresholds, num_cat-2L))
  } else
  {
    ## if Y is dichotomous, no constraints needed on thresholds
    U <- as.matrix(Matrix::bdiag(U_rhomu, U_sigma2)) # no constraints on mu!
    C <- c(-maxcor, -maxcor, tol_sigma2)
  }
  
  return(list(ui = U, ci = C))
}
