
## chi square test for bivariate normality
chisq_test <- function(x, y, theta, Kx, Ky)
{
  
  ## prepare input
  idxX <- 2:Kx
  idxY <- (Kx+1):(Kx+Ky-1)
  rho <- theta[1]
  thresX <- c(-Inf, theta[idxX], Inf) 
  thresY <- c(-Inf, theta[idxY], Inf)
  
  ## empirical and theoretical probabilities
  f <- fhat(x = x, y = y, Kx = Kx, Ky = Ky)
  p <- model_probabilities(rho = rho, thresX = thresX, 
                           thresY = thresY, Kx = Kx, Ky = Ky)
  N <- length(x)
  counts <- f * N # empirical counts per class
  
  ## calculate the summands
  summands <- counts * log(f / p)
  teststat <- 2.0 * sum(summands[f > 0])
  
  ## degrees of freedom
  #Kxx <- length(unique(x))
  #Kyy <- length(unique(y))
  #K <- sum(f > 0) # skip empty classes
  #df <- K - Kxx - Kyy
  df <- Kx * Ky - Kx - Ky
  
  ## p value
  pval <- stats::pchisq(q = teststat, df = df, lower.tail = FALSE)
  
  return(list(chisq = teststat, df = df, pval = pval))
}


# return asymptotic covariance matrix
polycor_variance <- function(theta, c1, c2, Kx, Ky, f, N)
{
  
  ## prepare input
  idxX <- 2:Kx
  idxY <- (Kx+1):(Kx+Ky-1)
  rho <- theta[1]
  thresX <- c(-Inf, theta[idxX], Inf) 
  thresY <- c(-Inf, theta[idxY], Inf) 
  
  # calculate M, W matrices
  tmp <- get_MW(rho = rho, 
                thresX = thresX, thresY = thresY,
                f = f, c1 = c1, c2 = c2, Kx = Kx, Ky = Ky)
  W <- tmp$W
  M <- tmp$M
  
  ## calculate Lambda matrix
  Lambda <- diag(f) - outer(f,f)
  
  ## calculate final matrix
  U <- W %*% Lambda %*% t(W)
  Minv <- solve(M)
  Sigma <- Minv %*% U %*% Minv
  colnames(Sigma) <- rownames(Sigma) <- theta_names(Kx = Kx, Ky = Ky)
  
  ## Sigma is the asy. covariance matrix of \sqrt(N)*thetahat
  # thus, the asymptotic covariance matrix of thetahat is given by Sigma / N
  return(Sigma / N)
}


# MLE variance (based on inverted fisher information )
polycor_variance_mle <- function(theta, Kx, Ky, N)
{
  ## prepare input
  idxX <- 2:Kx
  idxY <- (Kx+1):(Kx+Ky-1)
  rho <- theta[1]
  thresX <- c(-Inf, theta[idxX], Inf) 
  thresY <- c(-Inf, theta[idxY], Inf) 
  
  
  ## calculate fisher information matrix
  fisher <- 
    get_fisher(rho = rho, thresX = thresX, thresY = thresY, Kx = Kx, Ky = Ky)
  
  ## calculate asymptotic covariance matrix
  Sigma <- solve(fisher)
  
  ## dimension naming
  colnames(Sigma) <- rownames(Sigma) <- theta_names(Kx = Kx, Ky = Ky)
  
  ## Sigma is the asy. covariance matrix of \sqrt(N)*thetahat
  # thus, the asymptotic covariance matrix of thetahat is given by Sigma / N
  return(Sigma / N)
}


## main functionality for two-step-approach
# NOTE: init must be scalar
polycor_twostep <- function(f, contingency, Kx, Ky, N, init, maxcor)
{
  ## frequency counts (intentionally coded as a double)
  freq <- f * N
  
  ## calculate marginals from contingency table and construct threshold estimates
  x_marg <- rowSums(contingency)
  y_marg <- colSums(contingency)
  thresx <- stats::qnorm(cumsum(x_marg[seq_len(Kx-1L)]) / N)
  thresy <- stats::qnorm(cumsum(y_marg[seq_len(Ky-1L)]) / N)
  thresX <- unname(c(-Inf, thresx, Inf))
  thresY <- unname(c(-Inf, thresy, Inf))
  
  ## define local function for the optimization
  K <- Kx * Ky
  fn <- function(rho)
  {
    objective_mle_cpp(rho = rho, freq = freq, thresX = thresX, thresY = thresY,
                      Kx = Kx, Ky = Ky, K = K, mean = c(0.0, 0.0))
  }

  ## optimize
  opt <- stats::optimize(f = fn, interval = c(-1, 1))
  
  ## extract and name estimated parameters
  thetahat <- c(opt$minimum, thresx, thresy)
  names(thetahat) <- theta_names(Kx = Kx, Ky = Ky)
  
  ## extract model parameters
  theta_obj <- extractfromtheta(theta = thetahat, Kx = Kx, Ky = Ky)
  probs <- model_probabilities(rho = theta_obj$rho, 
                               thresX = theta_obj$thresX, thresY = theta_obj$thresY,
                               Kx = Kx, Ky = Ky)
  
  ## calculate pearson residuals
  pearson <- f / probs - 1.0
  
  return(structure(
    list(thetahat = thetahat, 
         stderr = NULL,
         vcov = NULL,
         residuals = vec2tab(pearson, Kx = Kx, Ky = Ky),
         probs = vec2tab(probs, Kx = Kx, Ky = Ky),
         f = vec2tab(f, Kx = Kx, Ky = Ky),
         chisq = NULL, 
         df = NULL,
         pval = NULL,
         objective = -opt$objective, # since likelihood is MAXimized
         optim = opt, 
         inputs = list(Kx = Kx, Ky = Ky, N = N)),
    class = c("robpolycor", "polycor")))
}
  


# main functionality for polychoric correlation
# assumes that c takes values in [1,Inf] (different than notation in paper, but advantageous for computing)
polycor_fast <- 
  function(f, 
           Kx, Ky,
           c1, c2,
           N,
           init,
           method = NULL,
           maxcor = 0.999,
           tol_thresholds = 0.01,
           constrained = "ifneeded",
           variance = FALSE) 
{
  ## initializations
  is_null_method <- is.null(method)
  K <- Kx * Ky
  logc1p1 <- log(c1) + 1.0
  logc2p1 <- log(c2) + 1.0
  
  ## local functions for optim
  # theta is a vector of dimension (1 + (Kx-1) + (Ky-1)) = Kx+Ky-1
  # theta[1] is rho
  # theta[(1+1):(Kx-1+1)] = theta[2:Kx] are X-thresholds (increasing)
  # theta[(Kx+1):(Kx+Ky-1)] are Y-thresholds (increasing)
  idxX <- 2:Kx
  idxY <- (Kx+1):(Kx+Ky-1)
  fn <- function(theta)
  {
    rho <- theta[1]
    thresX <- c(-Inf, theta[idxX], Inf) # first and last thresholds are fixed to be \mp Inf
    thresY <- c(-Inf, theta[idxY], Inf) # first and last thresholds are fixed to be \mp Inf
    objective_cpp_fast(rho = rho, f = f,
                       thresX = thresX, thresY = thresY, 
                       c1 = c1, c2 = c2,
                       Kx = Kx, Ky = Ky, K = K, 
                       logc1p1 = logc1p1, logc2p1 = logc2p1,
                       mean = c(0.0, 0.0),
                       maxcor = maxcor)
  }
  
  
  if (constrained == "ifneeded") {
    ## try unconstrained optimization first and use constrained optimization if 
    ## unconstrained fails (quick & dirty, mostly copy & paste from below)
    if (is_null_method) method <- "L-BFGS-B"
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
      lincon <- constrOptim_constraints(Kx = Kx, Ky = Ky, maxcor = maxcor, 
                                        tol_thresholds = tol_thresholds)
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
    lincon <- constrOptim_constraints(Kx = Kx, Ky = Ky, maxcor = maxcor, 
                                      tol_thresholds = tol_thresholds)
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
    if (is_null_method) method <- "L-BFGS-B"
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
  names(thetahat) <- theta_names(Kx = Kx, Ky = Ky)
  
  ## extract model parameters
  theta_obj <- extractfromtheta(theta = thetahat, Kx = Kx, Ky = Ky)
  probs <- model_probabilities(rho = theta_obj$rho, 
                               thresX = theta_obj$thresX, thresY = theta_obj$thresY,
                               Kx = Kx, Ky = Ky)
  
  ## calculate pearson residuals
  pearson <- f / probs - 1.0
  
  ## if requested, estimate asymptotic covariance matrix
  if(variance)
  {
    asv <- polycor_variance(theta = thetahat, c1 = c1, c2 = c2,
                             Kx = Kx, Ky = Ky, f = f, N = N)
    stderr <- sqrt(diag(asv))

  } else
  {
    asv <- stderr <- NULL
  }
  
  ## perform chisq test for bivariate normality
  chisq <- FALSE # TODO: this test needs to be implemented also for contingency tables, so do not yet make it accessible
  x <- y <- NULL # to avoid a warning in R CMD check 
  if(chisq)
  {
    chisqtest <- chisq_test(x = x, y = y, 
                            theta = thetahat, Kx = Kx, Ky = Ky)
    chisq <- chisqtest$chisq
    df <- chisqtest$df
    pval <- chisqtest$pval
  } else
  {
    chisq <- df <- pval <- NULL
  }
  
  return(structure(
    list(thetahat = thetahat, 
         stderr = stderr,
         vcov = asv,
         residuals = vec2tab(pearson, Kx = Kx, Ky = Ky),
         probs = vec2tab(probs, Kx = Kx, Ky = Ky),
         f = vec2tab(f, Kx = Kx, Ky = Ky),
         chisq = chisq, 
         df = df,
         pval = pval,
         objective = opt$value,
         optim = opt, 
         inputs = list(Kx = Kx, Ky = Ky, N = N, c = c2 - 1.0, method = method)), # re-center c to be in [0,Inf]
    class = c("robpolycor", "polycor")))
} # FUN



#' Robust estimation of polychoric correlation 
#' 
#' Implements to robust estimator of  Welz, Mair and Alfons (2024, \doi{10.48550/arXiv.2407.18835})  for the polychoric correlation model, based on the general theory of C-estimation proposed by Welz (2024, \doi{10.48550/arXiv.2403.11954}).
#' 
#' @param x Vector of integer-valued responses to first item, or contingency table (a \code{"\link[base]{table}"} object).
#' @param y Vector of integer-valued responses to second item; only required if \code{x} is not a contingency table.
#' @param c Tuning constant that governs robustness; must be in \code{[0, Inf]}. Defaults to 0.6.
#' @param variance Shall an estimated asymptotic covariance matrix be returned? Default is \code{TRUE}.
#' @param method Numerical optimization method, see \code{\link[stats]{optim}()} and \code{\link[stats]{constrOptim}()}. Default is to use \code{"L-BFGS-B"} in case of unconstrained optimization and \code{"Nelder-Mead"} in case of constrained optimization.
#' @param constrained Shall strict monotonicity of thresholds be explicitly enforced by linear constraints? This can be a logical (\code{TRUE} or \code{FALSE}), or \code{"ifneeded"} to first try unconstrained optimization and in case of an error perform constrained optimization. Default is \code{"ifneeded"}.
#' @param maxcor Maximum absolute correlation (to ensure numerical stability). Default is 0.999.
#' @param tol_thresholds Minimum distance between consecutive thresholds (to enforce strict monotonicity); only relevant in case of constrained optimization. Default is 0.01.
#' @param init Initialization of numerical optimization. Default is neutral.
#'
#' @return 
#' An object of class \code{"robpolycor"}, which is a list with the following components. 
#' \describe{
#'   \item{\code{thetahat}}{A vector of estimates for the polychoric correlation coefficient (\code{rho}) as well as thresholds for \code{x} (named \code{a1,a2,...,a_{Kx-1}}) and \code{y} (named \code{b1,b2,...,b_{Ky-1}}).}
#'   \item{\code{stderr}}{A vector of standard errors for each estimate in \code{thetahat}.}
#'   \item{\code{vcov}}{Estimated asymptotic covariance matrix of \code{thetahat}. The matrix \eqn{\Sigma} in the paper (asymptotic covariance matrix of \eqn{\sqrt{N} \hat{\theta}}) can be obtained via \code{vcov * N}, where \code{N} is the sample size.}
#'   \item{\code{chisq,pval,df}}{Currently \code{NULL}, will in a future release be the test statistic, p-value, and degrees of freedom of a test for bivariate normality.}
#'   \item{\code{objective}}{Value of minimized loss function.}
#'   \item{\code{optim}}{Object of class \code{optim}.}
#' }
#' 
#' @examples
#' ## example data
#' set.seed(123)
#' x <- sample(c(1,2,3), size = 100, replace = TRUE)
#' y <- sample(c(1,2,3), size = 100, replace = TRUE)
#' 
#' polycor(x,y)     # robust
#' polycor_mle(x,y) # non-robust MLE
#' 
#' @export
polycor <- function(x, y = NULL, c = 0.6, 
                    variance = TRUE,
                    constrained = "ifneeded",
                    method = NULL,
                    maxcor = 0.999,
                    tol_thresholds = 0.01,
                    init = initialize_param(x, y))
{
  ## initializations
  stopifnot(c >= 0)
  if (is.character(constrained)) 
  {
    constrained <- match.arg(constrained)
  } else 
  {
    constrained <- isTRUE(constrained)
  }
  if(is.table(x))
  {
    inputs <- input_table(x = x, y = y)
  } else
  {
    inputs <- input_vector(x = x, y = y)
  }
  if (sum(diag(inputs$contingency)) == inputs$N) {
    warning("observed variables are identical; polychoric model is not identified")
  }
  
  ## add 1 so that c is in [1,Inf] to comply with expectation of polycor_fast()
  c_reparam <- c + 1.0
  
  polycor_fast(f = inputs$f, Kx = inputs$Kx, Ky = inputs$Ky, N = inputs$N,
               c1 = 0.0, c2 = c_reparam, method = method, maxcor = maxcor, 
               tol_thresholds = tol_thresholds, constrained = constrained,
               init = init, variance = variance)
}


#' Maximum likelihood estimation of polychoric correlation coefficient
#' 
#' Implements the maximum likelihood estimator of Olsson (1979, Psychometrika, \doi{10.1007/BF02296207}) for the polychoric correlation model.
#' 
#' @inheritParams polycor
#' @param twostep Shall two-step estimation of Olsson (1979) <doi:10.1007/BF02296207> be performed? Default is \code{FALSE}.
#' 
#' @return An object of class \code{"robpolycor"}. See \code{\link{polycor}()} for details.
#' 
#' @examples
#' ## example data
#' set.seed(123)
#' x <- sample(c(1,2,3), size = 100, replace = TRUE)
#' y <- sample(c(1,2,3), size = 100, replace = TRUE)
#' 
#' polycor(x,y)     # robust
#' polycor_mle(x,y) # non-robust MLE
#' 
#' @export
polycor_mle <- function(x, y = NULL, 
                        variance = TRUE,
                        constrained = "ifneeded",
                        twostep = FALSE,
                        method = NULL,
                        maxcor = 0.999,
                        tol_thresholds = 0.01,
                        init = initialize_param(x, y))
{
  ## initializations
  if (is.character(constrained)) 
  {
    constrained <- match.arg(constrained)
  } else 
  {
    constrained <- isTRUE(constrained)
  }
  if(is.table(x))
  {
    inputs <- input_table(x = x, y = y)
  } else
  {
    inputs <- input_vector(x = x, y = y)
  }
  if (sum(diag(inputs$contingency)) == inputs$N) {
    warning("observed variables are identical; polychoric model is not identified")
  }
  
  if(twostep)
  {
    obj <- 
      polycor_twostep(f = inputs$f, contingency = inputs$contingency, 
                      Kx = inputs$Kx, Ky = inputs$Ky, N = inputs$N, 
                      init = init[1L], maxcor = maxcor)
  } else
  {
    # TODO: adjust to MLE loss function 
    obj <- 
      polycor_fast(f = inputs$f, Kx = inputs$Kx, Ky = inputs$Ky, N = inputs$N,
                   c1 = 0.0, c2 = Inf, method = method, maxcor = maxcor, 
                   tol_thresholds = tol_thresholds, constrained = constrained,
                   init = init, variance = FALSE)
  }
  
  
  # if variance requested, calculate it based on fisher information
  if(variance)
  {
    asv <- 
      polycor_variance_mle(theta = obj$thetahat, Kx = inputs$Kx, Ky = inputs$Ky, N = inputs$N)
    stderr <- sqrt(diag(asv))

    # if(twostep)
    # {
    #   # in twostep, variance of thresholds isn't estimated consistently, so drop
    #   stderr[2:length(stderr)] <- NA_real_
    #   asv <- sigma[1L,1L,drop = FALSE]
    # }
    
    obj$stderr <- stderr
    obj$vcov <- asv
  } # IF variance 
  
  return(structure(obj, class = c("robpolycor", "robpolycor_mle", "polycor")))
}
