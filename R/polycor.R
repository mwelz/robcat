# don't forget that you must manually change the NAMESPACE file to register functions & methods!
# TODO: allow for data table argument like in polycor package
# TODO: use objective_cpp() in polycor function and merge tertracor() with polycor()
# TODO: maybe profile to see where the bottlenecks are (prolly the integral)
# Rcpp::as<T>();
# This C++ package has an implementation of linearly constrained optimization. Ask AA if we're allowed to use header files of it (and if we even need it to begin with)
# TODO: implement feps() function like in the old R code; requires that 'mean' can be passed to model probs function DONE
# TODO: Ask AA if we can drop monotonicity constraint (like fox); I don't think we should
# TODO: make polycor() handle classes that are empirically zero
# perhaps implement MLE variance explicitly? (fisher-info-based?)
# TODO: chi square test for bivariate normality (see here: https://github.com/cran/polycor/blob/abb5d7b2ecc08ff282f2ee9dfeb89bce051229d9/R/polychor.R#L116). Implemented it, but I'm not sure how sensible this test is because it relies on cell empirical counts which are contaminated if there are outliers (so problem is in empirics, not the robustly estimated probabilities)
# TODO: ask AA how to incorporate inference when there is contamination: start with rho=0? correct sign?
# table(x, y) / n == matrix(fhat(x,y),5,5,byrow = T): clever trick to recover cell names!

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
  return(list(variance = Sigma / N, Sigma = Sigma))
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
  return(list(variance = Sigma / N, Sigma = Sigma))
}


# main functionality for polychoric correlation
polycor_fast <- 
  function(f, 
           Kx, Ky,
           c1, c2,
           N,
           init,
           method = "Nelder-Mead",
           tol = 0.001,
           constrained = TRUE,
           variance = FALSE) 
{
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
  feasible <- 1.0 - tol
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
                       feasible = feasible)
  }
  
  if(constrained)
  {
    ## the linear constraints
    lincon <- constrOptim_constraints(Kx = Kx, Ky = Ky, tol_rho = tol, tol_thresholds = tol)
    
    ## linearly constrained optimization
    opt <-
      stats::constrOptim(theta = init,
                         f = fn,
                         grad = NULL,
                         ui = lincon$ui,
                         ci = lincon$ci,
                         method = method,
                         hessian = FALSE)
  } else{
    ## without monotonicity constraints: runs faster, but less accurate
    opt <- 
      stats::optim(par = init, fn = fn, gr = NULL, 
                   lower = c(-1.+tol, rep(-Inf, length(Kx+Ky-2))),
                   upper = c(1.-tol, rep(Inf, length(Kx+Ky-2))), 
                   hessian = FALSE, method = method)
  } # IF
  
  ## extract and name estimated parameters
  thetahat <- opt$par
  names(thetahat) <- theta_names(Kx = Kx, Ky = Ky)
  
  ## extract model parameters
  theta_obj <- extractfromtheta(theta = thetahat, Kx = Kx, Ky = Ky)
  probs <- model_probabilities(rho = theta_obj$rho, 
                               thresX = theta_obj$thresX, thresY = theta_obj$thresY,
                               Kx = Kx, Ky = Ky)
  
  ## calculate pearson residuals
  pearson <- f / probs
  
  ## if requested, estimate covariance matrix
  if(variance)
  {
    tmp <- polycor_variance(theta = thetahat, c1 = c1, c2 = c2,
                            Kx = Kx, Ky = Ky, f = f, N = N)
    stderr <- sqrt(diag(tmp$variance))
    Sigma  <- tmp$Sigma
    
  } else
  {
    Sigma <- stderr <- NULL
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
         sigma = Sigma,
         residuals = vec2tab(pearson, Kx = Kx, Ky = Ky),
         probs = vec2tab(probs, Kx = Kx, Ky = Ky),
         f = vec2tab(f, Kx = Kx, Ky = Ky),
         chisq = chisq, 
         df = df,
         pval = pval,
         objective = opt$value,
         optim = opt, 
         inputs = list(Kx = Kx, Ky = Ky, N = N, c = c2)),
    class = c("robpolycor", "polycor")))
} # FUN



#' Robust estimation of polychoric correlation coefficient
#' 
#' @param x vector of integer-valued responses to first item or contingency table (a \code{table} object)
#' @param y vector of integer-valued responses to second item; only required if \code{x} is not a contingency table 
#' @param c tuning constant that governs robustness; defaults to 1.6
#' @param variance shall an estimated asymptotic covariance matrix be returned? Default is \code{TRUE}
#' @param method numerical optimization method
#' @param constrained shall strict monotonicity of thresholds be explicitly enforced by linear constraints? 
#' @param tol tolerance in numerical optimization
#' @param init initialization of numerical optimization. Default is neutral
#'
#' @return 
#' A list with the following components. 
#' \describe{
#'   \item{\code{theahat}}{A vector of estimates for the polychoric correlation coefficient (\code{rho}) as well as thresholds for \code{x} (named \code{a1,a2,...,a_{Kx-1}}) and \code{y} (named \code{b1,b2,...,b_{Ky-1}}).}
#'   \item{\code{stderr}}{A vector of standard errors for each estimate in \code{theahat}}
#'   \item{\code{sigma}}{Estimated asymptotic covariance matrix \eqn{\Sigma}, evaluated at the estimates in \code{theahat}.}
#'   \item{\code{chisq,pval,df}}{Test statistic, p-value, and degrees of freedom of a test for bivariate normality}
#'   \item{\code{objective}}{Value of minimized loss function}
#'   \item{\code{optim}}{Object of class \code{optim}}
#' }
#' @export
polycor <- function(x, y = NULL, c = 1.6, 
                    variance = TRUE,
                    constrained = TRUE,
                    method = ifelse(constrained, "Nelder-Mead", "L-BFGS-B"),
                    tol = 0.001,
                    init = initialize_param(x, y))
{
  stopifnot(c >= 1)
  if(is.table(x))
  {
    inputs <- input_table(x = x, y = y)
  } else
  {
    inputs <- input_vector(x = x, y = y)
  }
  
  polycor_fast(f = inputs$f, Kx = inputs$Kx, Ky = inputs$Ky, N = inputs$N,
               c1 = 0.0, c2 = c, 
               method = method, tol = tol, constrained = constrained,
               init = init, 
               variance = variance)
}


#' Maximum likelihood estimation of polychoric correlation coefficient
#' 
#' @param x vector of integer-valued responses to first item or contingency table (a \code{table} object)
#' @param y vector of integer-valued responses to second item; only required if \code{x} is not a contingency table 
#' @param variance shall an estimated asymptotic covariance matrix be returned? Default is \code{TRUE}
#' @param method numerical optimization method; default is Nelder-Mead
#' @param constrained shall strict monotonicity of thresholds be explicitly enforced by linear constraints? 
#' @param tol tolerance in numerical optimization
#' @param init initialization of numerical optimization. Default is neutral
#' @export
polycor_mle <- function(x, y = NULL, 
                        variance = TRUE,
                        constrained = TRUE,
                        method = ifelse(constrained, "Nelder-Mead", "L-BFGS-B"),
                        tol = 0.001,
                        init = initialize_param(x, y))
{
  if(is.table(x))
  {
    inputs <- input_table(x = x, y = y)
  } else
  {
    inputs <- input_vector(x = x, y = y)
  }
  
  obj <- 
    polycor_fast(f = inputs$f, Kx = inputs$Kx, Ky = inputs$Ky, N = inputs$N,
                 c1 = 0.0, c2 = Inf, 
                 method = method, tol = tol, constrained = constrained,
                 init = init, 
                 variance = FALSE)
  
  # if variance requested, calculate it based on fisher information
  if(variance)
  {
    tmp <- 
      polycor_variance_mle(theta = obj$thetahat, Kx = inputs$Kx, Ky = inputs$Ky, N = inputs$N)
    obj$stderr <- sqrt(diag(tmp$variance))
    obj$sigma <- tmp$Sigma
  } 
  
  return(obj)
}
