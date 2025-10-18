#' Obtain estimated asymptotic variance-covariance matrix 
#' 
#' Estimate asymptotic variance-covariance matrix of polychoric model.
#' 
#' Method for classes \code{"robpolycor"} and \code{"polycor"}. Returns the estimated asymptotic variance-covariance matrix of a point estimate \code{thetahat}. The matrix \eqn{\Sigma} in the paper (asymptotic variance-covariance matrix of \eqn{\sqrt{N} \hat{\theta}}) can be obtained via multiplying the returned matrix by the sample size.
#' 
#' @param object Object of class \code{"robpolycor"} or \code{"polycor"}.
#' @param ... Additional parameters to be passed down.
#' @return A numeric matrix, being the estimated asymptotic covariance matrix for the model parameters
#' @examples
#' set.seed(123)
#' x <- sample(c(1,2,3), size = 100, replace = TRUE)
#' y <- sample(c(1,2,3), size = 100, replace = TRUE)
#' fit <- polycor(x,y) 
#' 
#' vcov(fit)
#' 
#' @export
vcov.robpolycor <- function(object, ...)
{
  stopifnot(inherits(x = object, what = "polycor"))
  
  fhat     <- as.numeric(t(object$f)) ## transpose to ensure expected order
  N        <- object$inputs$N
  Kx       <- object$inputs$Kx
  Ky       <- object$inputs$Ky
  thetahat <- object$thetahat
  
  if(inherits(x = object, what = "robpolycor_mle"))
  {
    ## for MLE, compute fisher-information-based covariance matrix
    out <- polycor_variance_mle(theta = thetahat, Kx = Kx, Ky = Ky, N = N)
  } else
  {
    out <- polycor_variance(theta = thetahat, c1 = -Inf, 
                            c2 = object$inputs$c + 1.0, ## internally, c is in range [1, Inf], so add 1
                            Kx = Kx, Ky = Ky, f = fhat, N = N)
  }
  return(out)
} 


#' Obtain estimated asymptotic variance-covariance matrix 
#' 
#' Estimate asymptotic variance-covariance matrix of polyserial model.
#' 
#' Method for classes \code{"robpolyserial"} and \code{"polyserial"}. Returns the estimated asymptotic variance-covariance matrix of a point estimate \code{thetahat}. The matrix \eqn{\Sigma} in the paper (asymptotic variance-covariance matrix of \eqn{\sqrt{N} \hat{\theta}}) can be obtained via multiplying the returned matrix by the sample size.
#' 
#' @param object Object of class \code{"robpolyserial"} or \code{"polyserial"}.
#' @param ... Additional parameters to be passed down.
#' @return A numeric matrix, being the estimated asymptotic covariance matrix for the model parameters
#' @examples
#' ## example data
#' set.seed(123)
#' x <- rnorm(n = 100)
#' y <- sample(c(1,2), size = 100, replace = TRUE)
#' 
#' fit <- polyserial(x,y)
#' vcov(fit)
#' 
#' @export
vcov.robpolyserial <- function(object, ...)
{
  stopifnot(inherits(x = object, what = "polyserial"))
  
  thetahat <- object$thetahat
  x        <- object$inputs$x
  y        <- object$inputs$y

  if(inherits(x = object, what = "robpolyserial_mle"))
  {
    ## for MLE, compute fisher-information-based covariance matrix
    out <- variance_polyserial_mle(theta = thetahat, x = x, y = y)
  } else
  {
    out <-  variance_polyserial(theta = thetahat, x = x, y = y, alpha = object$inputs$alpha)
  }
  return(out$vcov)
} 
