#' perform test if individual cells are outlying
#' 
#' @param x object of class \code{polycor}
#' @param twosided shall test be one- or two-sided with right-tail alternative?
#' @param adjust method for adjustment of p-values for multiple comparisons; default is FDR
#' @param ... additional parameters to be passed down
#' 
#' @export
celltest <- function(x, twosided = FALSE, adjust = "fdr", ...)
{
  stopifnot(inherits(x, what = "polycor"))
  Kx <- x$inputs$Kx ; Ky <- x$inputs$Ky ; N <- x$inputs$N
  thetaobj <- extractfromtheta(theta = x$thetahat, Kx = Kx, Ky = Ky, includeInf = TRUE)
  obj <- 
    celltest_cpp(rho = thetaobj$rho, 
                 thresX = thetaobj$thresX,
                 thresY = thetaobj$thresY,
                 Kx = Kx, Ky = Ky,  
                 probs = as.matrix(x$probs), 
                 f = as.matrix(x$f), 
                 sigma = x$sigma, N = N, 
                 twosided = twosided)
  
  pval <- matrix(stats::p.adjust(as.numeric(obj$pval), method = adjust, ...), 
                 nrow = Kx, ncol = Ky, byrow = FALSE)
  
  return(
    list(
      teststat = mat2tab(obj$teststat),
      stderr = mat2tab(obj$stderr),
      pval_adjusted = mat2tab(pval),
      pval_raw = mat2tab(obj$pval) ))
}
