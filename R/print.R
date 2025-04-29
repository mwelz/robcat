#' Print method for classes \code{"robpolycor"} and \code{"polycor"}.
#' 
#' @param x Object of class \code{"robpolycor"} or \code{"polycor"}.
#' @param digits Number of digits to be printed.
#' @param ... Additional parameters to be passed down.
#' @return A print to the console.
#' 
#' @examples
#' set.seed(123)
#' x <- sample(c(1,2,3), size = 100, replace = TRUE)
#' y <- sample(c(1,2,3), size = 100, replace = TRUE)
#' fit <- polycor(x,y) 
#' 
#' print(fit)
#' fit # equivalent
#' 
#' @export
print.robpolycor <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  stopifnot(inherits(x = x, what = "polycor"))
  Kx <- x$inputs$Kx ; Ky <- x$inputs$Ky
  thetaobj <- extractfromtheta(theta = x$thetahat, Kx = Kx, Ky = Ky, includeInf = FALSE)
  stderrobj <- extractfromtheta(theta = x$stderr, Kx = Kx, Ky = Ky, includeInf = FALSE)
  
  ## rho
  cat("\nPolychoric Correlation\n")
  rho <- signif(cbind(Estimate = thetaobj$rho, Std.Err. = stderrobj$rho), 
                digits = digits)
  print(rho)
  
  ## thresX
  thresX <- signif(cbind(Estimate = thetaobj$thresX, Std.Err. = stderrobj$thresX), 
                   digits = digits)
  cat("\nX-thresholds\n") 
  print(thresX)
  
  ## thresY
  thresY <- signif(cbind(Estimate = thetaobj$thresY, Std.Err. = stderrobj$thresY), 
                   digits = digits)
  cat("\nY-thresholds\n") 
  print(thresY)
}
