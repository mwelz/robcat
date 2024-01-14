#' Print method for class "polycor"
#' 
#' @param x object of class \code{polycor}
#' @param digits number of digits to be printed
#' @param ... additional parameters to be passed down
#' @export
print.polycor <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  stopifnot(inherits(x = x, what = "polycor"))
  Kx <- x$inputs$Kx ; Ky <- x$inputs$Ky
  thetaobj <- extractfromtheta(theta = x$thetahat, Kx = Kx, Ky = Ky, includeInf = FALSE)
  stderrobj <- extractfromtheta(theta = x$stderr, Kx = Kx, Ky = Ky, includeInf = FALSE)
  
  ## rho
  cat("\nPolychoric Correlation\n")
  rho <- signif(cbind(Coefficient = thetaobj$rho, Std.Err. = stderrobj$rho), 
                digits = digits)
  print(rho)
  
  ## thresX
  thresX <- signif(cbind(Threshold = thetaobj$thresX, Std.Err. = stderrobj$thresX), 
                   digits = digits)
  cat("\nX-thresholds\n") 
  print(thresX)
  
  ## thresY
  thresY <- signif(cbind(Threshold = thetaobj$thresY, Std.Err. = stderrobj$thresY), 
                   digits = digits)
  cat("\nY-thresholds\n") 
  print(thresY)
}
