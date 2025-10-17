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



#' Print method for classes \code{"robpolyserial"} and \code{"polyserial"}.
#' 
#' @param x Object of class \code{"robpolyserial"} or \code{"polyserial"}.
#' @param digits Number of digits to be printed.
#' @param ... Additional parameters to be passed down.
#' @return A print to the console.
#' 
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#' y <- sample(c(1,2), size = 100, replace = TRUE)
#' fit <- polyserial(x,y) 
#' 
#' print(fit)
#' fit # equivalent
#' 
#' @export
print.robpolyserial <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  thetaobj  <- polyserial_extractfromtheta(theta = x$thetahat, num_y = x$inputs$num_y)
  stderrobj <- polyserial_extractfromtheta(theta = x$stderr, num_y = x$inputs$num_y)
  
  ## rho
  cat("\nPolyserial Correlation\n")
  rho <- signif(cbind(Estimate = thetaobj$rho, Std.Err. = stderrobj$rho), 
                digits = digits)
  print(rho)
  
  ## moments of X
  thresX <- signif(cbind(Estimate = thetaobj$momentsX, Std.Err. = stderrobj$momentsX), 
                   digits = digits)
  cat("\nMoments of X\n") 
  print(thresX)
  
  ## thresholds of Y
  thresY <- signif(cbind(Estimate = thetaobj$thres, Std.Err. = stderrobj$thres), 
                   digits = digits)
  cat("\nThresholds of Y\n") 
  print(thresY)
}
