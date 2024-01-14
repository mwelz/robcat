#' Plot method for class "polycor"
#' 
#' @param x object of class \code{polycor}
#' @param ... additional parameters to be passed down
#' @import ggplot2
#' @export
plot.polycor <- function(x, ...)
{
  stopifnot(inherits(x = x, what = "polycor"))
  resid <- x[["residuals"]]
  freq <- x[["f"]]
  Kx <- x$inputs$Kx; Ky <- x$inputs$Ky
  
  ## TODO: find better way to convert to long format
  df <- as.data.frame(matrix(NA, Kx * Ky, 4L))
  colnames(df) <- c("x", "y", "Frequency", "Residual")
  k <- 1L
  for(x in seq_len(Kx))
  {
    for(y in seq_len(Ky))
    {
      df[k,] <- c(x, y, freq[x,y], resid[x,y])
      k <- k + 1L
    }
  }
  
  ## initialize variables to avoid R CMD issue
  x <- y <- Residual <- Frequency <- NULL
  
  ## plot
  p <- 
    ggplot(df, aes(x = y, y = x, color = Residual, size = Frequency)) +
    geom_point() +
    scale_colour_gradient2(
      low = "orange",
      mid = "steelblue",
      high = "red",
      midpoint = 1) +
    theme_bw()
  return(p)
}
