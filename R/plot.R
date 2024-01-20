#' Plot method for class "polycor"
#' 
#' @param x object of class \code{polycor}
#' @param cutoff cutoff beyond which a Pearson residual is classified as outlying
#' @param ... additional parameters to be passed down
#' @import ggplot2
#' @export
plot.robpolycor <- function(x, cutoff = stats::qnorm(0.001, lower.tail = FALSE), ...)
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
  x <- y <- Residual <- Frequency <- Outlier <- NULL
  
  ## is a cell outlying?
  outlier_bool <- df$Residual > cutoff
  
  # prepare data frame
  df$x <- factor(df$x, levels = seq_len(Kx))
  df$y <- factor(df$y, levels = rev(seq_len(Ky)))
  df$Outlier <- paste0("PR > ", round(cutoff, 2))
  df$Outlier[df$Residual <= cutoff] <- NA 
  titlecolor <- paste0("Pearson\nResidual", ifelse(any(outlier_bool), "\n(PR)", ""))
  
  ## the idea to use separate scales for outliers is from here:
  # https://stackoverflow.com/a/9812648
  p <- 
    ggplot() +
    scale_color_gradient("Residual", low = "lightblue", high = "darkblue")+
    geom_point(data = subset(df, outlier_bool), 
               mapping = aes(x = y, y = x, 
                             fill = stringr::str_wrap(Outlier, 3),
                             size = Frequency), 
               color = "red") +
    geom_point(data = subset(df, !outlier_bool), 
               mapping = aes(x = y, y = x, 
                             color = Residual, 
                             size = Frequency)) +
    xlab("Y") + ylab("X") +
    scale_y_discrete(limits = rev) + # flip y-axis
    scale_x_discrete(limits = rev) + # flip x-axis
    guides(size = guide_legend(order = 1, title = "Relative\nFrequency"),
           color  = guide_colorbar(order = 2, title = titlecolor), 
           fill = guide_legend(order = 3, title = "Poor\nFit")) 
  return(p)
}
