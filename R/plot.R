#' Plot method for classes \code{"robpolycor"} and \code{"polycor"}.
#' 
#' @param x Object of class \code{"robpolycor"} or \code{"polycor"}.
#' @param cutoff Cutoff beyond which the color scale for Pearson residuals is truncated.
#' @param ... Additional parameters to be passed down.
#' @import ggplot2
#' @export
plot.robpolycor <- function(x, 
                            cutoff = 3, #stats::qnorm(0.001, lower.tail = FALSE),
                            #outliertest = FALSE,
                            #twosided = FALSE,
                            #adjust = "fdr",
                            #sig_level = 0.001,
                            ...)
{
  stopifnot(inherits(x = x, what = "polycor"))
  resid <- x[["residuals"]]
  freq <- x[["f"]]
  Kx <- x$inputs$Kx; Ky <- x$inputs$Ky
  
  outliertest <- FALSE # not yet reliable (also update plot_test())!
  if(!outliertest)
  {
    plot_truncated(resid = resid, freq = freq, Kx = Kx, Ky = Ky, cutoff = cutoff)
  } else{
    stop("not yet ready!")
    #pvals <- celltest(x, twosided = twosided, adjust = adjust)[["pval_adjusted"]]
    #plot_test(resid = resid, freq = freq, pval = pvals, Kx = Kx, Ky = Ky, sig_level = sig_level)
  }
  
  
}

#' make plot with truncated pearson residuals
#' @import ggplot2
#' @noRd
plot_truncated <- function(resid, freq, Kx, Ky, cutoff)
{
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
  is_outlier <- df$Residual > cutoff
  
  # prepare data frame
  Xseq <- rev(seq_len(Kx))
  Yseq <- seq_len(Ky)
  df$x <- factor(df$x, levels = Xseq)
  df$y <- factor(df$y, levels = Yseq)
  df$Residual[is_outlier] <- cutoff
  #titlecolor <- paste0("Pearson\nResidual", ifelse(any(outlier_bool), "\n(PR)", ""))
  
  ## the idea to use separate scales for outliers is from here:
  # https://stackoverflow.com/a/9812648
  p <-
    ggplot() +
    theme_bw() +
    scale_color_gradient2("Residual", 
                          low = "blue", 
                          mid = "gray", 
                          high = "red", 
                          midpoint = 0, 
                          transform = "log1p", # log(x+1), to take care of negative PR values 
                          limits = c(0.5, 4)-1, 
                          breaks = c(0.5, 1, 2, 4)-1, 
                          labels = c("-0.5", "0.0", "1.0", "\u2265 3.0")) +
    geom_point(data = df, 
               mapping = aes(x = y, y = x, 
                             color = Residual, 
                             size = Frequency)) +
    scale_y_discrete(limits = as.character(Xseq)) + 
    scale_x_discrete(limits = as.character(Yseq)) +
    coord_fixed() +
    guides(size = guide_legend(order = 1, title = "Relative\nFrequency"),
           color  = guide_colorbar(order = 2, title = "Pearson\nResidual")) +
    xlab("Y") + ylab("X")
  return(p)
}


#' make plot with outlyingness test 
#' @import ggplot2
#' @noRd
plot_test <- function(resid, freq, pval, Kx, Ky, sig_level)
{
  ## TODO: find better way to convert to long format
  df <- as.data.frame(matrix(NA, Kx * Ky, 5L))
  colnames(df) <- c("x", "y", "Frequency", "Residual", "pval")
  k <- 1L
  for(x in seq_len(Kx))
  {
    for(y in seq_len(Ky))
    {
      df[k,] <- c(x, y, freq[x,y], resid[x,y], pval[x,y])
      k <- k + 1L
    }
  }
  
  ## initialize variables to avoid R CMD issue
  x <- y <- Residual <- Frequency <- Outlier <- NULL
  
  ## is a cell outlying?
  outlier_bool <- df$pval < sig_level
  
  # prepare data frame
  Xseq <- rev(seq_len(Kx))
  Yseq <- seq_len(Ky)
  df$x <- factor(df$x, levels = Xseq)
  df$y <- factor(df$y, levels = Yseq)
  df$Outlier <- paste0("p < ", sig_level)
  df$Outlier[!outlier_bool] <- NA 
  titlecolor <- "Pearson\nResidual"
  
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
    scale_y_discrete(limits = as.character(Xseq)) + # flip y-axis
    scale_x_discrete(limits = as.character(Yseq)) +
    guides(size = guide_legend(order = 1, title = "Relative\nFrequency"),
           color  = guide_colorbar(order = 2, title = titlecolor), 
           fill = guide_legend(order = 3, title = "Outlier")) 
  return(p)
}