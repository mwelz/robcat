#' Plot method for class "polycor"
#' 
#' @param x object of class \code{polycor}
#' @param cutoff cutoff beyond which a Pearson residual is classified as outlying
#' @param outliertest shall the outlyingness test results be used in plot?
#' @param twosided shall outlyingness test be two sided?
#' @param adjust shall outlyingness tests be adjusted for multiple comparisons?
#' @param sig_level significance level which is used in plotting tets results 
#' @param ... additional parameters to be passed down
#' @import ggplot2
#' @export
plot.robpolycor <- function(x, 
                            cutoff = stats::qnorm(0.001, lower.tail = FALSE),
                            outliertest = TRUE,
                            twosided = FALSE,
                            adjust = "fdr",
                            sig_level = 0.001,
                            ...)
{
  stopifnot(inherits(x = x, what = "polycor"))
  resid <- x[["residuals"]]
  freq <- x[["f"]]
  Kx <- x$inputs$Kx; Ky <- x$inputs$Ky
  
  if(!outliertest)
  {
    plot_truncated(resid = resid, freq = freq, Kx = Kx, Ky = Ky, cutoff = cutoff)
  } else{
    pvals <- celltest(x, twosided = twosided, adjust = adjust)[["pval_adjusted"]]
    plot_test(resid = resid, freq = freq, pval = pvals, Kx = Kx, Ky = Ky, sig_level = sig_level)
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
  outlier_bool <- df$Residual > cutoff
  
  # prepare data frame
  Xseq <- rev(seq_len(Kx))
  Yseq <- seq_len(Ky)
  df$x <- factor(df$x, levels = Xseq)
  df$y <- factor(df$y, levels = Yseq)
  df$Outlier <- paste0("PR > ", round(cutoff, 2))
  df$Outlier[!outlier_bool] <- NA 
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
    scale_y_discrete(limits = as.character(Xseq)) + # flip y-axis
    scale_x_discrete(limits = as.character(Yseq)) +
    guides(size = guide_legend(order = 1, title = "Relative\nFrequency"),
           color  = guide_colorbar(order = 2, title = titlecolor), 
           fill = guide_legend(order = 3, title = "Poor\nFit")) 
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