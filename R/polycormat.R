## get unique distinct pairs in a pxp correlation matrix
# returns list of such pairs
unique_pairs <- function(p)
{
  ## initialize list (number of unique distinct item pairs)
  num_pairs <- p * (p - 1) / 2
  x <- rep(list(NA), num_pairs)
  nam <- rep(NA_character_, num_pairs)
  k <- 1L
  
  for(i in seq_len(p))
  {
    for(j in seq_len(p))
    {
      if(j < i){
        x[[k]] <- c(i, j)
        nam[k] <- paste0("(", i, ",", j, ")")
        k <- k + 1L
      } else{
        next
      }
    }
  }
  names(x) <- nam
  return(x)
}


## unique_pairs_ls is returned by 'unique_pairs'
## polycor_fn is a function object that is either polycor_mle() or polycor() and only has x,y as arguments 
## all other arguments in the call to polycor() or polycor_mle() is fixed because it's constant for all item pairs
pairwise_polycor <- function(polycor_fn, unique_pairs_ls, data, num_cores, parallel)
{
  if(parallel)
  {
    out <- 
      parallel::mclapply(seq_along(unique_pairs_ls), mc.cores = num_cores, FUN = function(i){
        
        ## get responses to the item pair
        pair <- unique_pairs_ls[[i]]
        # system(sprintf('echo "%s\n"', paste0("start ", pair[1], " ", pair[2])))
        x <- data[, pair[1L]]
        y <- data[, pair[2L]]
        
        ## run pairwise polycor
        pp <- polycor_fn(x = x, y = y)
        # system(sprintf('echo "%s\n"', paste0("finished ", pair[1], " ", pair[2])))
        return(pp)
      })
  } else
  {
    out <- 
      lapply(seq_along(unique_pairs_ls), FUN = function(i){
        
        ## get responses to the item pair
        pair <- unique_pairs_ls[[i]]
        x <- data[, pair[1L]]
        y <- data[, pair[2L]]
        
        ## run pairwise polycor
        pp <- polycor_fn(x = x, y = y)
        return(pp)
      })
  } # IF
 
  names(out) <- names(unique_pairs_ls)
  return(out)
}


get_correlations <- function(data, c, variance, constrained, method, maxcor, tol_thresholds, num_cores, parallel, mle)
{
  p <- ncol(data)
  unique_pairs_ls <- unique_pairs(p)
  
  ## prepare function
  if(mle)
  {
    polycor_fn <- function(x, y)
    {
      polycor_mle(x = x, y = y, variance = variance, constrained = constrained,
                  method = method, maxcor = maxcor, tol_thresholds = tol_thresholds,
                  init = initialize_param(x, y))
    }
  } else
  {
    polycor_fn <- function(x, y)
    {
      polycor(x = x, y = y, c = c, variance = variance, constrained = constrained,
              method = method, maxcor = maxcor, tol_thresholds = tol_thresholds,
              init = initialize_param(x, y))
    }
  }
  
  ## get the pairwise correlations in list form
  out <- 
    pairwise_polycor(polycor_fn = polycor_fn, unique_pairs_ls = unique_pairs_ls, 
                     data = data, num_cores = num_cores, parallel = parallel)
  return(out)
  
} # FUN


## 'correlations' is returned by get_correlations()
## transform it to an actual correlation matrix
correlations2cormat <- function(correlations, p)
{
  item_pairs <- names(correlations)
  cormat <- matrix(NA_real_, p, p)
  
  for(i in seq_len(p))
  {
    for(j in seq_len(p))
    {
      if(i == j)
      {
        cormat[i,i] <- 1.0
      } else if(j < i)
      {
        pair <- paste0("(", i, ",", j, ")")
        k <- which(item_pairs == pair)
        cormat[i,j] <- cormat[j,i] <- correlations[[k]][["thetahat"]]["rho"]
      } else next
    }
  }
  return(cormat)
}


## main fucntionality of polychoric correlation matrix
get_cormat <- function(data, c, variance, constrained, method, maxcor, tol_thresholds, num_cores, parallel, mle, return_polycor)
{
  correlations <- get_correlations(data = data, c = c, variance = variance, constrained = constrained, 
                                   method = method, maxcor = maxcor, tol_thresholds = tol_thresholds,
                                   num_cores = num_cores, parallel = parallel, mle = mle)
  mat <- correlations2cormat(correlations, p = ncol(data))
  colnames(mat) <- rownames(mat) <- colnames(data)
  
  if(return_polycor)
  {
    out <- list(cormat = mat, polycor_objects = correlations)
  } else{
    out <- mat
  }
  return(out)
}


#' Robust estimation of polychoric correlation matrix
#' 
#' @param data Data matrix of integer-valued responses, individual respondents are in rows
#' @param c tuning constant that governs robustness; takes values in [0,Inf]. Defaults to 0.6.
#' @param parallel Logical. Shall parallelization be used?
#' @param num_cores Number of cores to be used, only relevant if \code{parallel = TRUE}
#' @param return_polycor Logical. Shall the individual \code{polycor} objects for each item pair estimate be returned?
#' @param variance shall an estimated asymptotic covariance matrix be returned? Default is \code{TRUE}
#' @param method numerical optimization method
#' @param constrained shall strict monotonicity of thresholds be explicitly enforced by linear constraints? 
#' @param maxcor maximum absolute correlation (to insure numerical stability)
#' @param tol_thresholds minimum distance between consecutive thresholds (to enforce strict monotonicity); only relevant if \code{constrained = TRUE}
#' 
#' @export
polycormat <- function(data, c = 0.6, 
                       parallel = FALSE, 
                       num_cores = 1L,
                       return_polycor = TRUE,
                       variance = TRUE,
                       constrained = TRUE,
                       method = ifelse(constrained, "Nelder-Mead", "L-BFGS-B"),
                       maxcor = 0.999,
                       tol_thresholds = 0.01)
{
  get_cormat(data = data, c = c, variance = variance, constrained = constrained, 
             method = method, maxcor = maxcor, tol_thresholds = tol_thresholds,
             num_cores = num_cores, parallel = parallel, mle = FALSE, return_polycor = return_polycor)
}
  

#' Maximum likelihood estimation of polychoric correlation matrix
#' 
#' @param data Data matrix of integer-valued responses, individual respondents are in rows
#' @param parallel Logical. Shall parallelization be used?
#' @param num_cores Number of cores to be used, only relevant if \code{parallel = TRUE}
#' @param return_polycor Logical. Shall the individual \code{polycor} objects for each item pair estimate be returned?
#' @param variance shall an estimated asymptotic covariance matrix be returned? Default is \code{TRUE}
#' @param method numerical optimization method
#' @param constrained shall strict monotonicity of thresholds be explicitly enforced by linear constraints? 
#' @param maxcor maximum absolute correlation (to insure numerical stability)
#' @param tol_thresholds minimum distance between consecutive thresholds (to enforce strict monotonicity); only relevant if \code{constrained = TRUE}
#'
#' @export
polycormat_mle <- function(data,
                           parallel = FALSE, 
                           num_cores = 1L,
                           return_polycor = TRUE,
                           variance = TRUE,
                           constrained = TRUE,
                           method = ifelse(constrained, "Nelder-Mead", "L-BFGS-B"),
                           maxcor = 0.999,
                           tol_thresholds = 0.01)
{
  get_cormat(data = data, c = Inf, variance = variance, constrained = constrained, 
             method = method, maxcor = maxcor, tol_thresholds = tol_thresholds,
             num_cores = num_cores, parallel = parallel, mle = TRUE, return_polycor = return_polycor)
}

## TODO:  possible eigenvalue correction and windows; palallel in description & namespace, export funs