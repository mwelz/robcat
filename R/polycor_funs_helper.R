# get thresholds for a given number of Likert points (=K)
get_thresholds <- function(K)
{
  c(-Inf, 1:(K-1) - stats::median(1:(K-1)), Inf)
} # FUN


# Initialize thresholds for K answer categories
init_thresholds <- function(K)
{
  tmp <- get_thresholds(K)
  out <- tmp[-c(1, K+1)] # eliminate the infinite fixed values
  out
} # FUN


#' Neutral initialization of starting values
#' 
#' Initializes starting values for numerical optimization in a neutral way. The optimization problem itself is convex, so the initialization should not matter much.
#' 
#' @param x Vector of integer-valued responses to first rating variable, or contingency table (a \code{table} object).
#' @param y Vector of integer-valued responses to second rating variable; only required if \code{x} is not a contingency table. 
#' @export
initialize_param <- function(x, y)
{
  if(is.table(x))
  {
    Kx <- nrow(x)
    Ky <- ncol(x)
  } else
  {
    Kx <- max(x, na.rm = TRUE)
    Ky <- max(y, na.rm = TRUE)
  }
  c(0, init_thresholds(Kx), init_thresholds(Ky))
}


# get discretized versions of a latent variable, governed by the thresholds
# TODO: Consider deleting this because we only need it in simulations
get_discretization <- function(thresholds, latent)
{
  num_likert <- length(thresholds) - 1
  y <- rep(NA_integer_, length(latent))
  
  for(j in 2:length(thresholds))
  {
    response <- j - 1L
    y[thresholds[j-1] <= latent & latent < thresholds[j]] <- response
  }
  y
}


# get the matrix and vector that enforce the linear constraints in stats::constrOptim()
constrOptim_constraints <- function(Kx, Ky, 
                                    maxcor = 0.999,
                                    tol_thresholds = 0.001)
{
  ## number of constraints needed for each parameter
  # if an item is dichotomous -> no constraints required for its thresholds
  num_constr_rho <- 2L         # -1 < rho & rho < 1
  num_constr_thresX <- Kx - 2L # monotoncity of X-thresholds
  num_constr_thresY <- Ky - 2L # monotoncity of Y-thresholds
  num_constr <- num_constr_rho + num_constr_thresX + num_constr_thresY
  num_param <- Kx + Ky - 1L
  
  
  ## Matrix U and vector C (ui and ci, respectively, in ?stats::constrOptim)
  # U is a (Kx+Ky-2) x (Kx+Ky-1) block matrix of zeros as base
  # from top left to bottom right, we have the following three submatrices:
  # for \rho: [1:2, 1]
  # for X-submatrix: [(2+1):(Kx-2+2), (1+1):(Kx-1+1)] = [3:Kx,2:Kx] is (Kx-2) x (Kx-1)
  # for Y-submatrix: [(2+(Kx-2)+1):(Kx+Ky-2), (1+(Kx-1)+1):(Kx+Ky-1)] = [(Kx+1):Km2, (Kx+1):Km1] is (Ky-2) x (Ky-1)
  U <- matrix(0L, nrow = num_constr, ncol = num_param)
  C <- rep(NA_real_, num_constr)
  
  ## rho-constraints
  U[seq_len(num_constr_rho), 1L] <- c(1L, -1L)  # enforce \rho \in [-1,1]
  C[seq_len(num_constr_rho)] <- maxcor * (-1.0) # tolerance for \rho \in [-1,1]
  
  
  ## X-threshold constraints: enforce monotonicity of adjacent thresholds
  # only applicable if item isn't dichotomous
  if(num_constr_thresX > 0L)
  {
    U[1:num_constr_thresX + num_constr_rho, 2:Kx] <- constrOptim_submatrix(Kx)
    C[1:num_constr_thresX + num_constr_rho] <- tol_thresholds # tolerance
  } 
  
  ## Y-threshold constraints: enforce monotonicity of adjacent thresholds
  # only applicable if item isn't dichotomous
  if(num_constr_thresY > 0L)
  {
    U[1:num_constr_thresY + num_constr_thresX + num_constr_rho,
      (Kx+1L):num_param] <- constrOptim_submatrix(Ky)
    C[1:num_constr_thresY + num_constr_thresX + num_constr_rho] <- tol_thresholds # tolerance
  }
  
  return(list(ui = U, ci = C))
  
} # FUN


# submatrix that imposes the strictly increasing monotonicity of the thresholds
# downward stepwise pattern of c(-1,1)
constrOptim_submatrix <- function(K)
{
  
  # fails if K = 2 (but no constraints needed there anyway)
  stopifnot(K > 2L)
  
  # initialize
  Km1 <- K - 1L
  Km2 <- K - 2L
  A <- matrix(NA_integer_, nrow = Km2, ncol = Km1)
  
  for(i in seq_len(Km2))
  {
    tmp <- rep(0L, Km1)
    tmp[i:(i+1L)] <- c(-1L, 1L)
    A[i,] <- tmp
  } # FOR
  return(A)
}


# give names to the dimensions of theta vector
theta_names <- function(Kx, Ky)
{
  c("rho", paste0("a", 1:(Kx-1)), paste0("b", 1:(Ky-1)))
}


## extract parameters from theta parameter
extractfromtheta <- function(theta, Kx, Ky, includeInf = TRUE)
{
  idxX <- seq(from = 2, to = Kx)
  idxY <- seq(from = Kx+1L, to = Kx+Ky-1L)
  if(includeInf) {
    thresX <- c(-Inf, theta[idxX], Inf) 
    thresY <- c(-Inf, theta[idxY], Inf) 
  } else{
    thresX <- theta[idxX]
    thresY <- theta[idxY]
  }
  return(list(rho = theta[1L], thresX = thresX, thresY = thresY))
}


## cast matrix to table
mat2tab <- function(mat)
{
  Kx <- nrow(mat) ; Ky <- ncol(mat)
  tab <- as.table(mat)
  names(attr(tab, "dimnames")) = c("x", "y")
  attr(tab, "dimnames")$x <- as.character(seq_len(Kx))
  attr(tab, "dimnames")$y <- as.character(seq_len(Ky))
  return(tab)
}


## create contingency table based on a vector that was obtained by fhat() or model_probablities()
vec2tab <- function(probs_vec, Kx, Ky, byrow = TRUE)
{
  if(is.null(probs_vec))
  {
    out <- NULL
  } else{
    ## create table/matrix in the same order than they were created (that is, rows first)
    out <- mat2tab(matrix(probs_vec, nrow = Kx, ncol = Ky,
                          byrow = byrow))
  }
  return(out)
}
