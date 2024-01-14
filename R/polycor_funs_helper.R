# get thresholds for a given number of Likert points (=K)
get_thresholds <- function(K)
{
  c(-Inf, 1:(K-1) - stats::median(1:(K-1)), Inf)
} # FUN


#' Initialize thresholds for K answer categories
#' @param K number of answer categories
#' @export
init_thresholds <- function(K)
{
  tmp <- get_thresholds(K)
  out <- tmp[-c(1, K+1)] # eliminiate the infinite fixed values
  out
} # FUN


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
                                    tol_rho = 0.001,
                                    tol_thresholds = 0.001)
{
  Km1 <- Kx + Ky - 1L
  Km2 <- Km1 - 1L
  
  ### the U matrix (ui in ?stats::constrOptim)
  # U is a (Kx+Ky-2) x (Kx+Ky-1) block matrix of zeros as base
  # from top left to bottom right, we have the following three submatrices:
  # for \rho: [1:2, 1]
  # for X-submatrix: [(2+1):(Kx-2+2), (1+1):(Kx-1+1)] = [3:Kx,2:Kx] is (Kx-2) x (Kx-1)
  # for Y-submatrix: [(2+(Kx-2)+1):(Kx+Ky-2), (1+(Kx-1)+1):(Kx+Ky-1)] = [(Kx+1):Km2, (Kx+1):Km1] is (Ky-2) x (Ky-1)
  U    <- matrix(0L, nrow = Km2, ncol = Km1)
  subx <- constrOptim_submatrix(Kx)
  suby <- constrOptim_submatrix(Ky)
  U[1:2,1]                  <- c(1L, -1L) # enforce \rho \in [-1,1]
  U[3:Kx, 2:Kx]             <- subx       # enforce monotonicity of X-thresholds
  U[(Kx+1):Km2, (Kx+1):Km1] <- suby       # enforce monotonicity of Y-thresholds
  
  ### the C vector (ci in ?stats::constrOptim)
  C        <- rep(NA_real_, Km2)
  C[1:2]   <- -(1.0 - tol_rho) # tolerance for \rho \in [-1,1]
  C[3:Km2] <- tol_thresholds # tolerance for monotonicity of thresholds
  
  return(list(ui = U, ci = C))
  
} # FUN


# submatrix that imposes the strictly increasing monotonicity of the thresholds
# downward stepwise pattern of c(-1,1)
constrOptim_submatrix <- function(K)
{
  stopifnot(K > 2L) # fails if K = 2 (but no constraints needed there anyway)
  
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
extractfromtheta <- function(theta, Kx, Ky)
{
  idxX <- seq(from = 2, to = Kx)
  idxY <- seq(from = Kx+1L, to = Kx+Ky-1L)
  thresX <- c(-Inf, theta[idxX], Inf) 
  thresY <- c(-Inf, theta[idxY], Inf) 
  return(list(rho = theta[1L], thresX = thresX, thresY = thresY))
}


## create contingency table based on a vector that was obtained by fhat() or model_probablities()
vec2tab <- function(probs_vec, Kx, Ky)
{
  if(is.null(probs_vec))
  {
    out <- NULL
  } else{
    ## create table/matrix in the same order than they were created (that is, rows first)
    tab <- as.table(
      matrix(probs_vec, nrow = Kx, ncol = Ky,
             byrow = TRUE))
    
    ## adjust naming
    names(attr(tab, "dimnames")) = c("x", "y")
    attr(tab, "dimnames")$x <- as.character(seq_len(Kx))
    attr(tab, "dimnames")$y <- as.character(seq_len(Ky))
    out <- tab
  }
  return(out)
}
