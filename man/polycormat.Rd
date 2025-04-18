% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/polycormat.R
\name{polycormat}
\alias{polycormat}
\title{Robust estimation of polychoric correlation matrix}
\usage{
polycormat(
  data,
  c = 0.6,
  parallel = FALSE,
  num_cores = 1L,
  return_polycor = TRUE,
  variance = TRUE,
  constrained = TRUE,
  method = ifelse(constrained, "Nelder-Mead", "L-BFGS-B"),
  maxcor = 0.999,
  tol_thresholds = 0.01
)
}
\arguments{
\item{data}{Data matrix or \code{\link[base]{data.frame}} of integer-valued responses, individual respondents are in rows and responses to the items in the columns.}

\item{c}{tuning constant that governs robustness; takes values in \code{[0, Inf]}. Defaults to 0.6.}

\item{parallel}{Logical. Shall parallelization be used? Default is \code{FALSE}.}

\item{num_cores}{Number of cores to be used, only relevant if \code{parallel = TRUE}. Defaults to the number of system cores.}

\item{return_polycor}{Logical. Shall the individual \code{"\link{polycor}"} objects for each item pair estimate be returned? Deafult is \code{TRUE}.}

\item{variance}{Shall an estimated asymptotic covariance matrix be returned? Default is \code{TRUE}.}

\item{constrained}{Shall strict monotonicity of thresholds be explicitly enforced by linear constraints?}

\item{method}{Numerical optimization method.}

\item{maxcor}{Maximum absolute correlation (to insure numerical stability).}

\item{tol_thresholds}{Minimum distance between consecutive thresholds (to enforce strict monotonicity); only relevant if \code{constrained = TRUE}.}
}
\value{
If \code{return_polycor = TRUE}, returns a list with a polychoric correlation matrix and list of \code{"\link{polycor}"} objects. If \code{return_polycor = FALSE}, then only a correlation matrix is returned.
}
\description{
A useful wrapper of \code{\link{polycor}} to robustly estimate a polychoric correlation matrix by calculating all unique pairwise polychoric correlation coefficients.
}
\examples{
## example data
set.seed(123)
data <- matrix(sample(c(1,2,3), size = 3*100, replace = TRUE), nrow = 100)
polycormat(data)     # robust 
polycormat_mle(data) # non-robust MLE

}
