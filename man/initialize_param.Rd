% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/polycor_funs_helper.R
\name{initialize_param}
\alias{initialize_param}
\title{Neutral initialization of starting values}
\usage{
initialize_param(x, y)
}
\arguments{
\item{x}{Vector of integer-valued responses to first rating variable, or contingency table (a \code{table} object).}

\item{y}{Vector of integer-valued responses to second rating variable; only required if \code{x} is not a contingency table.}
}
\value{
A vector of initial values for the polychoric correlation coefficient, the X-threshold parameters, and the Y-threshold parameters
}
\description{
Initializes starting values for numerical optimization in a neutral way. The optimization problem itself is convex, so the initialization should not matter much.
}
\examples{
## example data
set.seed(123)
x <- sample(c(1,2,3), size = 100, replace = TRUE)
y <- sample(c(1,2,3), size = 100, replace = TRUE)
initialize_param(x, y)

}
