## check for NAs
check_NA <- function(x)
{
  if(any(is.na(x))) stop("Missing values detected in supplied data")
}

is_integer <- function(x)
{
  if(!is.vector(x)) stop("Responses must be supplied as vector or contingency table")
  if(!is.integer(x)) stop("Responses must be integer-valued for now")
}

## check if support of x is in {1,2,...,K}
check_support <- function(x, K)
{
  supp <- seq_len(K)
  all(unique(x) %in% supp)
}

## input check when x is a table
# y is something else that will be disregarded
input_table <- function(x, y)
{
  check_NA(x)
  if(!is.integer(x)) stop("Frequency counts in a contingency table must be integer-valued")
  if(!is.null(y)) message("In addition to a contingency table in x, data for y were also passed. The y-input is disregarded and only data from the contingency table will be used")
  N <- sum(x)
  f <- as.integer(t(x)) / as.numeric(N)
  if(!isTRUE(all.equal(target = 1, current = sum(f)))) stop("Relative frequencies do not sum up to 1")
  Kx <- nrow(x)
  Ky <- ncol(x)
  return(list(N = N, f = f, Kx = Kx, Ky = Ky))
}

## input check when x is a vector (enforces that y is also a table)
input_vector <- function(x, y)
{
  if(is.null(y)) stop("Since x is not a contingency table, you must also pass responses to the second item in y")
  check_NA(x) ; check_NA(y)
  is_integer(x); is_integer(y)
  N <- length(x)
  if(N != length(y)) stop("x and y must be of same length")
  Kx <- max(x)
  Ky <- max(y)
  check_support(x, Kx)
  check_support(y, Ky)
  f <- fhat(x = x, y = y, Kx = Kx, Ky = Ky)
  return(list(N = N, f = f, Kx = Kx, Ky = Ky))
}
