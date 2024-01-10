#' example
#' 
#' @param x number
#' @param y number
#' @export
add_r <- function(x, y)
{
  return(add(x,y)) # a C++ function
}

#' hi
#' 
#' @param x some number
#' @export
ident <- function(x) x