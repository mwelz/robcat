## extract parameters from theta parameter
polyserial_extractfromtheta <- function(theta, num_y)
{
  d        <- num_y + 2L
  idxX     <- seq(from = 2L, to = 3L)
  idxY     <- seq(from = 4L, to = d)
  momentsX <- theta[idxX]
  thres    <- theta[idxY]
  return(list(rho = theta[1L], momentsX = momentsX, thres = thres))
}
