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


### extract functions relevant for computing the A-matrix (theoretical part in J)
## get the function that takes derivative of p_x wrt to i-th parameter
get_px_d1_i <- function(i)
{
  ## marginal only depends on mu and sigma2, everything else is 0
  if(i == 2L)
  {
    f <- get("px_d1_mu")
  } else if(i == 3L)
  {
    f <- get("px_d1_sigma2")
  } else
  {
    f <- function(x, mu, sigma2){return(0.0)}
  }
  return(f)
}


## get the function that takes derivative of p_y|x wrt to i-th parameter
get_py_x_d1_i <- function(i)
{
  if(i == 1L)
  {
    f <- get("py_x_d1_rho")
  } else if(i == 2L)
  {
    f <- get("py_x_d1_mu")
  } else if(i == 3L)
  {
    f <- get("py_x_d1_sigma2")
  } else
  {
    f <- function(y, x, rho, mu, sigma2, thres, num_y)
    {
      py_x_d1_tauk(y = y, x = x, k = i - 3L, rho = rho, mu = mu,
                   sigma2 = sigma2, thres = thres, num_y = num_y)
    }
  }
  return(f)
}

## get the function that takes 2nd derivative of p_x wrt to (i,j)th parameter
# j >= i
get_px_d2_ij <- function(i, j)
{
  if(i == j)
  {
    if(i == 2L)
    {
      f <- get("px_d2_mu2")
    } else if(i == 3L)
    {
      f <- get("px_d2_sigma22")
    } else
    {
      ## all remaining 2nd derivatives are 0
      f <- function(x, mu, sigma2){return(0.0)}
    }
  } else
  {
    ## j > i
    if(i == 2L && j == 3L)
    {
      f <- get("px_d2_musigma2")
    } else
    {
      ## all remaining cross-derivatives are 0
      f <- function(x, mu, sigma2){return(0.0)}
    }
  }
  return(f)
}


## get the function that takes 2nd derivative of p_y|x wrt to (i,j)th parameter
# j >= i
get_py_x_d2_ij <- function(i, j)
{
  if(i == j)
  {
    ## CASE 1: 2nd order derivative
    if(i == 1L)
    {
      f <- get("py_x_d2_rho2")
    } else if(i == 2L)
    {
      f <- get("py_x_d2_mu2")
    } else if(i == 3L)
    {
      f <- get("py_x_d2_sigma22")
    } else
    {
      f <- function(y, x, rho, mu, sigma2, thres, num_y)
      {
        py_x_d2_tauk2(y = y, x = x, k = i - 3L, rho = rho, mu = mu,
                      sigma2 = sigma2, thres = thres, num_y = num_y)
      }
    }
  } else
  {
    ## CASE 2: cross-derivative: j > i
    if(i == 1L && j == 2L)
    {
      f <- get("py_x_d2_rhomu")
    } else if(i == 1L && j == 3L)
    {
      f <- get("py_x_d2_rhosigma2")
    } else if(i == 2L && j == 3L)
    {
      f <- get("py_x_d2_musigma2")
    } else if(i == 1L && j > 3L)
    {
      f <- function(y, x, rho, mu, sigma2, thres, num_y)
      {
        py_x_d2_rhotauk(y = y, x = x, k = j - 3L, rho = rho, mu = mu,
                        sigma2 = sigma2, thres = thres, num_y = num_y)
      }
    } else if(i == 2L && j > 3L)
    {
      f <- function(y, x, rho, mu, sigma2, thres, num_y)
      {
        py_x_d2_mutauk(y = y, x = x, k = j - 3L, rho = rho, mu = mu,
                       sigma2 = sigma2, thres = thres, num_y = num_y)
      }
    } else if(i == 3L && j > 3L)
    {
      f <- function(y, x, rho, mu, sigma2, thres, num_y)
      {
        py_x_d2_sigma2tauk(y = y, x = x, k = j - 3L, rho = rho, mu = mu,
                           sigma2 = sigma2, thres = thres, num_y = num_y)
      }
    } else
    {
      ## cross-derivative sof distinct thresholds: structural 0
      f <- function(y, x, rho, mu, sigma2, thres, num_y){return(0.0)}
    }
  }
  return(f)
}

