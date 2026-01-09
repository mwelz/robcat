## hessian matrix of p_{Y|X}; 2nd-order-derivatives

### 1.1 mu ----
py_x_d2_mu2 <- function(y, x, rho, mu, sigma2, thres, num_y)
{
  scaling <- rho / sqrt(sigma2 * (1 - rho^2))
  
  if(y == 1L)
  {
    deriv1 <- tstar_d1_mu(x = x, rho = rho, mu = mu, sigma2 = sigma2)
    phi1 <- phi_prime(get_taustar_y(x = x, rho = rho, mu = mu, sigma2 = sigma2, 
                                    tau_y = thres[y]))
    sum1 <- deriv1 * phi1
    sum2 <- 0.0
    
  } else if(y == num_y)
  {
    sum1 <- 0.0
    deriv2 <- tstar_d1_mu(x = x, rho = rho, mu = mu, sigma2 = sigma2)
    phi2 <- phi_prime(get_taustar_y(x = x, rho = rho, mu = mu, sigma2 = sigma2, 
                              tau_y = thres[y-1L]))
    sum2 <- deriv2 * phi2
  } else
  {
    deriv1 <- tstar_d1_mu(x = x, rho = rho, mu = mu, sigma2 = sigma2)
    phi1 <- phi_prime(get_taustar_y(x = x, rho = rho, mu = mu, sigma2 = sigma2, 
                                    tau_y = thres[y]))
    sum1 <- deriv1 * phi1
    
    deriv2 <- tstar_d1_mu(x = x, rho = rho, mu = mu, sigma2 = sigma2)
    phi2 <- phi_prime(get_taustar_y(x = x, rho = rho, mu = mu, sigma2 = sigma2, 
                                    tau_y = thres[y-1L]))
    sum2 <- deriv2 * phi2
  }
  
  diff <- sum1 - sum2
  return(scaling * diff)
}


### 1.2 sigma2 ----
py_x_d2_sigma22 <- function(y, x, rho, mu, sigma2, thres, num_y)
{
  if(y == 1L)
  {
    ## terms for 1st difference
    diff1_part1 <- phi(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y]))
    diff1_part2 <- 0.0
    
    ## terms for 2nd difference
    diff2_phi1 <- phi_prime(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y]))
    diff2_deriv1 <- tstar_d1_sigma2(x = x, rho = rho, mu = mu, sigma2 = sigma2)
    diff2_part1 <- diff2_phi1 * diff2_deriv1
    diff2_part2 <- 0.0
    
  } else if(y == num_y)
  {
    ## terms for 1st difference
    diff1_part1 <- 0.0
    diff1_part2 <- phi(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y-1L]))
    
    ## terms for 2nd difference
    diff2_part1 <- 0.0
    diff2_phi2 <- phi_prime(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y-1L]))
    diff2_deriv2 <- tstar_d1_sigma2(x = x, rho = rho, mu = mu, sigma2 = sigma2)
    diff2_part2 <- diff2_phi2 * diff2_deriv2
    
  } else
  {
    ## terms for 1st difference
    diff1_part1 <- phi(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y]))
    diff1_part2 <- phi(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y-1L]))
    
    ## terms for 2nd difference
    diff2_phi1 <- phi_prime(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y]))
    diff2_deriv1 <- tstar_d1_sigma2(x = x, rho = rho, mu = mu, sigma2 = sigma2)
    diff2_part1 <- diff2_phi1 * diff2_deriv1
    
    diff2_phi2 <- phi_prime(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y-1L]))
    diff2_deriv2 <- tstar_d1_sigma2(x = x, rho = rho, mu = mu, sigma2 = sigma2)
    diff2_part2 <- diff2_phi2 * diff2_deriv2
  }
  
  ## assemble the final term
  scaling <- rho * (x - mu) / (2 * sigma2^1.5 * sqrt(1 - rho^2))
  diff1 <- (-1.5) * (diff1_part1 - diff1_part2) / sigma2
  diff2 <- diff2_part1 - diff2_part2
  return(scaling * (diff1 + diff2))
}


### 1.3 tau_k ----
# defined for k \in [num_y-1]
py_x_d2_tauk2 <- function(y, x, k, rho, mu, sigma2, thres, num_y)
{
  stopifnot(k < num_y)
  ym1     <- y - 1L
  scaling <- 1.0 / sqrt(1.0 - rho^2)
  
  if(k == y)
  {
    a <- phi_prime(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y]))
    b <- tstar_d1_tauk(
      x = x, y = y, k = y, rho = rho, mu = mu, sigma2 = sigma2)
    out <- a * b * scaling
  } else if(k == ym1)
  {
    a <- phi_prime(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[ym1]))
    b <- tstar_d1_tauk(
      x = x, y = y, k = y, rho = rho, mu = mu, sigma2 = sigma2)
    out <- (-1) * a * b * scaling
  } else
  {
    out <- 0.0
  }
  return(out)
}


### 1.4 rho ----
py_x_d2_rho2 <- function(y, x, rho, mu, sigma2, thres, num_y)
{
  z <- (x - mu) / sqrt(sigma2)
  
  if(y == 1L)
  {
    # 1st parenthesis: part 1
    paren1_part1a <- phi(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y]))
    paren1_part1b <- rho * thres[y] - z
    paren1_part1 <- paren1_part1a * paren1_part1b
    
    # 1st parenthesis: part 2
    paren1_part2 <- 0.0
    
    # 2nd parenthesis: part 1
    paren2_part1a <- phi_prime(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y]))
    paren2_part1b <- rho * thres[y] - z
    paren2_part1c <- tstar_d1_rho(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y])
    paren2_part1d <- thres[y] * paren1_part1a
    paren2_part1 <- 
      paren2_part1a * paren2_part1b * paren2_part1c + paren2_part1d
    
    # 2nd parenthesis: part 2
    paren2_part2 <- 0.0
    
  } else if(y == num_y)
  {
    # 1st parenthesis: part1
    paren1_part1 <- 0.0
    
    # 1st parenthesis: part 2
    paren1_part2a <- phi(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y-1L]))
    paren1_part2b <- rho * thres[y-1L] - z
    paren1_part2 <- paren1_part2a * paren1_part2b
    
    # 2nd parenthesis: part1
    paren2_part1 <- 0.0
    
    # 2nd parenthesis: part 2
    paren2_part2a <- phi_prime(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y-1L]))
    paren2_part2b <- rho * thres[y-1L] - z
    paren2_part2c <- tstar_d1_rho(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y-1L])
    paren2_part2d <- thres[y-1L] * paren1_part2a
    paren2_part2 <- 
      paren2_part2a * paren2_part2b * paren2_part2c + paren2_part2d
    
  } else
  {
    # 1st parenthesis: part 1
    paren1_part1a <- phi(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y]))
    paren1_part1b <- rho * thres[y] - z
    paren1_part1 <- paren1_part1a * paren1_part1b
    
    # 1st parenthesis: part 2
    paren1_part2a <- phi(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y-1L]))
    paren1_part2b <- rho * thres[y-1L] - z
    paren1_part2 <- paren1_part2a * paren1_part2b
    
    # 2nd parenthesis: part 1
    paren2_part1a <- phi_prime(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y]))
    paren2_part1b <- rho * thres[y] - z
    paren2_part1c <- tstar_d1_rho(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y])
    paren2_part1d <- thres[y] * paren1_part1a
    paren2_part1 <- 
      paren2_part1a * paren2_part1b * paren2_part1c + paren2_part1d
    
    # 2nd parenthesis: part 2
    paren2_part2a <- phi_prime(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y-1L]))
    paren2_part2b <- rho * thres[y-1L] - z
    paren2_part2c <- tstar_d1_rho(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y-1L])
    paren2_part2d <- thres[y-1L] * paren1_part2a
    paren2_part2 <- 
      paren2_part2a * paren2_part2b * paren2_part2c + paren2_part2d
  }
  
  ## assemble the final term
  scaling1 <- 3 * rho / ((1 - rho^2)^2.5)
  diff1 <- scaling1 * (paren1_part1 - paren1_part2)
  
  scaling2 <- 1 / ((1 - rho^2)^1.5)
  diff2 <- scaling2 * (paren2_part1 - paren2_part2)
  
  return(diff1 + diff2)
}
