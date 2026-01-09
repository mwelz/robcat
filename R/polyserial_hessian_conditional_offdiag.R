## hessian matrix of p_{Y|X}; cross-derivatives

py_x_d2_rhomu <- function(y, x, rho, mu, sigma2, thres, num_y)
{
  if(y == 1L)
  {
    # 1st parenthesis
    paren1_part1 <- phi(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y]))
    paren1_part2 <- 0.0
    
    # 2nd parenthesis: part 1
    paren2_part1a <- phi_prime(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y]))
    paren2_part1b <- tstar_d1_rho(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y])
    paren2_part1 <- paren2_part1a * paren2_part1b
    
    # 2nd parenthesis: part 2
    paren2_part2 <- 0.0
    
  } else if(y == num_y)
  {
    # 1st parenthesis
    paren1_part1 <- 0.0
    paren1_part2 <- phi(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y-1L]))
    
    # 2nd parenthesis: part 1
    paren2_part1 <- 0.0
    
    # 2nd parenthesis: part 2
    paren2_part2a <- phi_prime(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y-1L]))
    paren2_part2b <- tstar_d1_rho(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y-1L])
    paren2_part2 <- paren2_part2a * paren2_part2b
    
  } else
  {
    # 1st parenthesis
    paren1_part1 <- phi(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y]))
    paren1_part2 <- phi(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y-1L]))
    
    # 2nd parenthesis: part 1
    paren2_part1a <- phi_prime(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y]))
    paren2_part1b <- tstar_d1_rho(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y])
    paren2_part1 <- paren2_part1a * paren2_part1b
    
    # 2nd parenthesis: part 2
    paren2_part2a <- phi_prime(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y-1L]))
    paren2_part2b <- tstar_d1_rho(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y-1L])
    paren2_part2 <- paren2_part2a * paren2_part2b
  }
  
  ## assemble the final term
  scale0 <- 1.0 / sqrt(sigma2 * (1 - rho^2))
  scale1 <- 1.0 + (rho^2) / (1.0 - rho^2)
  diff1  <- scale1 * (paren1_part1 - paren1_part2) 
  diff2  <- rho * (paren2_part1 - paren2_part2)
  return(scale0 * (diff1 + diff2))
}



py_x_d2_rhosigma2 <- function(y, x, rho, mu, sigma2, thres, num_y)
{
  if(y == 1L)
  {
    # 1st parenthesis
    paren1_part1 <- phi(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y]))
    paren1_part2 <- 0.0
    
    # 2nd parenthesis: part 1
    paren2_part1a <- phi_prime(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y]))
    paren2_part1b <- tstar_d1_rho(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y])
    paren2_part1 <- paren2_part1a * paren2_part1b
    
    # 2nd parenthesis: part 2
    paren2_part2 <- 0.0
    
  } else if(y == num_y)
  {
    # 1st parenthesis
    paren1_part1 <- 0.0
    paren1_part2 <- phi(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y-1L]))
    
    # 2nd parenthesis: part 1
    paren2_part1 <- 0.0
    
    # 2nd parenthesis: part 2
    paren2_part2a <- phi_prime(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y-1L]))
    paren2_part2b <- tstar_d1_rho(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y-1L])
    paren2_part2 <- paren2_part2a * paren2_part2b
    
  } else
  {
    # 1st parenthesis
    paren1_part1 <- phi(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y]))
    paren1_part2 <- phi(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y-1L]))
    
    # 2nd parenthesis: part 1
    paren2_part1a <- phi_prime(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y]))
    paren2_part1b <- tstar_d1_rho(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y])
    paren2_part1 <- paren2_part1a * paren2_part1b
    
    # 2nd parenthesis: part 2
    paren2_part2a <- phi_prime(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y-1L]))
    paren2_part2b <- tstar_d1_rho(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y-1L])
    paren2_part2 <- paren2_part2a * paren2_part2b
  }
  
  ## assemble the final term
  scale0 <- (x - mu) / (2.0 * sigma2^1.5 * sqrt(1.0 - rho^2))
  scale1 <- 1.0 + (rho^2) / (1.0 - rho^2)
  diff1  <- scale1 * (paren1_part1 - paren1_part2) 
  diff2  <- rho * (paren2_part1 - paren2_part2)
  return(scale0 * (diff1 + diff2))
}



py_x_d2_musigma2 <- function(y, x, rho, mu, sigma2, thres, num_y)
{
  if(y == 1L)
  {
    # 1st parenthesis: part 1
    paren1_part1a <- phi_prime(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y]))
    paren1_part1b <- tstar_d1_sigma2(
      x = x, rho = rho, mu = mu, sigma2 = sigma2)
    paren1_part1 <- paren1_part1a * paren1_part1b
    
    # 1st parenthesis: part 2
    paren1_part2 <- 0.0
    
    # 2nd parenthesis
    paren2_part1 <- phi(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y]))
    paren2_part2 <- 0.0
    
  } else if(y == num_y)
  {
    # 1st parenthesis: part 1
    paren1_part1 <- 0.0
    
    # 1st parenthesis: part 2
    paren1_part2a <- phi_prime(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y-1L]))
    paren1_part2b <- tstar_d1_sigma2(
      x = x, rho = rho, mu = mu, sigma2 = sigma2)
    paren1_part2 <- paren1_part2a * paren1_part2b
    
    # 2nd parenthesis
    paren2_part1 <- 0.0
    paren2_part2 <- phi(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y-1L]))
    
  } else
  {
    # 1st parenthesis: part 1
    paren1_part1a <- phi_prime(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y]))
    paren1_part1b <- tstar_d1_sigma2(
      x = x, rho = rho, mu = mu, sigma2 = sigma2)
    paren1_part1 <- paren1_part1a * paren1_part1b
    
    # 1st parenthesis: part 2
    paren1_part2a <- phi_prime(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y-1L]))
    paren1_part2b <- tstar_d1_sigma2(
      x = x, rho = rho, mu = mu, sigma2 = sigma2)
    paren1_part2 <- paren1_part2a * paren1_part2b
    
    # 2nd parenthesis
    paren2_part1 <- phi(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y]))
    paren2_part2 <- phi(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y-1L]))
  }
  
  ## assemble the final term
  scale0 <- rho / sqrt(sigma2 * (1.0 - rho^2))
  scale1 <- 1.0 / (2 * sigma2)
  diff1 <- paren1_part1 - paren1_part2
  diff2 <- scale1 * (paren2_part1 - paren2_part2)
  return(scale0 * (diff1 - diff2))
}


## defined for k \in [num_y-1]
py_x_d2_mutauk <- function(y, x, k, rho, mu, sigma2, thres, num_y)
{
  stopifnot(k < num_y)
  ym1     <- y - 1L
  scaling <- rho / (sqrt(sigma2) * (1.0 - rho^2))
  
  if(k == y)
  {
    tmp <- phi_prime(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y]))
    out <- scaling * tmp
    
  } else if(k == ym1)
  {
    tmp <- phi_prime(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[ym1]))
    out <- (-1.0) * scaling * tmp
  } else
  {
    out <- 0.0
  }
  return(out)
}


## defined for k \in [num_y-1]
py_x_d2_sigma2tauk <- function(y, x, k, rho, mu, sigma2, thres, num_y)
{
  stopifnot(k < num_y)
  ym1     <- y - 1L
  scaling <- rho * (x - mu) / (2.0 * sigma2^1.5 * (1.0 - rho^2))
  
  if(k == y)
  {
    tmp <- phi_prime(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[y]))
    out <- scaling * tmp
    
  } else if(k == ym1)
  {
    tmp <- phi_prime(get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[ym1]))
    out <- (-1.0) * scaling * tmp
  } else
  {
    out <- 0.0
  }
  return(out)
}


## defined for k \in [num_y-1]
py_x_d2_rhotauk <- function(y, x, k, rho, mu, sigma2, thres, num_y)
{
  stopifnot(k < num_y)
  ym1          <- y - 1L
  oneminusrho2 <- 1.0 - rho^2
  scaling1     <- rho / oneminusrho2^1.5
  scaling2     <- 1.0 / sqrt(oneminusrho2)
  
  if(k == y)
  {
     tstar <- get_taustar_y(
       x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[k])
     part1 <- phi(tstar)
     part2 <- phi_prime(tstar) * tstar_d1_rho(x = x, rho = rho, mu = mu,
                                              sigma2 = sigma2, tau_y = thres[k])
     out <- scaling1 * part1 + scaling2 * part2
  } else if(k == ym1)
  {
    tstar <- get_taustar_y(
      x = x, rho = rho, mu = mu, sigma2 = sigma2, tau_y = thres[k])
    part1 <- phi(tstar)
    part2 <- phi_prime(tstar) * tstar_d1_rho(x = x, rho = rho, mu = mu,
                                             sigma2 = sigma2, tau_y = thres[k])
    out <- (-1.0) * (scaling1 * part1 + scaling2 * part2)
  } else
  {
    out <- 0.0
  }
  return(out)
}
