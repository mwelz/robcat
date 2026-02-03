## point polyserial coefficient (calculate from polyserial correlation coefficient)

y_mean <- function(tau, ycat = 1:(length(tau)+1))
{
  py_pmf <- py(tau)
  sum(ycat * py_pmf)
}

y_var <- function(tau, ycat = 1:(length(tau)+1))
{
  m <- y_mean(tau = tau, ycat = ycat)
  m2 <- y_mean(tau = tau, ycat = ycat^2)
  m2 - m^2
}


xy_mean <- function(theta, ycat = 1:(length(theta)-2))
{
  ymax <- max(ycat)
  tau <- theta[startsWith(names(theta), "tau")]
  rho <- theta[1]
  mu <- theta[2]
  sigma <- sqrt(theta[3])
  numy <- length(ycat)
  
  Phisummands <- sapply(seq_len(numy - 1L), function(j){
    stats::pnorm(tau[j]) * (ycat[j+1L] - ycat[j])
  })
  
  phisummands <- sapply(seq_len(numy - 1L), function(j){
    stats::dnorm(tau[j]) * (ycat[j+1L] - ycat[j])
  })
  
  part1 <- mu * (ymax - sum(Phisummands))
  part2 <- rho * sigma * sum(phisummands)
  return(unname(part1 + part2))
}


pointpolyserial <- function(theta, ycat = 1:(length(theta)-2))
{
  theta <- as.numeric(theta)
  names(theta) <- polyserial_thetanames(max(ycat))
  tau <- theta[startsWith(names(theta), "tau")]
  mux <- theta[2]
  sigmax <- sqrt(theta[3])
  
  muxy <- xy_mean(theta = theta, ycat = ycat)
  muy <- y_mean(tau = tau, ycat = ycat)
  sigmay <- sqrt(y_var(tau = tau, ycat = ycat))
  
  numer <- muxy - mux * muy
  denom <- sigmax * sigmay
  
  return(unname(numer / denom))
}