# robord: Robust Ordinal Data Analysis

Example of robust estimation of polychoric correlation coefficient:

```R
## 5 answer categories each, define latent thresholds as follows
Kx <- Ky <- 5
thresX <- c(-Inf, -1.5, -1, -0.25, .75, Inf)
thresY <- c(-Inf, -1.5, -1, .5, 1.5, Inf)
rho_true <- 0.3 # true polychoric correlation

## simulate rating data
set.seed(20240111)
latent <- mvtnorm::rmvnorm(1000, c(0, 0), matrix(c(1, rho_true, rho_true, 1), 2, 2))
xi <- latent[,1]
eta <- latent[,2]
x <- as.integer(cut(xi, thresX))
y <- as.integer(cut(eta, thresY))

## robust estimation of rho and thresholds (no contamination here)
polycor <- robord::polycor(x = x, y = y, c = 1.5)
polycor$thetahat ## accurate
# > polycor$thetahat
#       rho         a1         a2         a3         a4         b1         b2         b3         b4 
# 0.3151414 -1.4868373 -0.9829001 -0.2316910  0.7888495 -1.5112325 -0.9889749  0.4277239  1.4582497

## MLE
mle <- robord::polycor_mle(x = x, y = y)
mle$thetahat ## accurate and equivalent to robust estimation (as expected)
# > mle$thetahat
#       rho         a1         a2         a3         a4         b1         b2         b3         b4 
# 0.3151109 -1.4868862 -0.9829336 -0.2318141  0.7887485 -1.5112743 -0.9889993  0.4276572  1.4582239 
```


## Authors
Max Welz (welz@ese.eur.nl)
