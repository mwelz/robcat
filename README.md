# robord: Robust Ordinal Data Analysis

## Example of robust estimation of polychoric correlation coefficient

### Generate simulated data

```R
## 5 answer categories each, define latent thresholds as follows
Kx <- Ky <- 5
thresX <- c(-Inf, -1.5, -1, -0.25, 0.75, Inf)
thresY <- c(-Inf, -1.5, -1, 0.5, 1.5, Inf)
rho_true <- 0.3 # true polychoric correlation

## simulate rating data
set.seed(20240111)
latent <- mvtnorm::rmvnorm(1000, c(0, 0), matrix(c(1, rho_true, rho_true, 1), 2, 2))
xi <- latent[,1]
eta <- latent[,2]
x <- as.integer(cut(xi, thresX))
y <- as.integer(cut(eta, thresY))
```

### Compare MLE and robust estimator without contamination

```R
## MLE
mle <- robord::polycor_mle(x = x, y = y)
mle$thetahat 
# > mle$thetahat
#       rho         a1         a2         a3         a4         b1         b2         b3         b4 
# 0.3151109 -1.4868862 -0.9829336 -0.2318141  0.7887485 -1.5112743 -0.9889993  0.4276572  1.4582239 

## robust
polycor <- robord::polycor(x = x, y = y, c = 1.5)
polycor$thetahat
# > polycor$thetahat
#       rho         a1         a2         a3         a4         b1         b2         b3         b4 
# 0.3151414 -1.4868373 -0.9829001 -0.2316910  0.7888495 -1.5112325 -0.9889749  0.4277239  1.4582497
```

Thus, in the absence of contamination, both estimators yield equivalent solutions. Next, we introduce 20% contamination.

### Compare MLE and robust estimator with contamination

```R
## replace 20% of observations with negative leverage points
x[1:200] <- 1
y[1:200] <- Ky

## MLE
mle <- robord::polycor_mle(x = x, y = y)
mle$thetahat 
# > mle$thetahat
#         rho          a1          a2          a3          a4          b1          b2          b3          b4 
# -0.34675954 -0.63244517 -0.39741400  0.10278048  0.93030935 -1.57214524 -1.12479616  0.03080319  0.63796166 

## robust
polycor <- robord::polycor(x = x, y = y, c = 1.5)
polycor$thetahat
# > polycor$thetahat
#       rho         a1         a2         a3         a4         b1         b2         b3         b4 
# 0.3180104 -1.4461457 -0.9605778 -0.2342293  0.7795890 -1.5299883 -0.9981569  0.4092214  1.4566111 

```

We see that 20% contamination leads to a substantial bias in the MLE, whereas the robust estimator is still accurate.

## Authors
Max Welz (welz@ese.eur.nl)
