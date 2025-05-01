# robcat: Robust Categorical Data Analysis

[![CRAN](https://www.R-pkg.org/badges/version/robcat)](https://CRAN.R-project.org/package=robcat) 

This package implements the methodology proposed in the working paper _Robust Estimation and Inference in Categorical Data_ by [Welz (2024)](https://arxiv.org/abs/2403.11954). Here we will demonstrate how the methodology can be used to robustly estimate polychoric correlation, which is described in detail in our companion paper _Robust Estimation of Polychoric Correlation_ by [Welz, Mair, and Alfons (2024)](https://arxiv.org/abs/2407.18835).

Package `robcat` is on CRAN (The Comprehensive R Archive Network), hence the latest release can be easily installed from the `R` command line via

```R
install.packages("robcat")
```

To install the latest development version from GitHub, you can pull this repository and install it from the `R` command line via
```R
install.packages("devtools")
devtools::install_github("mwelz/robcat")
```
If you already have the package `devtools` installed, you can skip the first line.


## Example of robust estimation of polychoric correlation coefficient

### Generate simulated data

```R
library("robcat")

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
mle <- polycor_mle(x = x, y = y)
mle$thetahat 
# > mle$thetahat
#       rho         a1         a2         a3         a4         b1         b2         b3         b4 
# 0.3151109 -1.4868862 -0.9829336 -0.2318141  0.7887485 -1.5112743 -0.9889993  0.4276572  1.4582239 

## robust
polycor <- polycor(x = x, y = y)
polycor$thetahat
# > polycor$thetahat
#       rho         a1         a2         a3         a4         b1         b2         b3         b4 
# 0.3151731 -1.4867730 -0.9828117 -0.2317529  0.7887506 -1.5110644 -0.9888898  0.4276510  1.4583212 
```

Thus, in the absence of contamination, both estimators yield equivalent solutions. Next, we introduce 20% contamination.

### Compare MLE and robust estimator with contamination

```R
## replace 20% of observations with negative leverage points
x[1:200] <- 1
y[1:200] <- Ky

## MLE
mle <- polycor_mle(x = x, y = y)
mle$thetahat 
# > mle$thetahat
#         rho          a1          a2          a3          a4          b1          b2          b3          b4 
# -0.34689341 -0.63243705 -0.39737803  0.10283550  0.93007017 -1.57234648 -1.12487769  0.03064128  0.63781488 

## robust
polycor <- polycor(x = x, y = y)
polycor$thetahat
# > polycor$thetahat
#       rho         a1         a2         a3         a4         b1         b2         b3         b4 
# 0.3170347 -1.4412686 -0.9580768 -0.2337379  0.7789224 -1.5291725 -0.9982777  0.4074513  1.4534537

```

We see that 20% contamination leads to a substantial bias in the MLE, whereas the robust estimator is still accurate.
The package also provides methods for printing and plotting:

```R
## print and plot method
polycor
# > polycor
# 
# Polychoric Correlation
# Estimate Std.Err.
# rho    0.317  0.03892
# 
# X-thresholds
#     Estimate Std.Err.
# a1  -1.4410  0.06601
# a2  -0.9581  0.05252
# a3  -0.2337  0.04446
# a4   0.7789  0.04958
# 
# Y-thresholds
#     Estimate Std.Err.
# b1  -1.5290  0.06921
# b2  -0.9983  0.05308
# b3   0.4075  0.04536
# b4   1.4530  0.06761

plot(polycor)
```

<img src="./inst/doc/readme_plots/READMEplot.svg" width="67%" style="display: block; margin: auto;" />

Indeed, the Pearson residual of contaminated cell `(x,y) = (1,5)` is excessively large compared to the others, which are all around the value 0.


## Report issues and request features

If you experience any bugs or issues or if you have any suggestions for additional features, please submit an issue via the [*Issues*](https://github.com/mwelz/robcat/issues) tab of this repository.  Please have a look at existing issues first to see if your problem or feature request has already been discussed.


## Contribute to the package

If you want to contribute to the package, you can fork this repository and create a pull request after implementing the desired functionality.


## Ask for help

If you need help using the package, or if you are interested in collaborations related to this project, please get in touch with the [package maintainer](https://mwelz.github.io/).
