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

## replace 20% of observations with negative leverage points
x[1:200] <- 1
y[1:200] <- Ky

## MLE
mle <- polycor_mle(x = x, y = y)
mle$thetahat 
# > mle$thetahat
#         rho          a1          a2          a3          a4          b1          b2          b3          b4 
# -0.34675954 -0.63244517 -0.39741400  0.10278048  0.93030935 -1.57214524 -1.12479616  0.03080319  0.63796166 

## robust
polycor <- polycor(x = x, y = y)
polycor$thetahat
# > polycor$thetahat
#       rho         a1         a2         a3         a4         b1         b2         b3         b4 
# 0.3170347 -1.4412686 -0.9580768 -0.2337379  0.7789224 -1.5291725 -0.9982777  0.4074513  1.4534537  


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

p <- plot(polycor)
dir_save <- "inst/doc/readme_plots"
ggsave(filename = "READMEplot.svg", device = "svg", path = dir_save, width = 6, height = 6, plot = p)

celltest(polycor)$pval_adjusted
