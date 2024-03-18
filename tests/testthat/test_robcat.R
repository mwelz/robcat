set.seed(1)
N <- 1000L
p <- 3L
K <- 5L
thres <- c(-Inf, -1.5, -1, -0.25, 0.75, Inf)
latent <- mvtnorm::rmvnorm(N, sigma = diag(p))
data <- sapply(seq_len(p), function(j) as.integer(cut(latent[,j], thres)))

mat_par <- robcat::polycormat(data, parallel = FALSE, constrained = FALSE)
san1 <- polycor(data[,2], data[,1], constrained = FALSE)
san2 <- polycor(data[,3], data[,1], constrained = FALSE)
san3 <- polycor(data[,3], data[,2], constrained = FALSE)

test_that("Input order shouldn't affect estimates, only labeling of thresholds", {
  tmp <-  polycor(data[,1], data[,2], constrained = FALSE)
  # flip thresholds
  hat <- c(tmp$thetahat[1], tmp$thetahat[paste0("b", 1:(K-1))], tmp$thetahat[paste0("a", 1:(K-1))])
  expect_equal(unname(san1$thetahat), unname(hat))
})

test_that("equivalence between table and vector input", {
  tmp <-  polycor(x = table(data[,2], data[,1]), constrained = FALSE)
  expect_equal(tmp$thetahat, san1$thetahat)
})

test_that("Correlation matrix computation behaves as expected", {
  expect_equal(mat_par$cormat[1,2], unname(san1$thetahat[1]))
  expect_equal(mat_par$cormat[1,3], unname(san2$thetahat[1]))
  expect_equal(mat_par$cormat[3,2], unname(san3$thetahat[1]))
})
