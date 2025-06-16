# load preprocessed data
load("../polycorstuff/applications/arias2020/data_preprocessed.Rdata")

# replace "S" in variable names with "N" 
# (Arias et al., 2020, refer to neuroticism trait as emotional stability)
colnames(responses) <- gsub("S", "N", colnames(responses))


# tuning parameter for robust estimator
c <- 0.6


# keep only items for neuroticism
items <- names(responses)
keep <- which(startsWith(items, "N"))
responses_N <- responses[, keep]  # no NAs in Big Five (only in Greenleaf scale)

# undo reverse coding of negatively worded items
items_N <- names(responses_N)
neg_worded <- which(endsWith(items_N, "N"))
for(j in neg_worded) {
  reversed <- responses_N[, j]
  responses_N[, j] <- abs(reversed - 5L - 1L)  # 5 answer categories
}

i <- 1
j <- 1

x <- responses_N[, i]
y <- responses_N[, j]

robcat::polycor(x = x, y = y, c = 0.6, variance = FALSE,
                constrained = FALSE)
robcat::polycor(x = x, y = y, c = 0.6, variance = FALSE,
                constrained = FALSE, method = "Nelder-Mead")
robcat::polycor(x = x, y = y, c = 0.6, variance = FALSE,
                constrained = TRUE)
robcat::polycor_mle(x = x, y = y, variance = FALSE, 
                    twostep = TRUE)
robcat::polycor_mle(x = x, y = y, variance = FALSE, 
                    constrained = FALSE)

microbenchmark::microbenchmark(
  "L-BFGS-B" = robcat::polycor(x = x, y = y, c = Inf, variance = FALSE,
                               constrained = FALSE),
  "Nelder-Mead" = robcat::polycor(x = x, y = y, c = Inf, variance = FALSE,
                               constrained = FALSE, method = "Nelder-Mead"),
  times = 10
)
