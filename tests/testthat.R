if(requireNamespace("testthat", quietly = TRUE)) {
  library("robord", quietly = TRUE)
  testthat::test_check("robord")
}
