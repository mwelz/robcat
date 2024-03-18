if(requireNamespace("testthat", quietly = TRUE)) {
  library("robcat", quietly = TRUE)
  testthat::test_check("robcat")
}
