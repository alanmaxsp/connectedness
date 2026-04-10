test_that("build_Ginv devuelve Ginv cuadrada para X válida", {
  set.seed(1)
  X <- matrix(
    sample(c(0L, 1L, 2L), 40, replace = TRUE),
    nrow = 8, ncol = 5
  )
  
  out <- connectedness::build_Ginv(
    X = X,
    maf_threshold = 0,
    blend = 0.05,
    n_threads = 1L,
    tunedG = 0L
  )
  
  expect_type(out, "list")
  expect_true("Ginv" %in% names(out))
  expect_equal(dim(out$Ginv), c(nrow(X), nrow(X)))
  expect_true(is.numeric(out$Ginv))
})

test_that("build_Ginv falla si X tiene NA", {
  X <- matrix(c(0,1,2,NA,1,0), nrow = 3)
  expect_error(
    connectedness::build_Ginv(X),
    "must not contain NA values"
  )
})

test_that("build_Ginv exige A22 cuando tunedG es 2 o 3", {
  X <- matrix(sample(c(0L,1L,2L), 30, replace = TRUE), nrow = 6)
  expect_error(
    connectedness::build_Ginv(X, tunedG = 2L, maf_threshold = 0),
    "A22.*must be supplied"
  )
})