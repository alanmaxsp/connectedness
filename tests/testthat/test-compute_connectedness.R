test_that("compute_connectedness con relación custom devuelve objeto connectedness", {
  data <- data.frame(
    animal_id = c("A", "B", "C", "D"),
    region    = c("MU1", "MU1", "MU2", "MU2"),
    sex       = c("F", "M", "F", "M"),
    stringsAsFactors = FALSE
  )
  
  # Inversa kernel simple 4x4 (identidad)
  Kinv <- Matrix::Diagonal(4)
  animal_index <- setNames(1:4, c("A", "B", "C", "D"))
  
  out <- connectedness::compute_connectedness(
    data = data,
    animal_col = "animal_id",
    mu_col = "region",
    fixed_formula = ~ 1 + sex,
    sigma2a = 1,
    sigma2e = 1,
    relationship = "custom",
    rel_matrix = Kinv,
    animal_index = animal_index,
    verbose = FALSE
  )
  
  expect_s3_class(out, "connectedness")
  expect_true(all(c("CD","PEVD","qK","qC","n_target","relationship") %in% names(out)))
  expect_equal(out$relationship, "custom")
})