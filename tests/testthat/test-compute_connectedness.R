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
  
  out <- compute_connectedness(
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

test_that("compute_connectedness dry_run devuelve diagnosticos sin resolver MME", {
  data <- data.frame(
    animal_id = c("A", "B", "C", "D"),
    region    = c("MU1", "MU1", "MU2", "MU2"),
    sex       = c("F", "M", "F", "M"),
    stringsAsFactors = FALSE
  )

  Kinv <- Matrix::Diagonal(4)
  animal_index <- setNames(1:4, c("A", "B", "C", "D"))

  out <- compute_connectedness(
    data = data,
    animal_col = "animal_id",
    mu_col = "region",
    fixed_formula = ~ 1 + sex,
    sigma2a = 1,
    sigma2e = 1,
    relationship = "custom",
    rel_matrix = Kinv,
    animal_index = animal_index,
    dry_run = TRUE,
    verbose = FALSE
  )

  expect_s3_class(out, "connectedness_diagnostics")
  expect_equal(out$n_relationship, 4)
  expect_equal(out$n_records, 4)
  expect_equal(out$n_management_units, 2)
  expect_equal(out$mme_dim, out$n_relationship + out$n_fixed_effect_columns)
})

test_that("compute_connectedness detiene problemas que superan max_mme_dim", {
  data <- data.frame(
    animal_id = c("A", "B", "C", "D"),
    region    = c("MU1", "MU1", "MU2", "MU2"),
    sex       = c("F", "M", "F", "M"),
    stringsAsFactors = FALSE
  )

  Kinv <- Matrix::Diagonal(4)
  animal_index <- setNames(1:4, c("A", "B", "C", "D"))

  expect_error(
    compute_connectedness(
      data = data,
      animal_col = "animal_id",
      mu_col = "region",
      fixed_formula = ~ 1 + sex,
      sigma2a = 1,
      sigma2e = 1,
      relationship = "custom",
      rel_matrix = Kinv,
      animal_index = animal_index,
      max_mme_dim = 3,
      verbose = FALSE
    ),
    "exceeds max_mme_dim"
  )
})
