test_that("Schur backend reproduces full MME backend for custom relationship", {
  data <- data.frame(
    animal_id = c("A", "B", "C", "D", "E", "F"),
    region    = c("MU1", "MU1", "MU2", "MU2", "MU3", "MU3"),
    sex       = c("F", "M", "F", "M", "F", "M"),
    year      = c(2020, 2020, 2021, 2021, 2022, 2022),
    stringsAsFactors = FALSE
  )

  Kinv <- Matrix::Diagonal(6)
  animal_index <- setNames(seq_len(6), data$animal_id)

  res_full <- compute_connectedness(
    data = data,
    animal_col = "animal_id",
    mu_col = "region",
    fixed_formula = ~ 1 + sex + year,
    sigma2a = 1.2,
    sigma2e = 2.4,
    relationship = "custom",
    rel_matrix = Kinv,
    animal_index = animal_index,
    mme_backend = "full_mme",
    verbose = FALSE
  )

  res_schur <- compute_connectedness(
    data = data,
    animal_col = "animal_id",
    mu_col = "region",
    fixed_formula = ~ 1 + sex + year,
    sigma2a = 1.2,
    sigma2e = 2.4,
    relationship = "custom",
    rel_matrix = Kinv,
    animal_index = animal_index,
    mme_backend = "schur",
    verbose = FALSE
  )

  expect_equal(res_schur$CD, res_full$CD, tolerance = 1e-8)
  expect_equal(res_schur$PEVD, res_full$PEVD, tolerance = 1e-8)
  expect_equal(res_schur$qC, res_full$qC, tolerance = 1e-8)
  expect_equal(res_schur$qK, res_full$qK, tolerance = 1e-8)
  expect_equal(dimnames(res_schur$CD), dimnames(res_full$CD))
  expect_equal(names(res_schur$n_target), names(res_full$n_target))
})

test_that("dry_run reports Schur diagnostics", {
  data <- data.frame(
    animal_id = c("A", "B", "C", "D"),
    region    = c("MU1", "MU1", "MU2", "MU2"),
    sex       = c("F", "M", "F", "M"),
    stringsAsFactors = FALSE
  )

  Kinv <- Matrix::Diagonal(4)
  animal_index <- setNames(seq_len(4), data$animal_id)

  diag <- compute_connectedness(
    data = data,
    animal_col = "animal_id",
    mu_col = "region",
    fixed_formula = ~ 1 + sex,
    sigma2a = 1,
    sigma2e = 1,
    relationship = "custom",
    rel_matrix = Kinv,
    animal_index = animal_index,
    mme_backend = "schur",
    dry_run = TRUE,
    verbose = FALSE
  )

  expect_s3_class(diag, "connectedness_diagnostics")
  expect_true(all(c("schur_W_storage_mb", "schur_S_storage_mb", "requested_backend") %in% names(diag)))
  expect_equal(diag$requested_backend, "schur")
})
