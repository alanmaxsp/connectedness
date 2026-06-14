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

test_that("Schur auto selects CHOLMOD for sparse custom relationship", {
  data <- data.frame(
    animal_id = c("A", "B", "C", "D"),
    region    = c("MU1", "MU1", "MU2", "MU2"),
    sex       = c("F", "M", "F", "M"),
    stringsAsFactors = FALSE
  )

  Kinv <- Matrix::Diagonal(4)
  animal_index <- setNames(seq_len(4), data$animal_id)

  res <- compute_connectedness(
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
    schur_solver = "auto",
    verbose = FALSE
  )

  expect_equal(res$schur_solver, "cholmod")
})

test_that("Schur dense solver reproduces full MME backend for dense custom relationship", {
  data <- data.frame(
    animal_id = c("A", "B", "C", "D", "E", "F"),
    region    = c("MU1", "MU1", "MU2", "MU2", "MU3", "MU3"),
    sex       = c("F", "M", "F", "M", "F", "M"),
    stringsAsFactors = FALSE
  )

  Kinv <- diag(6)
  animal_index <- setNames(seq_len(6), data$animal_id)

  res_full <- compute_connectedness(
    data = data,
    animal_col = "animal_id",
    mu_col = "region",
    fixed_formula = ~ 1 + sex,
    sigma2a = 1.2,
    sigma2e = 2.4,
    relationship = "custom",
    rel_matrix = Kinv,
    animal_index = animal_index,
    mme_backend = "full_mme",
    verbose = FALSE
  )

  res_dense <- compute_connectedness(
    data = data,
    animal_col = "animal_id",
    mu_col = "region",
    fixed_formula = ~ 1 + sex,
    sigma2a = 1.2,
    sigma2e = 2.4,
    relationship = "custom",
    rel_matrix = Kinv,
    animal_index = animal_index,
    mme_backend = "schur",
    schur_solver = "dense",
    verbose = FALSE
  )

  expect_equal(res_dense$schur_solver, "dense")
  expect_equal(res_dense$CD, res_full$CD, tolerance = 1e-8)
  expect_equal(res_dense$PEVD, res_full$PEVD, tolerance = 1e-8)
  expect_equal(res_dense$qC, res_full$qC, tolerance = 1e-8)
  expect_equal(res_dense$qK, res_full$qK, tolerance = 1e-8)
})

test_that("Schur auto selects dense for dense custom relationship", {
  data <- data.frame(
    animal_id = c("A", "B", "C", "D"),
    region    = c("MU1", "MU1", "MU2", "MU2"),
    sex       = c("F", "M", "F", "M"),
    stringsAsFactors = FALSE
  )

  Kinv <- diag(4)
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
    schur_solver = "auto",
    dry_run = TRUE,
    verbose = FALSE
  )

  expect_equal(diag$matrix_storage, "dense")
  expect_equal(diag$selected_schur_solver, "dense")
})

test_that("Schur dense solver reproduces full MME for non-diagonal dense custom relationship", {
  data <- data.frame(
    animal_id = c("A", "B", "C", "D", "E", "F"),
    region    = c("MU1", "MU1", "MU2", "MU2", "MU3", "MU3"),
    sex       = c("F", "M", "F", "M", "F", "M"),
    stringsAsFactors = FALSE
  )

  set.seed(123)
  M <- matrix(rnorm(36), 6, 6)
  Kinv <- crossprod(M) + diag(6) * 2
  animal_index <- setNames(seq_len(6), data$animal_id)

  res_full <- compute_connectedness(
    data = data,
    animal_col = "animal_id",
    mu_col = "region",
    fixed_formula = ~ 1 + sex,
    sigma2a = 1.2,
    sigma2e = 2.4,
    relationship = "custom",
    rel_matrix = Kinv,
    animal_index = animal_index,
    mme_backend = "full_mme",
    verbose = FALSE
  )

  res_dense <- compute_connectedness(
    data = data,
    animal_col = "animal_id",
    mu_col = "region",
    fixed_formula = ~ 1 + sex,
    sigma2a = 1.2,
    sigma2e = 2.4,
    relationship = "custom",
    rel_matrix = Kinv,
    animal_index = animal_index,
    mme_backend = "schur",
    schur_solver = "dense",
    verbose = FALSE
  )

  expect_equal(res_dense$CD, res_full$CD, tolerance = 1e-8)
  expect_equal(res_dense$PEVD, res_full$PEVD, tolerance = 1e-8)
  expect_equal(res_dense$qC, res_full$qC, tolerance = 1e-8)
  expect_equal(res_dense$qK, res_full$qK, tolerance = 1e-8)
})

test_that("CHOLMOD low-memory solver reproduces standard CHOLMOD solver", {
  data <- data.frame(
    animal_id = c("A", "B", "C", "D", "E", "F"),
    region    = c("MU1", "MU1", "MU2", "MU2", "MU3", "MU3"),
    sex       = c("F", "M", "F", "M", "F", "M"),
    year      = c(2020, 2020, 2021, 2021, 2022, 2022),
    stringsAsFactors = FALSE
  )

  Kinv <- Matrix::Diagonal(6)
  animal_index <- setNames(seq_len(6), data$animal_id)

  res_cholmod <- compute_connectedness(
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
    schur_solver = "cholmod",
    verbose = FALSE
  )

  res_lowmem <- compute_connectedness(
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
    schur_solver = "cholmod_lowmem",
    schur_x_block_size = 1L,
    verbose = FALSE
  )

  expect_equal(res_lowmem$schur_solver, "cholmod_lowmem")
  expect_equal(res_lowmem$CD, res_cholmod$CD, tolerance = 1e-8)
  expect_equal(res_lowmem$PEVD, res_cholmod$PEVD, tolerance = 1e-8)
  expect_equal(res_lowmem$qC, res_cholmod$qC, tolerance = 1e-8)
  expect_equal(res_lowmem$qK, res_cholmod$qK, tolerance = 1e-8)
})

test_that("dry_run reports low-memory Schur diagnostics", {
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
    schur_solver = "cholmod_lowmem",
    schur_x_block_size = 1L,
    dry_run = TRUE,
    verbose = FALSE
  )

  expect_equal(diag$selected_schur_solver, "cholmod_lowmem")
  expect_true(all(c(
    "schur_W_elements",
    "schur_W_exceeds_32bit",
    "schur_lowmem_X_block_size",
    "schur_lowmem_W_block_gib"
  ) %in% names(diag)))
  expect_equal(diag$schur_lowmem_X_block_size, 1L)
})

test_that("auto can select CHOLMOD low-memory when W is too large", {
  Kinv <- Matrix::Diagonal(2)
  Xsp <- Matrix::sparseMatrix(i = integer(), j = integer(), dims = c(1, 2))

  expect_equal(
    .choose_schur_solver(
      rel_matrix = Kinv,
      relationship = "custom",
      schur_solver = "auto",
      fixed_matrix = Xsp,
      lowmem_W_gib_threshold = 0
    ),
    "cholmod_lowmem"
  )
})
