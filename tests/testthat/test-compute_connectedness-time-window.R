test_that("min_records_per_year filters MUs active in at least half of year_window", {
  data <- data.frame(
    animal_id = paste0("A", 1:12),
    region = c(
      rep("MU1", 6),
      rep("MU2", 4),
      rep("MU3", 2)
    ),
    year = c(
      2020, 2020, 2021, 2021, 2022, 2022,
      2020, 2020, 2021, 2021,
      2020, 2021
    ),
    sex = rep(c("F", "M"), 6),
    stringsAsFactors = FALSE
  )

  Kinv <- Matrix::Diagonal(nrow(data))
  animal_index <- setNames(seq_len(nrow(data)), data$animal_id)

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
    year_col = "year",
    year_window = c(2020, 2022),
    min_records_per_year = 2,
    mme_backend = "schur",
    verbose = FALSE
  )

  expect_true(all(c("MU1", "MU2") %in% rownames(res$CD)))
  expect_false("MU3" %in% rownames(res$CD))
  expect_false(is.null(res$activity_summary))
  expect_equal(sum(res$activity_summary$eligible), 2)
  expect_equal(res$activity_summary$n_active_years[res$activity_summary$MU == "MU3"], 0)
})

test_that("min_records_per_year NULL does not filter MUs inside year_window", {
  data <- data.frame(
    animal_id = paste0("A", 1:12),
    region = c(
      rep("MU1", 6),
      rep("MU2", 4),
      rep("MU3", 2)
    ),
    year = c(
      2020, 2020, 2021, 2021, 2022, 2022,
      2020, 2020, 2021, 2021,
      2020, 2021
    ),
    sex = rep(c("F", "M"), 6),
    stringsAsFactors = FALSE
  )

  Kinv <- Matrix::Diagonal(nrow(data))
  animal_index <- setNames(seq_len(nrow(data)), data$animal_id)

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
    year_col = "year",
    year_window = c(2020, 2022),
    min_records_per_year = NULL,
    mme_backend = "schur",
    verbose = FALSE
  )

  expect_true("MU3" %in% rownames(res$CD))
  expect_null(res$activity_summary)
})

test_that("temporal activity filter errors when fewer than two MUs remain", {
  data <- data.frame(
    animal_id = paste0("A", 1:12),
    region = c(
      rep("MU1", 6),
      rep("MU2", 4),
      rep("MU3", 2)
    ),
    year = c(
      2020, 2020, 2021, 2021, 2022, 2022,
      2020, 2020, 2021, 2021,
      2020, 2021
    ),
    sex = rep(c("F", "M"), 6),
    stringsAsFactors = FALSE
  )

  Kinv <- Matrix::Diagonal(nrow(data))
  animal_index <- setNames(seq_len(nrow(data)), data$animal_id)

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
      year_col = "year",
      year_window = c(2020, 2022),
      min_records_per_year = 3,
      mme_backend = "schur",
      verbose = FALSE
    ),
    "Fewer than 2 MUs remain"
  )
})

test_that("year_window selects reported MUs while select_mus uses all records", {
  data <- data.frame(
    animal_id = paste0("A", 1:16),
    region = c(
      rep("MU1", 6),
      rep("MU2", 4),
      rep("MU3", 2),
      "MU1", "MU2", "MU3", "MU4"
    ),
    year = c(
      2020, 2020, 2021, 2021, 2022, 2022,
      2020, 2020, 2021, 2021,
      2020, 2021,
      2019, 2019, 2019, 2019
    ),
    sex = rep(c("F", "M"), 8),
    stringsAsFactors = FALSE
  )

  Kinv <- Matrix::Diagonal(nrow(data))
  animal_index <- setNames(seq_len(nrow(data)), data$animal_id)

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
    year_col = "year",
    year_window = c(2020, 2022),
    min_records_per_year = 2,
    mme_backend = "schur",
    verbose = FALSE
  )

  expect_equal(res$temporal_mode, "select_mus")
  expect_equal(rownames(res$CD), c("MU1", "MU2"))
  expect_equal(as.numeric(res$n_target), c(7, 5))
  expect_equal(res$report_mus, c("MU1", "MU2"))
  expect_true(all(res$overlap$MU1 %in% c("MU1", "MU2")))
  expect_true(all(res$overlap$MU2 %in% c("MU1", "MU2")))

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
    year_col = "year",
    year_window = c(2020, 2022),
    min_records_per_year = 2,
    mme_backend = "schur",
    dry_run = TRUE,
    verbose = FALSE
  )

  expect_equal(diag$n_records_activity_window, 12)
  expect_equal(diag$n_records_used_for_connectedness, 16)
  expect_equal(diag$n_mus_used_for_connectedness_records, 4)
  expect_equal(diag$n_management_units, 2)
})

test_that("temporal_mode filter_records keeps legacy record-filtered behavior", {
  data <- data.frame(
    animal_id = paste0("A", 1:16),
    region = c(
      rep("MU1", 6),
      rep("MU2", 4),
      rep("MU3", 2),
      "MU1", "MU2", "MU3", "MU4"
    ),
    year = c(
      2020, 2020, 2021, 2021, 2022, 2022,
      2020, 2020, 2021, 2021,
      2020, 2021,
      2019, 2019, 2019, 2019
    ),
    sex = rep(c("F", "M"), 8),
    stringsAsFactors = FALSE
  )

  Kinv <- Matrix::Diagonal(nrow(data))
  animal_index <- setNames(seq_len(nrow(data)), data$animal_id)

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
    year_col = "year",
    year_window = c(2020, 2022),
    min_records_per_year = 2,
    temporal_mode = "filter_records",
    mme_backend = "schur",
    verbose = FALSE
  )

  expect_equal(res$temporal_mode, "filter_records")
  expect_equal(rownames(res$CD), c("MU1", "MU2"))
  expect_equal(as.numeric(res$n_target), c(6, 4))
})
