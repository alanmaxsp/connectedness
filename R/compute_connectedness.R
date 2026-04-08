#' Compute connectedness metrics between management units
#'
#' Computes connectedness between management units (MUs) under the contrast
#' approach via Mixed Model Equations (MME). The function returns two
#' complementary metrics for all MU pairs: the Coefficient of Determination
#' (`CD`) and the Prediction Error Variance of Differences (`PEVD`).
#'
#' The connectedness analysis can be based on a pedigree-derived inverse
#' relationship matrix (`Ainv`), a genomic inverse relationship matrix (`Ginv`),
#' a single-step inverse relationship matrix (`Hinv`), or a user-supplied
#' inverse kernel.
#'
#' @param data A data frame containing the records used to define the analysis.
#'   It must include the animal identifier, the management-unit assignment, and
#'   all variables referenced in `fixed_formula`.
#' @param animal_col Character string giving the name of the animal ID column in
#'   `data`.
#' @param mu_col Character string giving the name of the management-unit column
#'   in `data`. Each animal must belong to exactly one MU.
#' @param fixed_formula A one-sided model formula describing the fixed-effects
#'   design matrix, for example `~ 1 + sex + birth_year`.
#' @param sigma2a Positive numeric scalar giving the additive genetic variance
#'   (or, more generally, the variance associated with the animal effect under
#'   the chosen kernel).
#' @param sigma2e Positive numeric scalar giving the residual variance.
#' @param relationship Character string indicating which inverse relationship
#'   matrix should underlie the connectedness analysis. Must be one of
#'   `"Ainv"`, `"Ginv"`, `"Hinv"`, or `"custom"`.
#' @param pedigree Optional pedigree data frame with three columns representing
#'   animal, sire, and dam. Required for `relationship = "Ainv"` and
#'   `relationship = "Hinv"`. It can also be used with `relationship = "Ginv"`
#'   when `animal_index` is derived from `genotyped_idx`.
#' @param X Optional genotype matrix (`n_gen x m`) coded as 0/1/2. Required for
#'   `relationship = "Ginv"` and `relationship = "Hinv"`.
#' @param genotyped_idx Optional integer vector of 1-based renumbered pedigree
#'   indices for the genotyped animals, in the same order as the rows of `X`.
#'   Required for `relationship = "Hinv"`, and also for `relationship = "Ginv"`
#'   when `animal_index` is not supplied directly.
#' @param animal_index Optional named integer vector mapping original animal IDs
#'   to row/column indices in the inverse relationship matrix actually used in
#'   the analysis. Required for `relationship = "custom"`; optional otherwise if
#'   it can be derived internally.
#' @param rel_matrix Optional user-supplied inverse relationship matrix of class
#'   `dgCMatrix`. Used only when `relationship = "custom"`.
#' @param maf_threshold Minor allele frequency threshold used when building
#'   `Ginv` or `Hinv` internally.
#' @param missing_code Integer code used to identify missing genotypes in `X`.
#' @param blend Blending factor applied to `G` before optional tuning.
#' @param chunk_size Number of SNP columns processed per chunk when building
#'   `G` internally.
#' @param n_threads Number of OpenMP threads used in the compiled code.
#' @param tunedG Integer tuning option for `G`. Use `0` for no tuning.
#' @param tau Scaling factor multiplying `G^{-1}` in the construction of
#'   `H^{-1}`.
#' @param omega Scaling factor multiplying `A22^{-1}` in the construction of
#'   `H^{-1}`.
#' @param scale_pevd Logical. If `TRUE`, the returned `PEVD` matrix is divided
#'   by `sigma2a` before being stored in the output object.
#' @param year_col Optional character string naming the year column in `data`.
#'   Required when `year_window` is specified.
#' @param year_window Optional numeric vector of length 2 specifying the time
#'   window to retain, for example `c(2003, 2022)`. If `NULL`, no temporal
#'   filtering is applied.
#' @param min_records_per_year Integer giving the minimum number of records per
#'   MU-year combination for a year to be considered valid when computing
#'   temporal overlap. Only used when `year_window` is not `NULL`.
#' @param verbose Logical. If `TRUE`, progress messages are printed.
#'
#' @return An object of class `"connectedness"`, which is a list with:
#' \describe{
#'   \item{CD}{Numeric matrix (U x U). Pairwise CD contrast values between MUs.}
#'   \item{PEVD}{Numeric matrix (U x U). Pairwise PEVD contrast values between MUs.
#'     If `scale_pevd = TRUE`, this matrix is returned on the scale `PEVD / sigma2a`.}
#'   \item{qK}{Numeric matrix (U x U). Denominator of the contrast under the
#'     kernel used in the analysis.}
#'   \item{qC}{Numeric matrix (U x U). Prediction error numerator of the contrast.}
#'   \item{n_target}{Named numeric vector. Number of target animals per MU.}
#'   \item{relationship}{Character string indicating the inverse relationship
#'     matrix used in the analysis.}
#'   \item{year_window}{The time window used, or `NULL` if no filtering was applied.}
#'   \item{overlap}{A data frame describing temporal overlap between MU pairs,
#'     or `NULL` if no temporal filtering was requested.}
#'   \item{call}{The matched function call.}
#' }
#'
#' @details
#' For each pair of management units \eqn{(i, j)}, the contrast assigns weight
#' \eqn{1/n_i} to animals in MU \eqn{i} and \eqn{-1/n_j} to animals in MU
#' \eqn{j}, where \eqn{n_i} and \eqn{n_j} are the numbers of target animals in
#' each unit.
#'
#' `PEVD` measures the prediction error variance of these contrasts, whereas
#' `CD` rescales the contrast information relative to the denominator implied by
#' the kernel actually used in the analysis. The object component `qK` therefore
#' refers generically to the denominator under `A`, `G`, `H`, or a custom kernel.
#'
#' When a time window is specified, the function also reports the years in which
#' MU pairs overlap according to the observed records and the chosen minimum
#' record threshold.
#'
#' @seealso [renum_pedigree()], [build_Ainv()], [build_Ginv()], [build_Hinv()],
#'   [print.connectedness()], [plot.connectedness()]
#'
#' @examples
#' \dontrun{
#' # Pedigree-based connectedness
#' res_A <- compute_connectedness(
#'   data          = my_data,
#'   animal_col    = "animal_id",
#'   mu_col        = "region",
#'   fixed_formula = ~ 1 + sex + birth_year,
#'   sigma2a       = 5.66,
#'   sigma2e       = 10.24,
#'   relationship  = "Ainv",
#'   pedigree      = my_pedigree
#' )
#'
#' # Genomic connectedness
#' res_G <- compute_connectedness(
#'   data          = my_genotyped_data,
#'   animal_col    = "animal_id",
#'   mu_col        = "region",
#'   fixed_formula = ~ 1 + sex,
#'   sigma2a       = 5.66,
#'   sigma2e       = 10.24,
#'   relationship  = "Ginv",
#'   X             = my_genotypes,
#'   animal_index  = my_index
#' )
#'
#' # Single-step connectedness
#' res_H <- compute_connectedness(
#'   data          = my_data,
#'   animal_col    = "animal_id",
#'   mu_col        = "region",
#'   fixed_formula = ~ 1 + sex,
#'   sigma2a       = 5.66,
#'   sigma2e       = 10.24,
#'   relationship  = "Hinv",
#'   pedigree      = my_pedigree,
#'   X             = my_genotypes,
#'   genotyped_idx = my_genotyped_idx
#' )
#' }
#'
#' @export
compute_connectedness <- function(
    data,
    animal_col,
    mu_col,
    fixed_formula,
    sigma2a,
    sigma2e,
    relationship         = c("Ainv", "Ginv", "Hinv", "custom"),
    pedigree             = NULL,
    X                    = NULL,
    genotyped_idx        = NULL,
    animal_index         = NULL,
    rel_matrix           = NULL,
    maf_threshold        = 0.05,
    missing_code         = 5L,
    blend                = 0.05,
    chunk_size           = 2000L,
    n_threads            = 1L,
    tunedG               = 0L,
    tau                  = 1.0,
    omega                = 1.0,
    scale_pevd           = FALSE,
    year_col             = NULL,
    year_window          = NULL,
    min_records_per_year = 10,
    verbose              = TRUE
) {

  cl <- match.call()
  relationship <- match.arg(relationship)

  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.")
  }
  if (!animal_col %in% names(data)) {
    stop(sprintf("Column '%s' not found in 'data'.", animal_col))
  }
  if (!mu_col %in% names(data)) {
    stop(sprintf("Column '%s' not found in 'data'.", mu_col))
  }
  if (!inherits(fixed_formula, "formula")) {
    stop("'fixed_formula' must be a formula, e.g. ~ 1 + sex + year.")
  }
  if (!is.numeric(sigma2a) || length(sigma2a) != 1L || is.na(sigma2a) || sigma2a <= 0) {
    stop("'sigma2a' must be a single positive number.")
  }
  if (!is.numeric(sigma2e) || length(sigma2e) != 1L || is.na(sigma2e) || sigma2e <= 0) {
    stop("'sigma2e' must be a single positive number.")
  }
  if (!is.logical(scale_pevd) || length(scale_pevd) != 1L || is.na(scale_pevd)) {
    stop("'scale_pevd' must be TRUE or FALSE.")
  }
  if (!is.numeric(min_records_per_year) || length(min_records_per_year) != 1L ||
      is.na(min_records_per_year) || min_records_per_year < 1 ||
      min_records_per_year != as.integer(min_records_per_year)) {
    stop("'min_records_per_year' must be a single positive integer.")
  }

  if (!is.null(year_window)) {
    if (is.null(year_col)) {
      stop("'year_col' is required when 'year_window' is specified.")
    }
    if (!year_col %in% names(data)) {
      stop(sprintf("Column '%s' (year_col) not found in 'data'.", year_col))
    }
    if (!is.numeric(year_window) || length(year_window) != 2L || year_window[1] > year_window[2]) {
      stop("'year_window' must be a numeric vector of length 2 with year_window[1] <= year_window[2].")
    }
  }
  if (is.null(year_window) && !is.null(year_col) && min_records_per_year != 10 && verbose) {
    message("'min_records_per_year' is ignored when 'year_window' is NULL.")
  }

  data[[animal_col]] <- as.character(data[[animal_col]])
  data[[mu_col]]     <- as.character(data[[mu_col]])

  mu_check <- unique(data[, c(animal_col, mu_col)])
  dup_animals <- mu_check[[animal_col]][duplicated(mu_check[[animal_col]])]
  if (length(dup_animals) > 0) {
    stop(sprintf(
      "%d animal(s) appear in more than one MU. Each animal must belong to exactly one MU.\n  First offenders: %s",
      length(dup_animals),
      paste(head(dup_animals, 5), collapse = ", ")
    ))
  }

  renum <- NULL

  if (relationship %in% c("Ainv", "Hinv") ||
      (relationship == "Ginv" && is.null(animal_index) && !is.null(pedigree))) {
    if (is.null(pedigree)) {
      stop(sprintf("'pedigree' is required when relationship = '%s'.", relationship))
    }
    if (verbose) message("Renumbering pedigree...")
    renum <- renum_pedigree(pedigree, verbose = verbose)
  }

  if (relationship == "Ainv") {
    if (verbose) message("Building A-inverse...")
    Ainv_res   <- build_Ainv(renum)
    rel_matrix <- Ainv_res$Ainv
    animal_index <- stats::setNames(renum$new_id, renum$animal)
    animal_index <- animal_index[!is.na(names(animal_index))]

  } else if (relationship == "Ginv") {
    if (is.null(X)) {
      stop("'X' is required when relationship = 'Ginv'.")
    }

    if (is.null(animal_index)) {
      if (is.null(renum) || is.null(genotyped_idx)) {
        stop("For relationship = 'Ginv', provide either 'animal_index' directly or both 'pedigree' and 'genotyped_idx'.")
      }
      if (length(genotyped_idx) != nrow(as.matrix(X))) {
        stop("When deriving 'animal_index' for Ginv, length(genotyped_idx) must equal nrow(X).")
      }
      animal_ids_gen <- renum$animal[match(as.integer(genotyped_idx), renum$new_id)]
      if (any(is.na(animal_ids_gen))) {
        stop("Failed to map some entries of 'genotyped_idx' back to original animal IDs.")
      }
      animal_index <- stats::setNames(seq_along(genotyped_idx), animal_ids_gen)
      storage.mode(animal_index) <- "integer"
    }

    if (verbose) message("Building G-inverse...")
    Ginv_res <- build_Ginv(
      X             = X,
      maf_threshold = maf_threshold,
      missing_code  = missing_code,
      blend         = blend,
      chunk_size    = chunk_size,
      n_threads     = n_threads,
      tunedG        = tunedG,
      A22           = NULL
    )
    rel_matrix <- Matrix::Matrix(Ginv_res$Ginv, sparse = TRUE)

  } else if (relationship == "Hinv") {
    if (is.null(X)) {
      stop("'X' is required when relationship = 'Hinv'.")
    }
    if (is.null(genotyped_idx)) {
      stop("'genotyped_idx' is required when relationship = 'Hinv'.")
    }
    if (length(genotyped_idx) != nrow(as.matrix(X))) {
      stop("length(genotyped_idx) must equal nrow(X) when relationship = 'Hinv'.")
    }

    if (verbose) message("Building H-inverse...")
    Hinv_res <- build_Hinv(
      renum                = renum,
      genotyped_idx        = genotyped_idx,
      X                    = X,
      maf_threshold        = maf_threshold,
      missing_code         = missing_code,
      blend                = blend,
      chunk_size           = chunk_size,
      n_threads            = n_threads,
      tunedG               = tunedG,
      tau                  = tau,
      omega                = omega,
      return_Ainv          = FALSE,
      return_F             = FALSE,
      return_A22           = FALSE,
      return_Ginv          = FALSE,
      return_allele_freqs  = FALSE
    )
    rel_matrix <- Hinv_res$Hinv

    animal_index <- stats::setNames(renum$new_id, renum$animal)
    animal_index <- animal_index[!is.na(names(animal_index))]

  } else if (relationship == "custom") {
    if (is.null(rel_matrix)) {
      stop("'rel_matrix' must be supplied when relationship = 'custom'.")
    }
    if (is.null(animal_index)) {
      stop("'animal_index' must be supplied when relationship = 'custom'.")
    }
  }

  if (!inherits(rel_matrix, "dgCMatrix")) {
    rel_matrix <- methods::as(rel_matrix, "dgCMatrix")
  }
  if (nrow(rel_matrix) != ncol(rel_matrix)) {
    stop("'rel_matrix' must be square.")
  }

  if (is.null(names(animal_index))) {
    stop("'animal_index' must be a named integer vector.")
  }
  animal_names <- names(animal_index)
  animal_index <- as.integer(animal_index)
  names(animal_index) <- animal_names
  if (any(is.na(names(animal_index))) || any(names(animal_index) == "")) {
    stop("'animal_index' names must be non-missing, non-empty animal IDs.")
  }
  if (anyDuplicated(names(animal_index))) {
    stop("'animal_index' names (animal IDs) must be unique.")
  }
  if (any(is.na(animal_index)) || any(animal_index < 1L) || any(animal_index > nrow(rel_matrix))) {
    stop("'animal_index' values must be integers in 1..N, where N = nrow(rel_matrix).")
  }

  N <- nrow(rel_matrix)

  ids_in_data <- data[[animal_col]]
  missing_ids <- setdiff(ids_in_data, names(animal_index))
  if (length(missing_ids) > 0) {
    stop(sprintf(
      "%d animal ID(s) in 'data' have no entry in the relationship matrix used for connectedness.\n  First missing: %s",
      length(missing_ids),
      paste(head(missing_ids, 5), collapse = ", ")
    ))
  }

  data$.new_id <- unname(animal_index[ids_in_data])
  if (any(is.na(data$.new_id))) {
    stop("Failed to map some records in 'data' to row indices in the relationship matrix. Check 'animal_index'.")
  }

  overlap_dt <- NULL

  if (!is.null(year_window)) {
    data[[year_col]] <- as.integer(data[[year_col]])
    Y1 <- year_window[1]
    Y2 <- year_window[2]

    if (verbose) message(sprintf("Filtering records to year window [%d, %d]...", Y1, Y2))
    data_window <- data[data[[year_col]] >= Y1 & data[[year_col]] <= Y2, ]

    if (nrow(data_window) == 0) {
      stop("No records remain after applying 'year_window'. Check the year range.")
    }

    overlap_dt <- .compute_overlap(data_window, mu_col, year_col, min_records_per_year)
  } else {
    data_window <- data
  }

  mu_levels <- sort(unique(data_window[[mu_col]]))
  U <- length(mu_levels)
  if (U < 2) {
    stop("At least 2 MUs with records are required to compute connectedness.")
  }

  mu_map <- stats::setNames(seq_along(mu_levels), mu_levels)

  mu_animal <- integer(N)
  target    <- logical(N)

  animal_mu <- unique(data_window[, c(".new_id", mu_col)])
  mu_animal[animal_mu$.new_id] <- mu_map[animal_mu[[mu_col]]]
  target[unique(data_window$.new_id)] <- TRUE

  tab_target <- table(mu_animal[target])
  if (verbose) {
    message("Target animals per MU:")
    print(tab_target)
  }

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required.")
  }
  Xsp <- Matrix::sparse.model.matrix(fixed_formula, data = data_window)

  id_rec <- as.integer(data_window$.new_id)

  if (verbose) {
    message(sprintf("Computing CD and PEVD via MME using %s...", relationship))
  }
  res_cpp <- .Call(
    `_connectedness_cd_contrast_mu_mme_sparse`,
    rel_matrix,
    id_rec,
    Xsp,
    as.integer(mu_animal),
    target,
    sigma2a,
    sigma2e,
    as.character(mu_levels),
    PACKAGE = "connectedness"
  )

  pevd_out <- res_cpp$PEVD
  if (scale_pevd) {
    pevd_out <- pevd_out / sigma2a
  }

  out <- structure(
    list(
      CD           = res_cpp$CD,
      PEVD         = pevd_out,
      qK           = res_cpp$qK,
      qC           = res_cpp$qC,
      n_target     = res_cpp$n_target_by_MU,
      relationship = relationship,
      year_window  = year_window,
      overlap      = overlap_dt,
      call         = cl
    ),
    class = "connectedness"
  )

  if (verbose) message("Done.")
  out
}

.compute_overlap <- function(data_window, mu_col, year_col, min_records_per_year) {

  df <- as.data.frame(data_window, stringsAsFactors = FALSE)
  df[[mu_col]]   <- as.character(df[[mu_col]])
  df[[year_col]] <- as.integer(df[[year_col]])

  counts <- stats::aggregate(
    rep(1L, nrow(df)),
    by = list(MU = df[[mu_col]], Year = df[[year_col]]),
    FUN = length
  )
  names(counts)[3] <- "N"

  counts <- counts[counts$N >= min_records_per_year, c("MU", "Year"), drop = FALSE]
  if (!nrow(counts)) {
    return(data.frame(MU1 = character(0), MU2 = character(0), Year = integer(0)))
  }

  years_by_mu <- split(counts$Year, counts$MU)
  mus <- sort(names(years_by_mu))

  out_list <- list()
  idx <- 1L

  for (i in seq_len(length(mus) - 1L)) {
    for (j in seq.int(i + 1L, length(mus))) {
      yy <- intersect(sort(unique(years_by_mu[[mus[i]]])), sort(unique(years_by_mu[[mus[j]]])))
      if (length(yy)) {
        out_list[[idx]] <- data.frame(
          MU1 = mus[i],
          MU2 = mus[j],
          Year = as.integer(yy),
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }

  if (!length(out_list)) {
    return(data.frame(MU1 = character(0), MU2 = character(0), Year = integer(0)))
  }

  do.call(rbind, out_list)
}
