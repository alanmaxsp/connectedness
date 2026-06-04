#' Compute connectedness metrics between management units
#'
#' Computes connectedness between management units (MUs) under the contrast
#' approach via Mixed Model Equations (MME). The function returns two
#' complementary metrics for all MU pairs: the Coefficient of Determination
#' (`CD`) and the Prediction Error Variance of Differences (`PEVD`).
#'
#' The connectedness analysis can be based on a pedigree-derived inverse
#' relationship matrix (`Ainv`), a genomic inverse relationship matrix (`Ginv`),
#' a combined pedigree-genomic inverse relationship matrix (`Hinv`), or a user-supplied
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
#' @param rel_matrix Optional user-supplied inverse relationship matrix. Sparse
#'   Matrix objects are kept sparse; dense matrices are kept dense for Schur
#'   solver selection when `relationship = "custom"`.
#' @param maf_threshold Minor allele frequency threshold used when building
#'   `Ginv` or `Hinv` internally.
#' @param missing_code Integer code used to identify missing genotypes in `X`.
#' @param blend Blending factor applied to `G` before optional tuning.
#' @param chunk_size Number of SNP columns processed per chunk when building
#'   `G` internally.
#' @param n_threads Number of OpenMP threads used in the genomic computations
#'   of the compiled backend. Currently, this argument affects `Ginv` and the
#'   genomic substeps of `Hinv`, but it is ignored for purely pedigree-based
#'   `Ainv` analyses and is not propagated to the final CD/PEVD MME solver.
#' @param tunedG Integer tuning option for `G`: `0` = no tuning; `1` = standardize by matching mean diagonal/off-diagonal contrast within `G`; `2` = affine tuning to match `A22` mean diagonal and off-diagonal; `3` = shift `G` by a constant so its global mean matches `A22`.
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
#' @param dry_run Logical. If `TRUE`, return problem-size diagnostics before the
#'   final MME solve instead of computing connectedness metrics.
#' @param max_mme_dim Positive integer or `Inf`. Maximum allowed dimension of
#'   the full MME system (`nrow(rel_matrix) + ncol(fixed_effects)`) before the
#'   final sparse direct solve is attempted. Set to `Inf` to disable this safety
#'   check. This safety check is applied only when `mme_backend = "full_mme"`.
#' @param mme_backend Character string indicating the numerical backend used for
#'   the MME solve. `"full_mme"` builds and factorizes the full MME system.
#'   `"schur"` factorizes the animal block and absorbs fixed effects through a
#'   Schur complement. The two backends are algebraically equivalent up to
#'   floating-point roundoff.
#' @param schur_solver Character string indicating the numerical solver used
#'   when `mme_backend = "schur"`. Options are `"auto"`, `"cholmod"`,
#'   `"dense"`, and `"eigen_sparse"`. `"cholmod"` uses CHOLMOD through the
#'   Matrix package for sparse kernels; `"dense"` uses a dense compiled solver;
#'   `"eigen_sparse"` keeps the Eigen sparse solver for diagnostics and small
#'   comparisons.
#' @param schur_block_size Positive integer. Number of MU right-hand-side
#'   columns solved per block in the Schur backend.
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
#'   \item{mme_backend}{Character string indicating the numerical backend used
#'     for the MME solve.}
#'   \item{schur_solver}{Character string indicating the selected Schur solver,
#'     or `NA` when `mme_backend = "full_mme"`.}
#'   \item{year_window}{The time window used, or `NULL` if no filtering was applied.}
#'   \item{overlap}{A data frame describing temporal overlap between MU pairs,
#'     or `NULL` if no temporal filtering was requested.}
#'   \item{call}{The matched function call.}
#' }
#' If `dry_run = TRUE`, the function returns a diagnostics list instead of a
#' connectedness result.
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
#' The `"full_mme"` backend directly factorizes the full MME matrix
#' `[X'X X'Z; Z'X Z'Z + lambda Kinv]`. The `"schur"` backend avoids this full
#' factorization by factorizing `Cuu = Z'Z + lambda Kinv` and using the fixed-
#' effect Schur complement `S = X'X - X'Z Cuu^{-1} Z'X`.
#'
#' The Schur formulation is algebraically independent of the relationship matrix
#' type. With `schur_solver = "auto"`, sparse inverse kernels such as `Ainv` are
#' solved through CHOLMOD via the Matrix package, while dense kernels such as
#' `Ginv` are solved through a dense compiled backend. The `"eigen_sparse"`
#' solver is retained mainly for diagnostics and small-scale comparisons.
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
    dry_run              = FALSE,
    max_mme_dim          = 500000L,
    mme_backend          = c("full_mme", "schur"),
    schur_solver         = c("auto", "cholmod", "dense", "eigen_sparse"),
    schur_block_size     = 16L,
    verbose              = TRUE
) {

  cl <- match.call()
  relationship <- match.arg(relationship)
  mme_backend <- match.arg(mme_backend)
  schur_solver <- match.arg(schur_solver)

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
  if (!is.logical(dry_run) || length(dry_run) != 1L || is.na(dry_run)) {
    stop("'dry_run' must be TRUE or FALSE.")
  }
  if ((!is.numeric(max_mme_dim) && !is.integer(max_mme_dim)) || length(max_mme_dim) != 1L ||
      is.na(max_mme_dim) || max_mme_dim <= 0) {
    stop("'max_mme_dim' must be a single positive number or Inf.")
  }
  if ((!is.numeric(schur_block_size) && !is.integer(schur_block_size)) ||
      length(schur_block_size) != 1L || is.na(schur_block_size) ||
      schur_block_size < 1 || schur_block_size != as.integer(schur_block_size)) {
    stop("'schur_block_size' must be a single positive integer.")
  }
  schur_block_size <- as.integer(schur_block_size)

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
      A22           = NULL,
	  verbose       = verbose
    )
    rel_matrix <- Ginv_res$Ginv

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
      return_allele_freqs  = FALSE,
      verbose              = verbose
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

  rel_matrix <- .normalize_rel_matrix(rel_matrix, relationship)
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

  selected_schur_solver <- .choose_schur_solver(
    rel_matrix = rel_matrix,
    relationship = relationship,
    schur_solver = schur_solver
  )

  diagnostics <- .connectedness_diagnostics(
    relationship = relationship,
    rel_matrix = rel_matrix,
    fixed_matrix = Xsp,
    data_window = data_window,
    mu_levels = mu_levels,
    target = target,
    year_window = year_window,
    scale_pevd = scale_pevd,
    mme_backend = mme_backend,
    schur_solver = schur_solver,
    selected_schur_solver = selected_schur_solver,
    max_mme_dim = max_mme_dim,
    call = cl
  )

  if (dry_run) {
    if (verbose) .print_connectedness_diagnostics(diagnostics)
    return(invisible(diagnostics))
  }

  if (mme_backend == "full_mme") {
    .check_connectedness_problem_size(diagnostics, max_mme_dim)
  }

  id_rec <- as.integer(data_window$.new_id)

  if (mme_backend == "full_mme") {
    if (verbose) {
      message(sprintf(
        "Computing CD and PEVD via full MME using %s relationship...",
        relationship
      ))
    }
    rel_matrix_sparse <- .as_dgCMatrix(rel_matrix)
    res_cpp <- .Call(
      `_connectedness_cd_contrast_mu_mme_sparse`,
      rel_matrix_sparse,
      id_rec,
      Xsp,
      as.integer(mu_animal),
      target,
      sigma2a,
      sigma2e,
      as.character(mu_levels),
      verbose,
      PACKAGE = "connectedness"
    )
  } else if (mme_backend == "schur") {
    if (verbose) {
      message(sprintf(
        "Computing CD and PEVD via Schur backend using %s solver...",
        selected_schur_solver
      ))
    }

    if (selected_schur_solver == "cholmod") {
      res_cpp <- .cd_contrast_mu_schur_cholmod_R(
        Kinv = rel_matrix,
        id_rec = id_rec,
        Xsp = Xsp,
        mu_animal = as.integer(mu_animal),
        target = target,
        sigma2a = sigma2a,
        sigma2e = sigma2e,
        mu_names = as.character(mu_levels),
        block_size = schur_block_size,
        verbose = verbose
      )
    } else if (selected_schur_solver == "dense") {
      dense_gb <- diagnostics$dense_matrix_gb
      if (is.finite(dense_gb) && dense_gb > 32 && verbose) {
        warning(sprintf(
          "Dense Schur solver selected. Dense Kinv storage alone is approximately %.1f GB.",
          dense_gb
        ), call. = FALSE)
      }
      Kinv_dense <- if (.is_sparse_matrix(rel_matrix)) as.matrix(rel_matrix) else as.matrix(rel_matrix)
      storage.mode(Kinv_dense) <- "double"
      res_cpp <- .Call(
        `_connectedness_cd_contrast_mu_schur_dense`,
        Kinv_dense,
        id_rec,
        Xsp,
        as.integer(mu_animal),
        target,
        sigma2a,
        sigma2e,
        as.character(mu_levels),
        schur_block_size,
        verbose,
        PACKAGE = "connectedness"
      )
    } else if (selected_schur_solver == "eigen_sparse") {
      if (!.is_sparse_matrix(rel_matrix)) {
        stop("schur_solver = 'eigen_sparse' requires a sparse relationship matrix.")
      }
      res_cpp <- .Call(
        `_connectedness_cd_contrast_mu_mme_schur_sparse`,
        .as_dgCMatrix(rel_matrix),
        id_rec,
        Xsp,
        as.integer(mu_animal),
        target,
        sigma2a,
        sigma2e,
        as.character(mu_levels),
        verbose,
        PACKAGE = "connectedness"
      )
    } else {
      stop("Unknown Schur solver.")
    }
  } else {
    stop("Unsupported 'mme_backend'.")
  }

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
      relationship  = relationship,
      mme_backend   = mme_backend,
      schur_solver  = if (mme_backend == "schur") selected_schur_solver else NA_character_,
      year_window   = year_window,
      overlap      = overlap_dt,
      call         = cl
    ),
    class = "connectedness"
  )

  if (verbose) message("Done.")
  invisible(out)
}



.is_sparse_matrix <- function(x) {
  inherits(x, "sparseMatrix")
}

.as_dgCMatrix <- function(x) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required.")
  }
  if (!.is_sparse_matrix(x)) {
    x <- Matrix::Matrix(x, sparse = TRUE)
  }
  if (!inherits(x, "dgCMatrix")) {
    x <- methods::as(methods::as(x, "generalMatrix"), "dgCMatrix")
  }
  x
}

.normalize_rel_matrix <- function(rel_matrix, relationship) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required.")
  }

  if (relationship == "Ginv") {
    if (.is_sparse_matrix(rel_matrix)) return(.as_dgCMatrix(rel_matrix))
    return(as.matrix(rel_matrix))
  }

  if (relationship == "Ainv") {
    return(.as_dgCMatrix(rel_matrix))
  }

  if (relationship == "Hinv") {
    if (.is_sparse_matrix(rel_matrix)) return(.as_dgCMatrix(rel_matrix))
    return(as.matrix(rel_matrix))
  }

  if (relationship == "custom") {
    if (.is_sparse_matrix(rel_matrix)) return(.as_dgCMatrix(rel_matrix))
    return(as.matrix(rel_matrix))
  }

  stop("Unsupported relationship.")
}

.choose_schur_solver <- function(rel_matrix,
                                 relationship,
                                 schur_solver = "auto",
                                 dense_density_threshold = 0.20) {
  if (schur_solver != "auto") return(schur_solver)

  is_sparse <- .is_sparse_matrix(rel_matrix)
  N <- nrow(rel_matrix)

  if (!is_sparse) return("dense")
  if (relationship == "Ainv") return("cholmod")
  if (relationship == "Ginv") return("dense")

  nnz <- Matrix::nnzero(rel_matrix)
  density <- nnz / (as.numeric(N) * as.numeric(N))

  if (relationship %in% c("Hinv", "custom")) {
    if (density <= dense_density_threshold) return("cholmod")
    return("dense")
  }

  "dense"
}

.connectedness_diagnostics <- function(relationship,
                                       rel_matrix,
                                       fixed_matrix,
                                       data_window,
                                       mu_levels,
                                       target,
                                       year_window,
                                       scale_pevd,
                                       mme_backend,
                                       schur_solver,
                                       selected_schur_solver,
                                       max_mme_dim,
                                       call) {

  rel_n <- nrow(rel_matrix)
  matrix_storage <- if (.is_sparse_matrix(rel_matrix)) "sparse" else "dense"
  rel_nnz <- if (.is_sparse_matrix(rel_matrix)) {
    Matrix::nnzero(rel_matrix)
  } else {
    as.numeric(nrow(rel_matrix)) * as.numeric(ncol(rel_matrix))
  }
  x_nnz <- Matrix::nnzero(fixed_matrix)
  p <- ncol(fixed_matrix)
  mme_dim <- rel_n + p
  n_records <- nrow(data_window)
  n_target <- sum(target)
  n_mu <- length(mu_levels)
  rel_density <- rel_nnz / (as.numeric(rel_n) * as.numeric(rel_n))

  # This is only the explicit sparse MME storage footprint before symbolic
  # factorization. Sparse direct factorization can require substantially more
  # memory depending on fill-in, so this is a lower-bound diagnostic.
  explicit_mme_nnz_lower_bound <- rel_nnz + 2 * x_nnz + p
  explicit_mme_storage_mb_lower_bound <-
    explicit_mme_nnz_lower_bound * (8 + 4) / 1024^2
  dense_matrix_gb <- as.numeric(rel_n) * as.numeric(rel_n) * 8 / 1024^3
  schur_W_storage_mb <- as.numeric(rel_n) * as.numeric(p) * 8 / 1024^2
  schur_S_storage_mb <- as.numeric(p) * as.numeric(p) * 8 / 1024^2
  schur_W_storage_gb <- schur_W_storage_mb / 1024
  schur_S_storage_gb <- schur_S_storage_mb / 1024
  recommended_backend <- if (is.finite(max_mme_dim) && mme_dim > max_mme_dim) {
    "schur"
  } else {
    "full_mme_or_schur"
  }

  structure(
    list(
      relationship = relationship,
      matrix_storage = matrix_storage,
      n_relationship = rel_n,
      n_records = n_records,
      n_target = n_target,
      n_management_units = n_mu,
      n_fixed_effect_columns = p,
      mme_dim = mme_dim,
      relationship_nonzeros = rel_nnz,
      relationship_density = rel_density,
      dense_matrix_gb = dense_matrix_gb,
      fixed_effect_nonzeros = x_nnz,
      explicit_mme_nonzeros_lower_bound = explicit_mme_nnz_lower_bound,
      explicit_mme_storage_mb_lower_bound = explicit_mme_storage_mb_lower_bound,
      schur_W_storage_mb = schur_W_storage_mb,
      schur_S_storage_mb = schur_S_storage_mb,
      schur_W_storage_gb = schur_W_storage_gb,
      schur_S_storage_gb = schur_S_storage_gb,
      recommended_backend = recommended_backend,
      requested_backend = mme_backend,
      requested_schur_solver = schur_solver,
      selected_schur_solver = selected_schur_solver,
      schur_solver_recommended = .choose_schur_solver(
        rel_matrix = rel_matrix,
        relationship = relationship,
        schur_solver = "auto"
      ),
      year_window = year_window,
      scale_pevd = scale_pevd,
      call = call
    ),
    class = "connectedness_diagnostics"
  )
}

.check_connectedness_problem_size <- function(diagnostics, max_mme_dim) {

  if (is.finite(max_mme_dim) && diagnostics$mme_dim > max_mme_dim) {
    stop(sprintf(
      paste0(
        "The connectedness MME system is too large for the default safe solve limit: ",
        "dimension p + N = %d exceeds max_mme_dim = %d. ",
        "Run compute_connectedness(..., dry_run = TRUE) to inspect problem-size diagnostics. ",
        "Then reduce the analysis scope, provide a smaller user-defined relationship matrix, ",
        "or set max_mme_dim = Inf to attempt the solve at your own risk."
      ),
      diagnostics$mme_dim,
      as.integer(max_mme_dim)
    ))
  }

  invisible(TRUE)
}

.print_connectedness_diagnostics <- function(x) {
  message("Connectedness dry-run diagnostics:")
  message(sprintf("  relationship              : %s", x$relationship))
  message(sprintf("  matrix storage            : %s", x$matrix_storage))
  message(sprintf("  relationship dimension    : %d", x$n_relationship))
  message(sprintf("  records                   : %d", x$n_records))
  message(sprintf("  target animals            : %d", x$n_target))
  message(sprintf("  management units          : %d", x$n_management_units))
  message(sprintf("  fixed-effect columns      : %d", x$n_fixed_effect_columns))
  message(sprintf("  MME dimension (p + N)     : %d", x$mme_dim))
  message(sprintf("  relationship nonzeros     : %d", x$relationship_nonzeros))
  message(sprintf("  relationship density      : %.6g", x$relationship_density))
  message(sprintf("  dense equivalent size     : %.2f GB", x$dense_matrix_gb))
  message(sprintf("  fixed-effect nonzeros     : %d", x$fixed_effect_nonzeros))
  message(sprintf(
    "  explicit MME storage lower bound: %.1f MB",
    x$explicit_mme_storage_mb_lower_bound
  ))
  message(sprintf("  Schur W storage           : %.1f MB (%.2f GB)", x$schur_W_storage_mb, x$schur_W_storage_gb))
  message(sprintf("  Schur S storage           : %.1f MB (%.2f GB)", x$schur_S_storage_mb, x$schur_S_storage_gb))
  message(sprintf("  requested backend         : %s", x$requested_backend))
  message(sprintf("  recommended backend       : %s", x$recommended_backend))
  message(sprintf("  requested Schur solver    : %s", x$requested_schur_solver))
  message(sprintf("  selected Schur solver     : %s", x$selected_schur_solver))
  message(sprintf("  recommended Schur solver  : %s", x$schur_solver_recommended))
  invisible(x)
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
