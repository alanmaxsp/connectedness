#' Compute connectedness metrics between management units
#'
#' Calculates the Coefficient of Determination (CD) and the Prediction Error
#' Variance of Differences (PEVD) between all pairs of management units (MUs)
#' using the contrast approach via Mixed Model Equations (MME).
#'
#' The user can choose whether connectedness is computed from a pedigree-based
#' inverse relationship matrix (A-inverse), a genomic inverse relationship
#' matrix (G-inverse), a single-step inverse relationship matrix (H-inverse),
#' or a custom inverse kernel supplied directly.
#'
#' @param data A data frame containing phenotypic records. Must include columns
#'   for animal ID, MU membership, and all variables in `fixed_formula`.
#' @param animal_col Character. Name of the column in `data` with animal IDs.
#' @param mu_col Character. Name of the column in `data` with MU labels.
#'   Each animal must belong to exactly one MU.
#' @param fixed_formula A one-sided formula for the fixed effects of the MME,
#'   e.g. `~ 1 + sex + birth_year`. The intercept is included by default.
#' @param sigma2a Positive numeric. Additive genetic variance (or, more
#'   generally, the variance associated with the animal effect).
#' @param sigma2e Positive numeric. Residual variance.
#' @param relationship Character string indicating which inverse relationship
#'   matrix to use. One of `"Ainv"`, `"Ginv"`, `"Hinv"`, or `"custom"`.
#' @param pedigree Optional pedigree data frame with exactly three columns
#'   (animal, sire, dam). Required for `relationship = "Ainv"` and `"Hinv"`.
#'   Also useful for `relationship = "Ginv"` if `animal_index` must be derived
#'   from `genotyped_idx`.
#' @param X Optional genotype matrix (`n_gen x m`, coded 0/1/2). Required for
#'   `relationship = "Ginv"` and `"Hinv"`.
#' @param genotyped_idx Optional integer vector of 1-based renumbered pedigree
#'   indices for the genotyped animals, in the same order as the rows of `X`.
#'   Required for `relationship = "Hinv"`, and also for `relationship = "Ginv"`
#'   when `animal_index` is not supplied explicitly.
#' @param animal_index Optional named integer vector mapping original animal IDs
#'   to row/column indices in the inverse relationship matrix actually used for
#'   connectedness. Required for `relationship = "custom"`; optional otherwise
#'   when it can be derived internally.
#' @param rel_matrix Optional user-supplied inverse relationship matrix
#'   (`dgCMatrix`). Used only when `relationship = "custom"`.
#' @param maf_threshold Minor allele frequency threshold used when building
#'   `Ginv` or `Hinv`.
#' @param missing_code Integer value indicating missing genotypes in `X`.
#' @param blend Blending factor applied to `G` before optional tuning.
#' @param chunk_size Number of SNP columns processed per chunk in G construction.
#' @param n_threads Number of OpenMP threads.
#' @param tunedG Integer tuning option for `G`. Use 0 for no tuning.
#' @param tau Scaling factor multiplying `G^{-1}` in `H^{-1}`.
#' @param omega Scaling factor multiplying `A22^{-1}` in `H^{-1}`.
#' @param year_col Character or `NULL` (default). Name of the column in `data`
#'   with the year of record. Required if `year_window` is specified.
#' @param year_window Numeric vector of length 2, e.g. `c(2003, 2022)`, or
#'   `NULL` (default). If `NULL`, all records are used without temporal
#'   filtering.
#' @param min_records_per_year Integer. Minimum number of records per MU-year
#'   combination to consider a year as valid within a MU when computing temporal
#'   overlap. Only used if `year_window` is not `NULL`.
#' @param verbose Logical. If `TRUE` (default), prints progress messages.
#'
#' @return An object of class `"connectedness"`, which is a list with:
#' \describe{
#'   \item{CD}{Numeric matrix (U x U). Coefficient of Determination for each
#'     pair of MUs. Diagonal is `NA`. Symmetric.}
#'   \item{PEVD}{Numeric matrix (U x U). Prediction Error Variance of
#'     Differences (in units of `sigma2e`). Diagonal is `NA`. Symmetric.}
#'   \item{qK}{Numeric matrix (U x U). Genetic denominator of the contrast
#'     based on the relationship matrix actually used (A, G, H, or custom K).}
#'   \item{qC}{Numeric matrix (U x U). PEV numerator of the contrast, used
#'     internally to compute CD and PEVD.}
#'   \item{n_target}{Named integer vector. Number of target animals per MU.}
#'   \item{relationship}{Character string indicating the inverse relationship
#'     matrix used.}
#'   \item{year_window}{The year window used (or `NULL` if not applied).}
#'   \item{overlap}{A `data.table` with columns MU1, MU2, Year describing the
#'     temporal overlap between MU pairs. `NULL` if no year filtering was used.}
#'   \item{call}{The matched call.}
#' }
#'
#' @details
#' For each pair of MUs \eqn{(i, j)}, the contrast vector assigns weight
#' \eqn{1/n_i} to animals in MU \eqn{i} and \eqn{-1/n_j} to animals in MU
#' \eqn{j}, where \eqn{n_i} and \eqn{n_j} are the numbers of target animals.
#'
#' The metric `qK` always refers to the kernel actually used in the analysis.
#' Therefore, when the analysis is based on `Ainv`, `Ginv`, or `Hinv`, the
#' denominator is computed with the corresponding `A`, `G`, or `H` block.
#'
#' @seealso [renum_pedigree()], [build_Ainv()], [build_Ginv()], [build_Hinv()],
#'   [print.connectedness()], [plot.connectedness()]
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
    year_col             = NULL,
    year_window          = NULL,
    min_records_per_year = 10,
    verbose              = TRUE
) {

  cl <- match.call()
  relationship <- match.arg(relationship)

  # ------------------------------------------------------------------
  # 1. Input validation
  # ------------------------------------------------------------------
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
  if (!is.numeric(min_records_per_year) || length(min_records_per_year) != 1L ||
      is.na(min_records_per_year) || min_records_per_year < 1 ||
      min_records_per_year != as.integer(min_records_per_year)) {
    stop("'min_records_per_year' must be a single positive integer.")
  }

  # year_window validation
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

  # ------------------------------------------------------------------
  # 2. Check one animal - one MU
  # ------------------------------------------------------------------
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

  # ------------------------------------------------------------------
  # 3. Build or validate the inverse relationship matrix
  # ------------------------------------------------------------------
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

  # Validate rel_matrix / animal_index
  if (!inherits(rel_matrix, "dgCMatrix")) {
    rel_matrix <- methods::as(rel_matrix, "dgCMatrix")
  }
  if (nrow(rel_matrix) != ncol(rel_matrix)) {
    stop("'rel_matrix' must be square.")
  }

  if (is.null(names(animal_index))) {
    stop("'animal_index' must be a named integer vector.")
  }
  animal_index <- as.integer(animal_index)
  names(animal_index) <- names(animal_index)
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

  # ------------------------------------------------------------------
  # 4. Assign integer index to each record in data
  # ------------------------------------------------------------------
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

  # ------------------------------------------------------------------
  # 5. Temporal filtering (optional)
  # ------------------------------------------------------------------
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

  # ------------------------------------------------------------------
  # 6. Build mu_animal and target vectors
  # ------------------------------------------------------------------
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

  # ------------------------------------------------------------------
  # 7. Build sparse model matrix for fixed effects
  # ------------------------------------------------------------------
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required.")
  }
  Xsp <- Matrix::sparse.model.matrix(fixed_formula, data = data_window)

  # ------------------------------------------------------------------
  # 8. Call C++ core
  # ------------------------------------------------------------------
  id_rec <- as.integer(data_window$.new_id)

  if (verbose) {
    message(sprintf("Computing CD and PEVD via MME using %s...", relationship))
  }
  res_cpp <- cd_contrast_mu_mme_sparse(
    Kinv              = rel_matrix,
    id_rec            = id_rec,
    X                 = Xsp,
    mu_animal         = as.integer(mu_animal),
    target_nullable   = target,
    sigma2a           = sigma2a,
    sigma2e           = sigma2e,
    mu_names_nullable = as.character(mu_levels)
  )

  # ------------------------------------------------------------------
  # 9. Assemble output
  # ------------------------------------------------------------------
  out <- structure(
    list(
      CD           = res_cpp$CD,
      PEVD         = res_cpp$PEVD,
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


# ----------------------------------------------------------------------
# Internal helper: compute temporal overlap table
# ----------------------------------------------------------------------
.compute_overlap <- function(data_window, mu_col, year_col, min_records_per_year) {

  dt <- data.table::as.data.table(data_window)
  data.table::setnames(dt, c(mu_col, year_col), c(".mu", ".year"))

  tab <- dt[, .N, by = .(.mu, .year)][N >= min_records_per_year]

  yrs <- tab[, .(Years = list(sort(unique(.year)))), by = .mu]
  mus <- sort(unique(yrs$.mu))

  pairs <- data.table::CJ(MU1 = mus, MU2 = mus)
  pairs[, W := Map(
    intersect,
    yrs[match(MU1, .mu), Years],
    yrs[match(MU2, .mu), Years]
  )]
  pairs[, nW := lengths(W)]

  overlap <- pairs[MU1 < MU2 & nW > 0,
                   .(Year = as.integer(unlist(W))),
                   by = .(MU1, MU2)]
  overlap
}
