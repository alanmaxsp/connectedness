#' Compute connectedness metrics between management units
#'
#' Calculates the Coefficient of Determination (CD) and the Prediction Error
#' Variance of Differences (PEVD) between all pairs of management units (MUs),
#' using the contrast approach via Mixed Model Equations (MME).
#'
#' @param data A data frame containing phenotypic records. Must include columns
#'   for animal ID, MU membership, and all variables in `fixed_formula`.
#' @param animal_col Character. Name of the column in `data` with animal IDs.
#'   These IDs must match those in `pedigree` (if used) or in `animal_index`.
#' @param mu_col Character. Name of the column in `data` with MU labels.
#'   Each animal must belong to exactly one MU.
#' @param fixed_formula A one-sided formula for the fixed effects of the MME,
#'   e.g. `~ 1 + sex + year_of_birth`. The intercept is included by default.
#' @param sigma2a Positive numeric. Additive genetic variance (from a prior
#'   variance components analysis, e.g. via REML or Bayesian methods).
#' @param sigma2e Positive numeric. Residual variance.
#' @param pedigree A data frame with exactly three columns (animal, sire, dam).
#'   Required if `rel_matrix` is `NULL`. Should contain the **complete**
#'   pedigree including ancestors without phenotypic records.
#' @param rel_matrix Optional. A sparse matrix of class `dgCMatrix` (N x N),
#'   the **inverse** of a relationship matrix (A-inverse, G-inverse, H-inverse,
#'   etc.). If provided, `pedigree` is ignored and `animal_index` is required.
#'   See Details for format requirements.
#' @param animal_index Optional named integer vector. Required when `rel_matrix`
#'   is provided. Names are original animal IDs (as in `data[[animal_col]]`);
#'   values are the corresponding row/column indices (1..N) in `rel_matrix`.
#' @param year_col Character or `NULL` (default). Name of the column in `data`
#'   with the year of record. Required if `year_window` is specified.
#' @param year_window Numeric vector of length 2, e.g. `c(2003, 2022)`, or
#'   `NULL` (default). If `NULL`, all records are used without temporal
#'   filtering.
#' @param min_records_per_year Integer. Minimum number of records per
#'   MU-year combination to consider a year as valid within a MU when
#'   computing the temporal overlap. Only used if `year_window` is not `NULL`.
#'   Default is 10.
#' @param verbose Logical. If `TRUE` (default), prints progress messages.
#'
#' @return An object of class `"connectedness"`, which is a list with:
#' \describe{
#'   \item{CD}{Numeric matrix (U x U). Coefficient of Determination for each
#'     pair of MUs. Diagonal is `NA`. Symmetric.}
#'   \item{PEVD}{Numeric matrix (U x U). Prediction Error Variance of
#'     Differences (in units of `sigma2e`). Diagonal is `NA`. Symmetric.}
#'   \item{qA}{Numeric matrix (U x U). Genetic denominator of the contrast
#'     (based on A). Used internally to compute CD.}
#'   \item{qC}{Numeric matrix (U x U). PEV numerator of the contrast
#'     (based on C-inverse). Used internally to compute CD and PEVD.}
#'   \item{n_target}{Named integer vector. Number of target animals per MU.}
#'   \item{year_window}{The year window used (or `NULL` if not applied).}
#'   \item{overlap}{A `data.table` with columns MU1, MU2, Year describing the
#'     temporal overlap between MU pairs. `NULL` if no year filtering was used.}
#'   \item{call}{The matched call.}
#' }
#'
#' @section Using a custom relationship matrix:
#' When passing `rel_matrix` directly (e.g. a genomic G-inverse):
#' \itemize{
#'   \item Must be class `dgCMatrix`, symmetric, and positive definite.
#'   \item Dimensions must be N x N where N equals the number of animals
#'     indexed by `animal_index`.
#'   \item `animal_index` must be a named integer vector mapping every animal
#'     ID that appears in `data[[animal_col]]` to its row/column in `rel_matrix`.
#'   \item `data` should be restricted to the animals present in `rel_matrix`.
#'   \item The user is responsible for ensuring correct scaling (same units as
#'     `sigma2a`).
#' }
#'
#' @details
#' The CD and PEVD are computed following the contrast approach of
#' Laloë (1993) and Laloë et al. (1996), as implemented in the MME framework.
#' For each pair of MUs \eqn{(i, j)}, the contrast vector \eqn{b} assigns
#' weight \eqn{1/n_i} to animals in MU i and \eqn{-1/n_j} to animals in MU j,
#' where \eqn{n_i} and \eqn{n_j} are the number of target animals.
#'
#' \deqn{CD_{ij} = 1 - \lambda \cdot \frac{b'C^{-1}b}{b'Ab}}
#' \deqn{PEVD_{ij} = \sigma^2_e \cdot b'C^{-1}b}
#'
#' where \eqn{\lambda = \sigma^2_e / \sigma^2_a} and C is the coefficient
#' matrix of the MME.
#'
#' @references
#' Laloë, D. (1993). Precision and information in linear models of genetic
#' evaluation. *Genetics Selection Evolution*, 25, 557–576.
#'
#' Laloë, D., Phocas, F., & Ménissier, F. (1996). Considerations on measures
#' of precision and connectedness in mixed linear models of genetic evaluation.
#' *Genetics Selection Evolution*, 28, 359–378.
#'
#' @examples
#' \dontrun{
#' # With pedigree (builds A-inverse internally)
#' res <- compute_connectedness(
#'   data          = my_data,
#'   pedigree      = my_pedigree,
#'   animal_col    = "animal_id",
#'   mu_col        = "region",
#'   fixed_formula = ~ 1 + sex + birth_year,
#'   sigma2a       = 5.66,
#'   sigma2e       = 10.24,
#'   year_col      = "birth_year",
#'   year_window   = c(2003, 2022)
#' )
#' print(res)
#' plot(res)
#'
#' # With a pre-built G-inverse (genomic)
#' res_g <- compute_connectedness(
#'   data          = my_genotyped_data,
#'   animal_col    = "animal_id",
#'   mu_col        = "region",
#'   fixed_formula = ~ 1 + sex,
#'   sigma2a       = 5.66,
#'   sigma2e       = 10.24,
#'   rel_matrix    = my_Ginv,        # dgCMatrix
#'   animal_index  = my_index        # named integer vector
#' )
#' }
#'
#' @seealso [renum_pedigree()], [build_Ainv()], [print.connectedness()],
#'   [plot.connectedness()]
#'
#' @export
compute_connectedness <- function(
    data,
    animal_col,
    mu_col,
    fixed_formula,
    sigma2a,
    sigma2e,
    pedigree             = NULL,
    rel_matrix           = NULL,
    animal_index         = NULL,
    year_col             = NULL,
    year_window          = NULL,
    min_records_per_year = 10,
    verbose              = TRUE
) {

  cl <- match.call()

  # ------------------------------------------------------------------
  # 1. Input validation
  # ------------------------------------------------------------------
  if (!is.data.frame(data))
    stop("'data' must be a data frame.")
  if (!animal_col %in% names(data))
    stop(sprintf("Column '%s' not found in 'data'.", animal_col))
  if (!mu_col %in% names(data))
    stop(sprintf("Column '%s' not found in 'data'.", mu_col))
  if (!inherits(fixed_formula, "formula"))
    stop("'fixed_formula' must be a formula, e.g. ~ 1 + sex + year.")
  if (!is.numeric(sigma2a) || length(sigma2a) != 1 || sigma2a <= 0)
    stop("'sigma2a' must be a single positive number.")
  if (!is.numeric(sigma2e) || length(sigma2e) != 1 || sigma2e <= 0)
    stop("'sigma2e' must be a single positive number.")
  if (!is.numeric(min_records_per_year) || length(min_records_per_year) != 1 ||
      is.na(min_records_per_year) || min_records_per_year < 1 ||
      min_records_per_year != as.integer(min_records_per_year))
    stop("'min_records_per_year' must be a single positive integer.")

  # rel_matrix / pedigree exclusion logic
  if (is.null(rel_matrix) && is.null(pedigree))
    stop("Either 'pedigree' or 'rel_matrix' must be provided.")
  if (!is.null(rel_matrix) && !is.null(pedigree) && verbose)
    message("Both 'rel_matrix' and 'pedigree' provided. Using 'rel_matrix'; pedigree is ignored.")
  if (!is.null(rel_matrix) && is.null(animal_index))
    stop("'animal_index' is required when 'rel_matrix' is provided. ",
         "It must be a named integer vector mapping animal IDs to row/column ",
         "indices in 'rel_matrix'.")
  if (!is.null(rel_matrix)) {
    if (!inherits(rel_matrix, "dgCMatrix"))
      stop("'rel_matrix' must be of class 'dgCMatrix' (sparse matrix from the Matrix package).")
    if (nrow(rel_matrix) != ncol(rel_matrix))
      stop("'rel_matrix' must be square.")
    if (!is.integer(animal_index) || is.null(names(animal_index)))
      stop("'animal_index' must be a named integer vector.")
    if (any(is.na(names(animal_index))) || any(names(animal_index) == ""))
      stop("'animal_index' names must be non-missing, non-empty animal IDs.")
    if (anyDuplicated(names(animal_index)))
      stop("'animal_index' names (animal IDs) must be unique.")
    if (any(is.na(animal_index)) || any(animal_index < 1L) ||
        any(animal_index > nrow(rel_matrix)))
      stop("'animal_index' values must be integers in 1..N, where N = nrow(rel_matrix).")
  }

  # year_window validation
  if (!is.null(year_window)) {
    if (is.null(year_col))
      stop("'year_col' is required when 'year_window' is specified.")
    if (!year_col %in% names(data))
      stop(sprintf("Column '%s' (year_col) not found in 'data'.", year_col))
    if (!is.numeric(year_window) || length(year_window) != 2 || year_window[1] > year_window[2])
      stop("'year_window' must be a numeric vector of length 2 with year_window[1] <= year_window[2].")
  }
  if (is.null(year_window) && !is.null(year_col) && min_records_per_year != 10 && verbose)
    message("'min_records_per_year' is ignored when 'year_window' is NULL.")

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
  # 3. Build relationship matrix if not provided
  # ------------------------------------------------------------------
  if (is.null(rel_matrix)) {
    if (verbose) message("Renumbering pedigree...")
    renum <- renum_pedigree(pedigree, verbose = verbose)

    if (verbose) message("Building A-inverse...")
    Ainv_res    <- build_Ainv(renum)
    rel_matrix  <- Ainv_res$Ainv

    # animal_index: maps original animal ID -> new integer index
    animal_index <- stats::setNames(renum$new_id, renum$animal)
    animal_index <- animal_index[!is.na(names(animal_index))]
  }

  N <- nrow(rel_matrix)

  # ------------------------------------------------------------------
  # 4. Assign integer index to each record in data
  # ------------------------------------------------------------------
  ids_in_data <- data[[animal_col]]
  missing_ids <- setdiff(ids_in_data, names(animal_index))
  if (length(missing_ids) > 0)
    stop(sprintf(
      "%d animal ID(s) in 'data' have no entry in the relationship matrix / pedigree.\n  First missing: %s",
      length(missing_ids),
      paste(head(missing_ids, 5), collapse = ", ")
    ))

  data$.new_id <- animal_index[ids_in_data]
  if (any(is.na(data$.new_id)))
    stop("Failed to map some records in 'data' to row indices in the relationship matrix. ",
         "Check 'animal_index' names and values.")

  # ------------------------------------------------------------------
  # 5. Temporal filtering (optional)
  # ------------------------------------------------------------------
  overlap_dt <- NULL

  if (!is.null(year_window)) {
    data[[year_col]] <- as.integer(data[[year_col]])
    Y1 <- year_window[1]; Y2 <- year_window[2]

    if (verbose) message(sprintf("Filtering records to year window [%d, %d]...", Y1, Y2))
    data_window <- data[data[[year_col]] >= Y1 & data[[year_col]] <= Y2, ]

    if (nrow(data_window) == 0)
      stop("No records remain after applying 'year_window'. Check the year range.")

    # compute temporal overlap table
    overlap_dt <- .compute_overlap(
      data_window, mu_col, year_col, min_records_per_year
    )
  } else {
    data_window <- data
  }

  # ------------------------------------------------------------------
  # 6. Build mu_animal and target vectors
  # ------------------------------------------------------------------
  mu_levels <- sort(unique(data_window[[mu_col]]))
  U <- length(mu_levels)
  if (U < 2)
    stop("At least 2 MUs with records are required to compute connectedness.")

  mu_map <- stats::setNames(seq_along(mu_levels), mu_levels)

  mu_animal <- integer(N)
  target    <- logical(N)

  # one row per animal in the window (already validated one-animal-one-MU)
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
  if (!requireNamespace("Matrix", quietly = TRUE))
    stop("Package 'Matrix' is required.")

  Xsp <- Matrix::sparse.model.matrix(fixed_formula, data = data_window)

  # ------------------------------------------------------------------
  # 8. Call C++ core
  # ------------------------------------------------------------------
  id_rec <- as.integer(data_window$.new_id)

  if (verbose) message("Computing CD and PEVD via MME...")
  res_cpp <- cd_contrast_mu_mme_sparse(
    Ainv             = rel_matrix,
    id_rec           = id_rec,
    X                = Xsp,
    mu_animal        = as.integer(mu_animal),
    target_nullable  = target,
    sigma2a          = sigma2a,
    sigma2e          = sigma2e,
    mu_names_nullable = as.character(mu_levels)
  )

  # ------------------------------------------------------------------
  # 9. Assemble output
  # ------------------------------------------------------------------
  out <- structure(
    list(
      CD          = res_cpp$CD,
      PEVD        = res_cpp$PEVD,
      qA          = res_cpp$qA,
      qC          = res_cpp$qC,
      n_target    = res_cpp$n_target_by_MU,
      year_window = year_window,
      overlap     = overlap_dt,
      call        = cl
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

  # valid MU-year combinations (above threshold)
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
