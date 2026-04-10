#' Build the sparse inverse of the H matrix
#'
#' Computes \eqn{H^{-1}} from a renumbered pedigree and a genotype matrix using
#' a combined pedigree-genomic relationship structure (H-kernel). Internally,
#' the function performs the full pipeline:
#' \enumerate{
#'   \item build \eqn{A^{-1}},
#'   \item compute individual inbreeding coefficients,
#'   \item construct \eqn{A_{22}},
#'   \item build and optionally tune \eqn{G}, then invert it,
#'   \item combine all components to obtain \eqn{H^{-1}}.
#' }
#'
#' @param renum A data frame as returned by [renum_pedigree()], with columns
#'   `new_id`, `new_sire`, and `new_dam`.
#' @param genotyped_idx Integer vector with the 1-based renumbered pedigree
#'   indices of the genotyped animals, in the same order as the rows of `X`.
#' @param X Genotype matrix of dimension `n_gen x m`, coded as 0/1/2.
#'   Missing genotypes must be coded with `missing_code`.
#' @param maf_threshold Minor allele frequency threshold used to filter SNPs.
#' @param missing_code Integer code representing missing genotypes.
#' @param blend Blending factor applied to \eqn{G} before optional tuning.
#' @param chunk_size Number of SNP columns processed per chunk.
#' @param n_threads Number of OpenMP threads.
#' @param tunedG Integer tuning option for \eqn{G}: `0` = no tuning; `1` = standardize by matching mean diagonal/off-diagonal contrast within \eqn{G}; `2` = affine tuning to match \eqn{A_{22}} mean diagonal and off-diagonal; `3` = shift \eqn{G} by a constant so its global mean matches \eqn{A_{22}}.
#' @param tau Scaling factor multiplying \eqn{G^{-1}} in the final expression.
#' @param omega Scaling factor multiplying \eqn{A_{22}^{-1}} in the final expression.
#' @param return_Ainv Logical; return `Ainv` in the output list.
#' @param return_F Logical; return inbreeding coefficients `F` in the output list.
#' @param return_A22 Logical; return `A22` in the output list.
#' @param return_Ginv Logical; return `Ginv` in the output list.
#' @param return_allele_freqs Logical; return retained SNP allele frequencies.
#'
#' @return A list containing `Hinv` and, optionally, intermediate matrices and
#' diagnostics.
#'
#' @seealso [renum_pedigree()], [build_Ainv()], [build_Ginv()]
#'
#' @export
build_Hinv <- function(renum,
                       genotyped_idx,
                       X,
                       maf_threshold       = 0.05,
                       missing_code        = 5L,
                       blend               = 0.05,
                       chunk_size          = 2000L,
                       n_threads           = 1L,
                       tunedG              = 0L,
                       tau                 = 1.0,
                       omega               = 1.0,
                       return_Ainv         = FALSE,
                       return_F            = FALSE,
                       return_A22          = FALSE,
                       return_Ginv         = FALSE,
                       return_allele_freqs = FALSE) {

  required_cols <- c("new_id", "new_sire", "new_dam")
  if (!all(required_cols %in% names(renum))) {
    stop("'renum' must have columns: new_id, new_sire, new_dam. ",
         "Use renum_pedigree() to generate it.")
  }

  renum <- renum[order(renum$new_id), ]

  sire_vec <- as.integer(renum$new_sire)
  dam_vec  <- as.integer(renum$new_dam)
  if (any(is.na(sire_vec)) || any(is.na(dam_vec))) {
    stop("NA values found in new_sire or new_dam. Check pedigree renumbering.")
  }

  genotyped_idx <- as.integer(genotyped_idx)
  if (length(genotyped_idx) < 1L) {
    stop("'genotyped_idx' must contain at least one genotyped animal.")
  }
  if (any(is.na(genotyped_idx))) {
    stop("'genotyped_idx' contains NA values.")
  }
  if (anyDuplicated(genotyped_idx)) {
    stop("'genotyped_idx' must not contain duplicated pedigree indices.")
  }

  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  if (!is.numeric(X) && !is.integer(X)) {
    stop("'X' must be a numeric or integer matrix with genotypes coded as 0/1/2.")
  }
  if (length(dim(X)) != 2L) {
    stop("'X' must be a two-dimensional matrix.")
  }
  if (nrow(X) != length(genotyped_idx)) {
    stop("Number of rows in 'X' must equal length(genotyped_idx).")
  }
  if (ncol(X) < 2L) {
    stop("'X' must have at least 2 SNP columns.")
  }
  if (anyNA(X)) {
    stop("'X' must not contain NA values. Use 'missing_code' for missing genotypes.")
  }
  storage.mode(X) <- "integer"

  .Call(
    `_connectedness_compute_Hinv_from_X`,
    sire_vec,
    dam_vec,
    genotyped_idx,
    X,
    maf_threshold,
    as.integer(missing_code),
    blend,
    as.integer(chunk_size),
    as.integer(n_threads),
    as.integer(tunedG),
    tau,
    omega,
    return_Ainv,
    return_F,
    return_A22,
    return_Ginv,
    return_allele_freqs,
    PACKAGE = "connectedness"
  )
}
