#' Build the dense inverse of the genomic relationship matrix (G-inverse)
#'
#' Computes a dense \eqn{G^{-1}} for the genotyped animals using a VanRaden
#' method-1 style construction of \eqn{G}, followed by matrix inversion.
#' Optional affine tuning of \eqn{G} can be applied before inversion.
#'
#' @param X Genotype matrix of dimension `n_gen x m`, coded as 0/1/2.
#'   Missing genotypes must be coded with `missing_code`.
#' @param maf_threshold Minor allele frequency threshold used to filter SNPs.
#' @param missing_code Integer code representing missing genotypes.
#' @param blend Blending factor applied to \eqn{G} before optional tuning.
#' @param chunk_size Number of SNP columns processed per chunk.
#' @param n_threads Number of OpenMP threads.
#' @param tunedG Integer tuning option for \eqn{G}. Use 0 for no tuning.
#' @param A22 Optional dense pedigree-based relationship matrix for the same
#'   genotyped animals and in the same order as the rows of `X`. Required when
#'   `tunedG` is 2 or 3.
#'
#' @return A list containing `Ginv` and several diagnostics describing SNP
#'   filtering and tuning.
#'
#' @details
#' The matrix `X` is coerced to an integer matrix before being passed to the
#' compiled backend. When `tunedG` equals 2 or 3, `A22` must be supplied and
#' must have dimensions `n_gen x n_gen`.
#'
#' @seealso [build_Ainv()], [build_Hinv()]
#'
#' @export
build_Ginv <- function(X,
                       maf_threshold = 0.05,
                       missing_code  = 5L,
                       blend         = 0.05,
                       chunk_size    = 2000L,
                       n_threads     = 1L,
                       tunedG        = 0L,
                       A22           = NULL) {

  X <- as.matrix(X)
  if (!is.numeric(X) && !is.integer(X)) {
    stop("'X' must be a numeric or integer matrix with genotypes coded as 0/1/2.")
  }
  if (length(dim(X)) != 2L) {
    stop("'X' must be a two-dimensional matrix.")
  }
  if (nrow(X) < 2L || ncol(X) < 2L) {
    stop("'X' must have at least 2 rows and 2 columns.")
  }
  if (any(is.na(X))) {
    stop("'X' must not contain NA values. Use 'missing_code' for missing genotypes.")
  }
  if (any(X != round(X))) {
    stop("All entries of 'X' must be integer-coded genotypes (e.g. 0/1/2 or missing_code).")
  }
  storage.mode(X) <- "integer"

  if (!is.null(A22)) {
    A22 <- as.matrix(A22)
    if (!is.numeric(A22)) {
      stop("'A22' must be a numeric matrix.")
    }
    if (!all(dim(A22) == c(nrow(X), nrow(X)))) {
      stop("'A22' must be a square matrix with dimensions nrow(X) x nrow(X).")
    }
  }

  if (tunedG %in% c(2L, 3L) && is.null(A22)) {
    stop("'A22' must be supplied when 'tunedG' is 2 or 3.")
  }

  compute_Ginv(
    X             = X,
    maf_threshold = maf_threshold,
    missing_code  = as.integer(missing_code),
    blend         = blend,
    chunk_size    = as.integer(chunk_size),
    n_threads     = as.integer(n_threads),
    tunedG        = as.integer(tunedG),
    A22           = A22
  )
}
