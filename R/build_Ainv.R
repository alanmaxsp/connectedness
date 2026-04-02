#' Build the sparse inverse of the numerator relationship matrix (A-inverse)
#'
#' Constructs A-inverse directly from a renumbered pedigree using the
#' Henderson rules, without ever forming A explicitly. Inbreeding coefficients
#' are computed via the recursive algorithm of Aguilar & Misztal.
#'
#' @param renum A data frame as returned by [renum_pedigree()], with columns
#'   `new_id`, `new_sire`, and `new_dam` (integer, 0 = unknown).
#'
#' @return A list with:
#' \describe{
#'   \item{Ainv}{A sparse matrix of class `dgCMatrix` (N x N), the inverse of
#'     the numerator relationship matrix.}
#'   \item{F}{Numeric vector of length N with individual inbreeding coefficients
#'     (0 = non-inbred).}
#' }
#'
#' @details
#' The pedigree must be in topological order (parents before offspring), as
#' produced by [renum_pedigree()]. The function calls compiled C++ code via
#' Rcpp and RcppEigen for efficiency with large pedigrees.
#'
#' @seealso [renum_pedigree()], [compute_connectedness()]
#'
#' @examples
#' ped <- data.frame(
#'   animal = c("A", "B", "C", "D"),
#'   sire   = c("0", "0", "A", "A"),
#'   dam    = c("0", "0", "B", "B")
#' )
#' renum <- renum_pedigree(ped, verbose = FALSE)
#' Ainv_res <- build_Ainv(renum)
#' Ainv_res$Ainv
#'
#' @export
build_Ainv <- function(renum) {

  required_cols <- c("new_id", "new_sire", "new_dam")
  if (!all(required_cols %in% names(renum)))
    stop("'renum' must have columns: new_id, new_sire, new_dam. ",
         "Use renum_pedigree() to generate it.")

  renum <- renum[order(renum$new_id), ]

  sire_vec <- as.integer(renum$new_sire)
  dam_vec  <- as.integer(renum$new_dam)

  if (any(is.na(sire_vec)) || any(is.na(dam_vec)))
    stop("NA values found in new_sire or new_dam. Check pedigree renumbering.")

  build_Ainv_sparse_RA(sire = sire_vec, dam = dam_vec, cache_parent_pairs = TRUE)
}
