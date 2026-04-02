#' Renumber a pedigree in topological order
#'
#' Assigns consecutive integer IDs (1..N) to animals in a pedigree, ensuring
#' that parents always receive a lower ID than their offspring. This is a
#' prerequisite for building the numerator relationship matrix inverse (A-inverse)
#' via [build_Ainv()].
#'
#' @param pedigree A data frame or data.table with exactly three columns in order:
#'   animal ID, sire ID, dam ID. Column names are ignored; only position matters.
#'   IDs can be character or numeric.
#' @param unknown_parent Character vector of codes used to denote unknown parents.
#'   Defaults to `c("0", "00000000000000")`. All codes are recoded to `"0"` internally.
#' @param verbose Logical. If `TRUE` (default), prints iteration progress to stderr.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{new_id}{Integer. New consecutive ID for the animal (1..N).}
#'   \item{new_sire}{Integer. New ID of the sire (0 = unknown).}
#'   \item{new_dam}{Integer. New ID of the dam (0 = unknown).}
#'   \item{animal}{Character. Original animal ID.}
#'   \item{sire}{Character. Original sire ID.}
#'   \item{dam}{Character. Original dam ID.}
#' }
#' Rows are sorted by `new_id`. Animals that could not be renumbered (due to
#' missing parents or cycles) are excluded with a warning.
#'
#' @details
#' The algorithm performs iterative topological sorting: in each iteration it
#' includes all animals whose parents have already been assigned an ID. This
#' handles arbitrary pedigree depths without recursion.
#'
#' @examples
#' ped <- data.frame(
#'   animal = c("A", "B", "C", "D"),
#'   sire   = c("0", "0", "A", "A"),
#'   dam    = c("0", "0", "B", "B")
#' )
#' renum_pedigree(ped)
#'
#' @export
renum_pedigree <- function(pedigree,
                           unknown_parent = c("0", "00000000000000"),
                           verbose = TRUE) {

  if (!requireNamespace("data.table", quietly = TRUE))
    stop("Package 'data.table' is required. Install it with install.packages('data.table').")
  if (!is.data.frame(pedigree) || ncol(pedigree) < 3)
    stop("'pedigree' must be a data.frame-like object with at least 3 columns: animal, sire, dam.")

  ped <- data.table::as.data.table(pedigree)
  data.table::setnames(ped, 1:3, c("animal", "sire", "dam"))
  ped[, c("animal", "sire", "dam") := lapply(.SD, as.character), .SDcols = 1:3]

  # Recode unknown parents
  for (code in unknown_parent) {
    ped[animal == code, animal := "0"]
    ped[sire   == code, sire   := "0"]
    ped[dam    == code, dam    := "0"]
  }
  ped[is.na(animal), animal := "0"]
  ped[is.na(sire),   sire   := "0"]
  ped[is.na(dam),    dam    := "0"]
  if (anyDuplicated(ped$animal[ped$animal != "0"]))
    stop("Duplicated animal IDs found in 'pedigree'. Each animal must appear only once.")

  ids      <- ped$animal
  all_ids  <- unique(c("0", ids, ped$sire, ped$dam))
  included <- stats::setNames(rep(FALSE, length(all_ids)), all_ids)
  included["0"] <- TRUE

  included_row <- logical(nrow(ped))
  new_id       <- integer(length(ids))
  names(new_id) <- ids

  cur  <- 0L
  iter <- 0L

  repeat {
    iter <- iter + 1L
    sire_inc <- ped$sire == "0" | included[ped$sire]
    dam_inc  <- ped$dam  == "0" | included[ped$dam]
    can      <- !included_row & sire_inc & dam_inc
    if (!any(can)) break

    idx   <- which(can)
    k     <- length(idx)
    codes <- seq.int(cur + 1L, cur + k)
    cur   <- cur + k

    animals_iter <- ped$animal[idx]
    new_id[animals_iter]   <- codes
    included_row[idx]      <- TRUE
    included[animals_iter] <- TRUE

    if (verbose)
      message(sprintf("  Iteration %4d: %6d animals included (total: %8d)", iter, k, cur))
  }

  if (any(!included_row)) {
    warning(sprintf(
      "%d animals could not be renumbered (missing parents or cycles in pedigree).",
      sum(!included_row)
    ))
  }

  sire_code <- new_id[ped$sire]; sire_code[is.na(sire_code)] <- 0L
  dam_code  <- new_id[ped$dam];  dam_code[is.na(dam_code)]   <- 0L

  res <- data.frame(
    new_id   = new_id[ped$animal],
    new_sire = sire_code,
    new_dam  = dam_code,
    animal   = ped$animal,
    sire     = ped$sire,
    dam      = ped$dam,
    stringsAsFactors = FALSE
  )
  res[order(res$new_id), ]
}
