#' Print method for connectedness objects
#'
#' @param x An object of class `"connectedness"`.
#' @param digits Integer. Number of decimal places for CD and PEVD matrices.
#'   Default is 3.
#' @param ... Further arguments (currently ignored).
#'
#' @export
print.connectedness <- function(x, digits = 3, ...) {
  cat("Connectedness analysis\n")
  cat("Call: "); print(x$call); cat("\n")

  if (!is.null(x$year_window))
    cat(sprintf("Year window: [%d, %d]\n\n", x$year_window[1], x$year_window[2]))

  cat("Target animals per MU:\n")
  print(x$n_target)

  cat("\nCoefficient of Determination (CD):\n")
  print(round(x$CD, digits))

  cat("\nPEV of Differences (PEVD):\n")
  print(round(x$PEVD, digits))

  invisible(x)
}


#' Plot method for connectedness objects
#'
#' Produces two plots: (1) a heatmap of the CD matrix between MU pairs, and
#' (2) if temporal data is available, a dot plot showing the years of overlap
#' between each MU pair.
#'
#' @param x An object of class `"connectedness"`.
#' @param which Character. Which plot to produce: `"CD"` (default), `"PEVD"`,
#'   `"overlap"`, or `"all"`.
#' @param ... Further graphical parameters passed to [graphics::image()] or
#'   [graphics::plot()].
#'
#' @export
plot.connectedness <- function(x, which = c("CD", "PEVD", "overlap", "all"), ...) {

  which <- match.arg(which)

  do_cd      <- which %in% c("CD", "all")
  do_pevd    <- which %in% c("PEVD", "all")
  do_overlap <- which %in% c("overlap", "all") && !is.null(x$overlap)

  n_plots <- sum(c(do_cd, do_pevd, do_overlap))
  if (n_plots == 0) {
    message("Nothing to plot. Use which = 'CD', 'PEVD', or 'overlap'.")
    return(invisible(x))
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  if (n_plots > 1)
    graphics::par(mfrow = c(1, n_plots))

  mu_names <- rownames(x$CD)
  U <- nrow(x$CD)

  # -- CD heatmap --
  if (do_cd) {
    mat <- x$CD
    mat[is.na(mat)] <- 0
    graphics::image(
      1:U, 1:U, mat,
      col  = grDevices::colorRampPalette(c("white", "#2171b5"))(50),
      xaxt = "n", yaxt = "n",
      xlab = "", ylab = "",
      main = "Coefficient of Determination (CD)",
      zlim = c(0, 1),
      ...
    )
    graphics::axis(1, at = 1:U, labels = mu_names, las = 2)
    graphics::axis(2, at = 1:U, labels = mu_names, las = 1)
    for (i in 1:U) for (j in 1:U) {
      v <- x$CD[i, j]
      if (!is.na(v))
        graphics::text(i, j, sprintf("%.2f", v), cex = 0.8)
    }
  }

  # -- PEVD heatmap --
  if (do_pevd) {
    mat <- x$PEVD
    mat[is.na(mat)] <- 0
    graphics::image(
      1:U, 1:U, mat,
      col  = grDevices::colorRampPalette(c("white", "#cb181d"))(50),
      xaxt = "n", yaxt = "n",
      xlab = "", ylab = "",
      main = "PEV of Differences (PEVD)",
      ...
    )
    graphics::axis(1, at = 1:U, labels = mu_names, las = 2)
    graphics::axis(2, at = 1:U, labels = mu_names, las = 1)
    for (i in 1:U) for (j in 1:U) {
      v <- x$PEVD[i, j]
      if (!is.na(v))
        graphics::text(i, j, sprintf("%.2f", v), cex = 0.8)
    }
  }

  # -- Overlap dot plot --
  if (do_overlap) {
    dt <- x$overlap
    pairs <- unique(paste(dt$MU1, dt$MU2, sep = "\u2013"))
    pairs <- pairs[order(pairs)]
    y_pos <- match(paste(dt$MU1, dt$MU2, sep = "\u2013"), pairs)

    graphics::par(mar = c(5, 10, 3, 2))
    graphics::plot(
      dt$Year, y_pos,
      pch  = 15, cex = 0.8,
      yaxt = "n",
      xlab = "Year", ylab = "",
      main = "Temporal overlap between MU pairs",
      ...
    )
    graphics::axis(2, at = seq_along(pairs), labels = pairs, las = 2)
    graphics::grid(nx = NA, ny = NULL, lty = 3)
  }

  invisible(x)
}
