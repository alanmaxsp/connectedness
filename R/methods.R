#' Print a connectedness object
#'
#' @param x An object of class `"connectedness"`.
#' @param digits Number of digits used in printed summaries.
#' @param ... Unused.
#'
#' @return The input object, invisibly.
#' @export
print.connectedness <- function(x, digits = 3, ...) {

  if (!inherits(x, "connectedness")) {
    stop("'x' must be an object of class 'connectedness'.")
  }

  .fmt_range <- function(mat, digits = 3) {
    vals <- as.numeric(mat)
    vals <- vals[is.finite(vals)]
    if (!length(vals)) return("all NA")
    paste0(
      formatC(min(vals), digits = digits, format = "f"),
      " to ",
      formatC(max(vals), digits = digits, format = "f")
    )
  }

  cat("Connectedness analysis\n")
  cat("----------------------\n")
  cat("Call: ")
  print(x$call)
  cat("\n")

  if (!is.null(x$relationship)) {
    cat("Relationship matrix:", x$relationship, "\n")
  }

  if (!is.null(x$year_window)) {
    cat(sprintf("Year window        : [%d, %d]\n", x$year_window[1], x$year_window[2]))
  }

  if (!is.null(x$n_target)) {
    nt <- as.numeric(x$n_target)
    nt <- nt[is.finite(nt)]
    if (length(nt)) {
      cat(sprintf("Management units   : %d\n", length(x$n_target)))
      cat(sprintf(
        "Target animals/MU  : min=%s, median=%s, max=%s\n",
        formatC(min(nt), digits = 0, format = "f"),
        formatC(stats::median(nt), digits = 1, format = "f"),
        formatC(max(nt), digits = 0, format = "f")
      ))
    }
    cat("\nTarget animals per MU:\n")
    print(x$n_target)
  }

  cat("\nCD range   : ", .fmt_range(x$CD, digits = digits), "\n", sep = "")
  cat("PEVD range : ", .fmt_range(x$PEVD, digits = digits), "\n", sep = "")

  invisible(x)
}


#' Plot a connectedness object
#'
#' Draws heatmaps of `CD` or `PEVD`, and optionally a temporal overlap plot.
#'
#' @param x An object of class `"connectedness"`.
#' @param which Which plot to produce: `"CD"`, `"PEVD"`, `"overlap"`, or `"all"`.
#' @param show_values Logical; if `TRUE`, print values inside the heatmap cells.
#' @param digits Number of decimal places shown inside cells.
#' @param triangle Which part of the matrix to display: `"full"`, `"upper"`, or `"lower"`.
#' @param ... Unused.
#'
#' @return The input object, invisibly.
#' @export
plot.connectedness <- function(x,
                               which = c("CD", "PEVD", "overlap", "all"),
                               show_values = TRUE,
                               digits = 2,
                               triangle = c("full", "upper", "lower"),
                               ...) {

  if (!inherits(x, "connectedness")) {
    stop("'x' must be an object of class 'connectedness'.")
  }

  which <- match.arg(which)
  triangle <- match.arg(triangle)

  do_cd      <- which %in% c("CD", "all")
  do_pevd    <- which %in% c("PEVD", "all")
  do_overlap <- which %in% c("overlap", "all") && !is.null(x$overlap)

  n_plots <- sum(c(do_cd, do_pevd, do_overlap))
  if (n_plots == 0) {
    message("Nothing to plot. Use which = 'CD', 'PEVD', or 'overlap'.")
    return(invisible(x))
  }

  .mat_to_long <- function(mat, metric, triangle = "full") {
    df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
    names(df) <- c("MU1", "MU2", "value")

    rn <- rownames(mat)
    cn <- colnames(mat)
    if (is.null(rn)) rn <- seq_len(nrow(mat))
    if (is.null(cn)) cn <- seq_len(ncol(mat))

    df$MU1 <- factor(df$MU1, levels = rn)
    df$MU2 <- factor(df$MU2, levels = cn)
    df$i <- as.integer(df$MU1)
    df$j <- as.integer(df$MU2)

    if (triangle == "upper") {
      df <- df[df$i <= df$j, , drop = FALSE]
    } else if (triangle == "lower") {
      df <- df[df$i >= df$j, , drop = FALSE]
    }

    df$label <- ifelse(is.na(df$value), "", formatC(df$value, digits = digits, format = "f"))
    df$metric <- metric
    df
  }

  .plot_one_base <- function(mat, main, palette_fun, zlim = NULL) {
    mu_names <- rownames(mat)
    U <- nrow(mat)
    if (is.null(mu_names)) mu_names <- seq_len(U)

    vals <- mat
    vals[is.na(vals)] <- NA_real_

    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mar = c(8, 8, 4, 2))

    if (is.null(zlim)) {
      finite_vals <- vals[is.finite(vals)]
      zlim <- if (length(finite_vals)) range(finite_vals) else c(0, 1)
    }

    image_mat <- t(vals[U:1, , drop = FALSE])
    graphics::image(
      x = seq_len(U),
      y = seq_len(U),
      z = image_mat,
      col = palette_fun(100),
      xaxt = "n",
      yaxt = "n",
      xlab = "",
      ylab = "",
      main = main,
      zlim = zlim
    )
    graphics::axis(1, at = seq_len(U), labels = mu_names, las = 2)
    graphics::axis(2, at = seq_len(U), labels = rev(mu_names), las = 2)
    graphics::box()

    if (show_values) {
      for (i in seq_len(U)) {
        for (j in seq_len(U)) {
          v <- vals[i, j]
          if (is.na(v)) next
          keep <- switch(
            triangle,
            full  = TRUE,
            upper = i <= j,
            lower = i >= j
          )
          if (keep) {
            graphics::text(
              x = j,
              y = U - i + 1,
              labels = formatC(v, digits = digits, format = "f"),
              cex = 0.8
            )
          }
        }
      }
    }
  }

  .plot_overlap_base <- function(dt) {
    dt <- as.data.frame(dt, stringsAsFactors = FALSE)
    if (!all(c("MU1", "MU2", "Year") %in% names(dt))) {
      stop("'overlap' must contain columns MU1, MU2, and Year.")
    }
    if (!nrow(dt)) {
      message("No temporal overlap to plot.")
      return(invisible(NULL))
    }

    dt$MU1 <- as.character(dt$MU1)
    dt$MU2 <- as.character(dt$MU2)
    dt$Year <- as.integer(dt$Year)

    pair <- paste(dt$MU1, dt$MU2, sep = "\u2013")
    pairs <- sort(unique(pair))
    y_pos <- as.numeric(match(pair, pairs))
    years <- as.numeric(dt$Year)

    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mar = c(5, 10, 4, 2))
    graphics::plot(
      x = years,
      y = y_pos,
      pch = 19,
      cex = 0.9,
      yaxt = "n",
      xlab = "Year",
      ylab = "",
      main = "Temporal overlap between MU pairs"
    )
    graphics::axis(2, at = seq_along(pairs), labels = pairs, las = 2)
    graphics::grid(nx = NA, ny = NULL, lty = 3)
  }

  use_ggplot <- requireNamespace("ggplot2", quietly = TRUE)

  if (use_ggplot) {
    plot_list <- list()

    if (do_cd) {
      df_cd <- .mat_to_long(x$CD, metric = "CD", triangle = triangle)
      p_cd <- ggplot2::ggplot(df_cd, ggplot2::aes(x = MU2, y = MU1, fill = value)) +
        ggplot2::geom_tile(color = "white", linewidth = 0.6) +
        ggplot2::scale_fill_gradient(
          low = "#f7fbff",
          high = "#08519c",
          limits = c(0, 1),
          na.value = "grey92",
          name = "CD"
        ) +
        ggplot2::coord_equal() +
        ggplot2::labs(
          title = "Coefficient of Determination (CD)",
          subtitle = paste("Relationship:", x$relationship),
          x = NULL,
          y = NULL
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          panel.grid = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
          legend.position = "right",
          plot.title = ggplot2::element_text(face = "bold")
        )

      if (show_values) {
        p_cd <- p_cd + ggplot2::geom_text(ggplot2::aes(label = label), size = 3.2)
      }
      plot_list[[length(plot_list) + 1L]] <- p_cd
    }

    if (do_pevd) {
      df_pevd <- .mat_to_long(x$PEVD, metric = "PEVD", triangle = triangle)
      p_pevd <- ggplot2::ggplot(df_pevd, ggplot2::aes(x = MU2, y = MU1, fill = value)) +
        ggplot2::geom_tile(color = "white", linewidth = 0.6) +
        ggplot2::scale_fill_gradient(
          low = "#fff5eb",
          high = "#b30000",
          na.value = "grey92",
          name = "PEVD"
        ) +
        ggplot2::coord_equal() +
        ggplot2::labs(
          title = "Prediction Error Variance of Differences (PEVD)",
          subtitle = paste("Relationship:", x$relationship),
          x = NULL,
          y = NULL
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          panel.grid = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
          legend.position = "right",
          plot.title = ggplot2::element_text(face = "bold")
        )

      if (show_values) {
        p_pevd <- p_pevd + ggplot2::geom_text(ggplot2::aes(label = label), size = 3.2)
      }
      plot_list[[length(plot_list) + 1L]] <- p_pevd
    }

    if (do_overlap) {
      dt <- as.data.frame(x$overlap, stringsAsFactors = FALSE)
      if (nrow(dt)) {
        dt$MU1 <- as.character(dt$MU1)
        dt$MU2 <- as.character(dt$MU2)
        dt$Year <- as.integer(dt$Year)
        dt$pair <- paste(dt$MU1, dt$MU2, sep = "\u2013")
        p_overlap <- ggplot2::ggplot(dt, ggplot2::aes(x = Year, y = stats::reorder(pair, Year))) +
          ggplot2::geom_point(size = 2.4, alpha = 0.9) +
          ggplot2::labs(
            title = "Temporal overlap between MU pairs",
            x = "Year",
            y = NULL
          ) +
          ggplot2::theme_minimal(base_size = 12) +
          ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(face = "bold")
          )
        plot_list[[length(plot_list) + 1L]] <- p_overlap
      }
    }

    if (length(plot_list) == 1L) {
      print(plot_list[[1L]])
    } else {
      if (requireNamespace("gridExtra", quietly = TRUE)) {
        gridExtra::grid.arrange(grobs = plot_list, ncol = length(plot_list))
      } else {
        for (p in plot_list) print(p)
      }
    }

  } else {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)

    if (n_plots > 1) graphics::par(mfrow = c(1, n_plots))

    if (do_cd) {
      .plot_one_base(
        mat = x$CD,
        main = paste0("CD (", x$relationship, ")"),
        palette_fun = function(n) grDevices::colorRampPalette(c("#f7fbff", "#08519c"))(n),
        zlim = c(0, 1)
      )
    }

    if (do_pevd) {
      .plot_one_base(
        mat = x$PEVD,
        main = paste0("PEVD (", x$relationship, ")"),
        palette_fun = function(n) grDevices::colorRampPalette(c("#fff5eb", "#b30000"))(n)
      )
    }

    if (do_overlap) {
      .plot_overlap_base(x$overlap)
    }
  }

  invisible(x)
}
