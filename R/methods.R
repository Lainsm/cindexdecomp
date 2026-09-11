fmt_n <- function(x) formatC(x, format = "d", big.mark = ",")
fmt_level <- function(p) paste0(format(100 * p, trim = TRUE), "%")

#' @param x A `cindex_decomp` object.
#' @param ... Ignored.
#' @rdname decompose_cindex
#' @export
print.cindex_decomp <- function(x, ...) {
  ci <- if (!is.null(x$boot)) stats::confint(x) else NULL
  total_n <- x$N_ee + x$N_ec

  cat("C-Index Decomposition\n")
  cat(sprintf(
    "Weighting: %s | n = %d, events = %d (%.1f%%), censoring %.1f%%\n\n",
    x$weighting, x$n, x$n_events,
    100 * x$n_events / x$n, 100 * x$censoring_rate
  ))

  row <- function(label, value, npairs, key) {
    share <- if (total_n > 0) 100 * npairs / total_n else NA_real_
    line <- sprintf(
      "  %-16s %8.4f %11s %7.1f%%",
      label, value, fmt_n(npairs), share
    )
    if (!is.null(ci) && key %in% rownames(ci)) {
      line <- paste0(line, sprintf("   [%.4f, %.4f]", ci[key, 1], ci[key, 2]))
    }
    cat(line, "\n", sep = "")
  }

  header <- sprintf("  %-16s %8s %11s %8s", "", "C-index", "pairs", "share")
  if (!is.null(ci)) {
    header <- paste0(header, sprintf("   %s CI", fmt_level(x$conf_level)))
  }
  cat(header, "\n", sep = "")
  row("Event-Event", x$C_ee, x$N_ee, "C_ee")
  row("Event-Censored", x$C_ec, x$N_ec, "C_ec")
  cat("  ", strrep("-", if (is.null(ci)) 46 else 65), "\n", sep = "")
  row("Global", x$C_global, total_n, "C_global")

  cat(sprintf("\n  Masking gap (C_ec - C_ee): %+.3f", x$gap))
  if (!is.null(ci)) {
    cat(sprintf("   [%+.3f, %+.3f]", ci["gap", 1], ci["gap", 2]))
  }
  cat("\n")

  if (!is.na(x$C_ee)) {
    if (abs(x$C_ee - 0.5) < 0.05) {
      cat("  ! Event-Event concordance is near chance (0.50).\n")
    } else if (x$C_ee < 0.45) {
      cat("  ! Event-Event concordance is well below chance, indicating\n")
      cat("    systematically reversed ranking rather than absent signal.\n")
    }
  }
  if (!is.na(x$C_global) && x$C_global < 0.5) {
    cat("  ! Global C-index is below 0.50. If the risk score is oriented\n")
    cat(
      "    so that higher means lower hazard, set",
      "higher_is_riskier = FALSE.\n"
    )
  }
  if (!is.null(x$boot) && !is.na(x$n_boot_valid) &&
    x$n_boot_valid < 0.9 * x$n_boot) {
    cat(sprintf(
      "  ! Only %d of %d bootstrap replicates were usable (%.0f%%).\n",
      x$n_boot_valid, x$n_boot, 100 * x$n_boot_valid / x$n_boot
    ))
    cat("    Intervals above are unreliable; consider more events or n_boot.\n")
  }
  if (is.null(ci)) {
    cat(
      "  Pair counts describe composition, not precision.",
      "Set n_boot > 0 for intervals.\n"
    )
  }
  invisible(x)
}

#' @param object A `cindex_decomp` object.
#' @rdname decompose_cindex
#' @export
summary.cindex_decomp <- function(object, ...) {
  structure(list(fit = object), class = "summary.cindex_decomp")
}

#' @export
print.summary.cindex_decomp <- function(x, ...) {
  print(x$fit)
  cat("\nBootstrap: ")
  if (is.null(x$fit$boot)) {
    cat("none (n_boot = 0)\n")
  } else {
    n_ok <- x$fit$n_boot_valid
    cat(sprintf(
      "%d replicates (%d usable), subject-level resampling, %s level\n",
      x$fit$n_boot, n_ok, fmt_level(x$fit$conf_level)
    ))
  }
  cat(sprintf("Orientation: higher_is_riskier = %s\n", x$fit$higher_is_riskier))
  invisible(x)
}

#' @param row.names Ignored.
#' @param optional Ignored.
#' @rdname decompose_cindex
#' @export
as.data.frame.cindex_decomp <- function(x, row.names = NULL,
                                        optional = FALSE, ...) {
  data.frame(
    weighting = x$weighting,
    n = x$n,
    n_events = x$n_events,
    censoring_rate = x$censoring_rate,
    C_ee = x$C_ee,
    C_ec = x$C_ec,
    C_global = x$C_global,
    gap = x$gap,
    N_ee = x$N_ee,
    N_ec = x$N_ec,
    stringsAsFactors = FALSE
  )
}
