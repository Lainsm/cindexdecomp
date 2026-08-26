fmt_n <- function(x) formatC(x, format = "d", big.mark = ",")

#' @param x A `cindex_decomp` object.
#' @param ... Ignored.
#' @rdname decompose_cindex
#' @export
print.cindex_decomp <- function(x, ...) {
  ci <- if (!is.null(x$boot)) confint(x) else NULL
  total_n <- x$N_ee + x$N_ec

  cat("C-Index Decomposition\n")
  cat(sprintf(
    "Weighting: %s | n = %d, events = %d (%.1f%%), censoring %.1f%%\n\n",
    x$weighting, x$n, x$n_events,
    100 * x$n_events / x$n, 100 * x$censoring_rate
  ))

  row <- function(label, value, npairs, key) {
    share <- if (total_n > 0) 100 * npairs / total_n else NA_real_
    line <- sprintf("  %-16s %8.4f %11s %7.1f%%",
                    label, value, fmt_n(npairs), share)
    if (!is.null(ci) && key %in% rownames(ci)) {
      line <- paste0(line, sprintf("   [%.4f, %.4f]", ci[key, 1], ci[key, 2]))
    }
    cat(line, "\n")
  }

  header <- sprintf("  %-16s %8s %11s %8s", "", "C-index", "pairs", "share")
  if (!is.null(ci)) {
    header <- paste0(header, sprintf("   %s CI",
                                     format(100 * x$conf_level, trim = TRUE)))
  }
  cat(header, "\n")
  row("Event-Event", x$C_ee, x$N_ee, "C_ee")
  row("Event-Censored", x$C_ec, x$N_ec, "C_ec")
  cat("  ", strrep("-", if (is.null(ci)) 46 else 66), "\n", sep = "")
  row("Global", x$C_global, total_n, "C_global")

  cat(sprintf("\n  Masking gap (C_ec - C_ee): %+.3f", x$gap))
  if (!is.null(ci)) {
    cat(sprintf("   [%+.3f, %+.3f]", ci["gap", 1], ci["gap", 2]))
  }
  cat("\n")

  if (!is.na(x$C_ee) && x$C_ee < 0.55) {
    cat("  ! Event-Event concordance is near chance (0.50).\n")
  }
  if (!is.na(x$C_global) && x$C_global < 0.5) {
    cat("  ! Global C-index is below 0.50. If the risk score is oriented\n")
    cat("    so that higher means lower hazard, set higher_is_riskier = FALSE.\n")
  }
  if (is.null(ci)) {
    cat("  Pair counts describe composition, not precision.",
        "Set n_boot > 0 for intervals.\n")
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
    cat(sprintf("%d replicates (%d usable), subject-level resampling, %.0f%% level\n",
                x$fit$n_boot, n_ok, 100 * x$fit$conf_level))
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
