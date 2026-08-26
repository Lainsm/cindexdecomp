#' Compare the decomposition across several models
#'
#' @description
#' Runs [decompose_cindex()] for each of several risk scores measured on
#' the same cohort, and assembles the results into one table with bootstrap
#' standard deviations. This is the input consumed by
#' `autoplot()` to draw the dumbbell comparison.
#'
#' @param risks A **named** list or data frame of risk-score vectors, each
#'   the same length as `time` and `status`. The names label the models.
#' @param time,status Numeric vectors describing the shared cohort.
#' @param weights A `cindex_weights` object. Defaults to
#'   [weights_harrell()].
#' @param higher_is_riskier Logical; see [decompose_cindex()]. Applied to
#'   every model.
#' @param n_boot Bootstrap replicates per model. Defaults to 1000, since
#'   the comparison table exists to carry uncertainty.
#' @param conf_level Confidence level stored on each fit.
#'
#' @return An object of class `cindex_comparison`.
#'
#' @examples
#' set.seed(42)
#' n <- 200
#' time   <- rexp(n, rate = 0.1)
#' status <- rbinom(n, 1, 0.6)
#' risks  <- list(model_a = rnorm(n), model_b = rnorm(n))
#' compare_decompositions(risks, time, status, n_boot = 0)
#'
#' @export
compare_decompositions <- function(risks, time, status,
                                   weights = weights_harrell(),
                                   higher_is_riskier = TRUE,
                                   n_boot = 1000,
                                   conf_level = 0.95) {
  if (is.data.frame(risks)) risks <- as.list(risks)
  if (!is.list(risks) || length(risks) == 0L) {
    stop("`risks` must be a non-empty named list or data frame of risk scores.",
         call. = FALSE)
  }
  if (is.null(names(risks)) || any(!nzchar(names(risks)))) {
    stop("`risks` must be fully named; the names label the models.",
         call. = FALSE)
  }
  bad <- vapply(risks, function(r) length(r) != length(time), logical(1))
  if (any(bad)) {
    stop("Every risk score must be the same length as `time` and `status`. ",
         "Offending: ", paste(names(risks)[bad], collapse = ", "), ".",
         call. = FALSE)
  }

  fits <- lapply(risks, function(r) {
    decompose_cindex(
      time, status, r,
      weights = weights,
      higher_is_riskier = higher_is_riskier,
      n_boot = n_boot,
      conf_level = conf_level
    )
  })
  names(fits) <- names(risks)

  boot_sd <- function(fit, col) {
    if (is.null(fit$boot)) NA_real_ else stats::sd(fit$boot[, col], na.rm = TRUE)
  }

  tab <- data.frame(
    model = names(fits),
    ci_ee = vapply(fits, function(f) f$C_ee, numeric(1)),
    ci_ec = vapply(fits, function(f) f$C_ec, numeric(1)),
    global_c = vapply(fits, function(f) f$C_global, numeric(1)),
    gap = vapply(fits, function(f) f$gap, numeric(1)),
    sd_ee = vapply(fits, boot_sd, numeric(1), col = "C_ee"),
    sd_ec = vapply(fits, boot_sd, numeric(1), col = "C_ec"),
    sd_global = vapply(fits, boot_sd, numeric(1), col = "C_global"),
    sd_gap = vapply(fits, boot_sd, numeric(1), col = "gap"),
    n_ee = vapply(fits, function(f) f$N_ee, numeric(1)),
    n_ec = vapply(fits, function(f) f$N_ec, numeric(1)),
    stringsAsFactors = FALSE
  )
  rownames(tab) <- NULL

  structure(
    list(
      table = tab,
      fits = fits,
      weighting = weights$name,
      n_boot = n_boot
    ),
    class = "cindex_comparison"
  )
}

#' @param x A `cindex_comparison` object.
#' @param ... Ignored.
#' @rdname compare_decompositions
#' @export
print.cindex_comparison <- function(x, ...) {
  cat("C-Index Decomposition Comparison\n")
  cat(sprintf("Weighting: %s | %d model(s) | n_boot = %d\n\n",
              x$weighting, nrow(x$table), x$n_boot))
  show <- x$table[, c("model", "ci_ee", "ci_ec", "global_c", "gap")]
  names(show) <- c("model", "C_ee", "C_ec", "C_global", "gap")
  print(show, digits = 4, row.names = FALSE)
  if (x$n_boot == 0) {
    cat("\nNo bootstrap: sd_* columns are NA. Set n_boot > 0 for uncertainty.\n")
  }
  invisible(x)
}

#' @param row.names Ignored.
#' @param optional Ignored.
#' @rdname compare_decompositions
#' @export
as.data.frame.cindex_comparison <- function(x, row.names = NULL,
                                            optional = FALSE, ...) {
  x$table
}
