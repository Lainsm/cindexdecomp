#' Subject-level bootstrap for a concordance decomposition
#'
#' Resamples **subjects** with replacement, not pairs. Pairs are strongly
#' dependent -- each subject appears in hundreds of them -- so resampling
#' pairs directly produces intervals several times too narrow.
#'
#' All four quantities are recomputed inside each replicate, so the
#' correlation between `C_ee` and `C_ec` is carried through to the gap
#' automatically.
#'
#' @param time,status,risk Numeric vectors of equal length.
#' @param weights A `cindex_weights` object.
#' @param n_boot Number of replicates.
#' @return A numeric matrix, `n_boot` rows, columns `C_ee`, `C_ec`,
#'   `C_global`, `gap`. Replicates that yield no comparable pairs of one
#'   kind are left as `NA`. The count of usable (complete) replicates is
#'   attached as the `"n_valid"` attribute.
#' @keywords internal
#' @noRd
bootstrap_decomp <- function(time, status, risk, weights, n_boot) {
  n <- length(time)
  out <- matrix(
    NA_real_, nrow = n_boot, ncol = 4L,
    dimnames = list(NULL, c("C_ee", "C_ec", "C_global", "gap"))
  )
  for (b in seq_len(n_boot)) {
    k <- sample.int(n, n, replace = TRUE)
    pc <- tryCatch(
      pair_counts(time[k], status[k], risk[k], weights),
      error = function(e) NULL
    )
    if (is.null(pc) || pc$W_ee <= 0 || pc$W_ec <= 0) next
    ee <- pc$S_ee / pc$W_ee
    ec <- pc$S_ec / pc$W_ec
    out[b, ] <- c(ee, ec, (pc$S_ee + pc$S_ec) / (pc$W_ee + pc$W_ec), ec - ee)
  }
  attr(out, "n_valid") <- sum(stats::complete.cases(out))
  out
}

#' Bootstrap confidence intervals for a concordance decomposition
#'
#' @param object A `cindex_decomp` object created with `n_boot > 0`.
#' @param parm Character vector selecting quantities. Defaults to all of
#'   `"C_ee"`, `"C_ec"`, `"C_global"`, `"gap"`.
#' @param level Confidence level. Defaults to the object's `conf_level`.
#' @param ... Ignored.
#' @return A matrix of percentile bounds, one row per quantity.
#' @examples
#' set.seed(1)
#' time   <- rexp(150, rate = 0.1)
#' status <- rbinom(150, 1, 0.6)
#' risk   <- rnorm(150)
#' fit <- decompose_cindex(time, status, risk, n_boot = 50)
#' confint(fit)
#' @export
confint.cindex_decomp <- function(object,
                                  parm = c("C_ee", "C_ec", "C_global", "gap"),
                                  level = NULL, ...) {
  if (is.null(object$boot)) {
    stop(
      "No bootstrap replicates are stored on this object. ",
      "Re-run decompose_cindex() with n_boot > 0 (e.g. n_boot = 1000).",
      call. = FALSE
    )
  }
  if (is.null(level)) level <- object$conf_level
  parm <- match.arg(parm, c("C_ee", "C_ec", "C_global", "gap"),
                    several.ok = TRUE)
  n_valid <- sum(stats::complete.cases(object$boot))
  if (n_valid < 0.9 * nrow(object$boot)) {
    warning(
      sprintf(
        "Only %d of %d bootstrap replicates were usable (%.0f%%). Intervals are based on the usable ones; with few events, consider raising n_boot.",
        n_valid, nrow(object$boot), 100 * n_valid / nrow(object$boot)
      ),
      call. = FALSE
    )
  }
  a <- (1 - level) / 2
  probs <- c(a, 1 - a)
  out <- t(vapply(
    parm,
    function(p) stats::quantile(object$boot[, p], probs = probs,
                                na.rm = TRUE, names = FALSE),
    numeric(2)
  ))
  colnames(out) <- paste0(format(100 * probs, trim = TRUE), " %")
  rownames(out) <- parm
  out
}
