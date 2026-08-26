#' Trace the decomposition as administrative censoring increases
#'
#' @description
#' Sweeps a series of administrative cut-off times over the cohort. At each
#' cut-off, subjects whose event falls after the cut-off are treated as
#' censored there, and the decomposition is recomputed. This traces how the
#' global C-index behaves as censoring rises, while `C_ee` -- the harder
#' ranking task -- is tracked separately.
#'
#' The cut-off is applied to the **whole cohort**, so the reported
#' censoring rate is the cohort's actual censoring rate at that cut-off.
#'
#' @param time,status,risk Numeric vectors of equal length.
#' @param weights A `cindex_weights` object. Defaults to
#'   [weights_harrell()].
#' @param higher_is_riskier Logical; see [decompose_cindex()].
#' @param n_thresholds Number of cut-offs to evaluate. Ignored if `probs`
#'   is supplied.
#' @param probs Quantiles of the observed event times at which to place the
#'   cut-offs.
#' @param min_pairs Minimum number of event-event pairs for an estimate to
#'   be considered precise. Rows below it are **flagged** in the
#'   `low_precision` column, never removed and never set to `NA`. `NULL`
#'   (default) selects `max(50, round(n_events^1.5 * 0.1))`.
#'
#' @return An object of class `cindex_curve`. Its `$data` element has one
#'   row per threshold, with columns `threshold`, `censoring`, `C_ee`,
#'   `C_ec`, `C_global`, `gap`, `N_ee`, `N_ec`, `low_precision`, `W_ee` and
#'   `W_ec`. `W_ee`/`W_ec` are the weighted totals behind `C_ee`/`C_ec`
#'   (`C_ee = S_ee / W_ee`, etc.), so the decomposition identity is
#'   verifiable at every row:
#'   `C_global == (W_ee * C_ee + W_ec * C_ec) / (W_ee + W_ec)`. `N_ee`/
#'   `N_ec` are pair *counts* and only equal `W_ee`/`W_ec` under a weighting
#'   that assigns every comparable pair weight 1 (Harrell's); under
#'   [weights_uno()] or another non-unit weighting the identity does not
#'   reconstruct from counts alone.
#'
#' @examples
#' set.seed(42)
#' time   <- rexp(300, rate = 0.1)
#' status <- rbinom(300, 1, 0.6)
#' risk   <- rnorm(300)
#' censoring_curve(time, status, risk, n_thresholds = 8)
#'
#' @export
censoring_curve <- function(time, status, risk,
                            weights = weights_harrell(),
                            higher_is_riskier = TRUE,
                            n_thresholds = 20,
                            probs = seq(0.05, 0.70, length.out = n_thresholds),
                            min_pairs = NULL) {
  validate_survival_inputs(time, status, risk)
  if (!inherits(weights, "cindex_weights")) {
    stop("`weights` must be a `cindex_weights` object.", call. = FALSE)
  }
  if (!is.logical(higher_is_riskier) || length(higher_is_riskier) != 1L) {
    stop("`higher_is_riskier` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!higher_is_riskier) risk <- -risk

  event_idx <- which(status == 1)
  n_events <- length(event_idx)
  if (is.null(min_pairs)) {
    min_pairs <- max(50, round(n_events^1.5 * 0.1))
  }

  taus <- unique(stats::quantile(time[event_idx], probs = probs, names = FALSE))

  rows <- lapply(taus, function(tau) {
    # The cut-off applies to EVERYONE. Subjects already censored stay
    # censored; events after the cut-off become censored at it.
    sim_time <- pmin(time, tau)
    sim_status <- status * as.numeric(time <= tau)
    if (sum(sim_status) < 2) return(NULL)

    pc <- pair_counts(sim_time, sim_status, risk, weights)
    total_w <- pc$W_ee + pc$W_ec
    if (total_w <= 0) return(NULL)

    C_ee <- if (pc$W_ee > 0) pc$S_ee / pc$W_ee else NA_real_
    C_ec <- if (pc$W_ec > 0) pc$S_ec / pc$W_ec else NA_real_

    data.frame(
      threshold = tau,
      censoring = mean(sim_status == 0),
      C_ee = C_ee,
      C_ec = C_ec,
      C_global = (pc$S_ee + pc$S_ec) / total_w,
      gap = C_ec - C_ee,
      N_ee = pc$N_ee,
      N_ec = pc$N_ec,
      low_precision = pc$N_ee < min_pairs,
      W_ee = pc$W_ee,
      W_ec = pc$W_ec
    )
  })

  dat <- do.call(rbind, rows)
  if (is.null(dat) || nrow(dat) == 0L) {
    stop("No threshold produced enough events to decompose. ",
         "Try a smaller `probs` range.", call. = FALSE)
  }
  dat <- dat[order(dat$censoring), , drop = FALSE]
  rownames(dat) <- NULL

  structure(
    list(
      data = dat,
      weighting = weights$name,
      higher_is_riskier = higher_is_riskier,
      min_pairs = min_pairs,
      n = length(time),
      n_events = n_events
    ),
    class = "cindex_curve"
  )
}

#' @param x A `cindex_curve` object.
#' @param ... Ignored.
#' @rdname censoring_curve
#' @export
print.cindex_curve <- function(x, ...) {
  cat("Censoring Curve\n")
  cat(sprintf("Weighting: %s | n = %d, events = %d | %d thresholds\n",
              x$weighting, x$n, x$n_events, nrow(x$data)))
  cat(sprintf("Censoring range: %.1f%% to %.1f%%\n\n",
              100 * min(x$data$censoring), 100 * max(x$data$censoring)))
  show <- x$data[, c("censoring", "C_ee", "C_ec", "C_global", "N_ee",
                     "low_precision")]
  show$censoring <- sprintf("%.1f%%", 100 * show$censoring)
  print(show, digits = 4, row.names = FALSE)
  if (any(x$data$low_precision)) {
    cat(sprintf("\n%d row(s) flagged low_precision (N_ee < %d).",
                sum(x$data$low_precision), x$min_pairs))
    cat(" Estimates are shown, not removed.\n")
  }
  invisible(x)
}

#' @param row.names Ignored.
#' @param optional Ignored.
#' @rdname censoring_curve
#' @export
as.data.frame.cindex_curve <- function(x, row.names = NULL,
                                       optional = FALSE, ...) {
  x$data
}
