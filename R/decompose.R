#' Decompose a concordance index into Event-Event and Event-Censored parts
#'
#' @description
#' Splits a concordance index into the two ranking tasks it pools: pairs
#' where both subjects had the event (`C_ee`), and pairs where one had the
#' event and one was censored (`C_ec`). As censoring rises the event-event
#' pairs vanish and `C_ec` dominates the global score, so a model close to
#' chance on true events can still report a stable global C-index.
#'
#' The decomposition applies to any concordance index of the form
#' `C = sum(w_ij * c_ij) / sum(w_ij)`. The weighting is chosen with
#' [weights_harrell()], [weights_uno()], [weights_truncated()] or
#' [weights_custom()], and the identity
#' `C = (W_ee * C_ee + W_ec * C_ec) / (W_ee + W_ec)` holds for each.
#'
#' @param x Either a numeric vector of observed times (default method) or a
#'   formula of the form `Surv(time, status) ~ risk` (formula method).
#' @param status Numeric vector of event indicators (1 = event,
#'   0 = censored). Default method only.
#' @param risk Numeric vector of risk scores. Default method only.
#' @param data A data frame in which to evaluate the formula. Formula
#'   method only.
#' @param weights A `cindex_weights` object. Defaults to
#'   [weights_harrell()].
#' @param higher_is_riskier Logical. `TRUE` (default) if larger values of
#'   `risk` indicate higher hazard. If `FALSE`, `risk` is negated once
#'   before any pair counting, so every component stays mutually
#'   consistent. Values below 0.5 are reported as they are: a
#'   worse-than-chance model is a finding, not an error.
#' @param n_boot Number of bootstrap replicates. `0` (default) skips
#'   resampling, in which case [confint()] is unavailable. Replicates
#'   resample subjects, not pairs.
#' @param conf_level Confidence level stored on the object and used as the
#'   default for [confint()].
#' @param ... Passed between methods.
#'
#' @return An object of class `cindex_decomp`.
#'
#' @examples
#' set.seed(42)
#' time   <- rexp(200, rate = 0.1)
#' status <- rbinom(200, 1, 0.6)
#' risk   <- rnorm(200)
#'
#' decompose_cindex(time, status, risk)
#' decompose_cindex(time, status, risk, weights = weights_uno())
#'
#' @export
decompose_cindex <- function(x, ...) {
  UseMethod("decompose_cindex")
}

#' @rdname decompose_cindex
#' @export
decompose_cindex.default <- function(x, status, risk,
                                     weights = weights_harrell(),
                                     higher_is_riskier = TRUE,
                                     n_boot = 0, conf_level = 0.95, ...) {
  # `x` is required for S3 dispatch, so `decompose_cindex(time = t, ...)`
  # leaves `x` unmatched and `time` falls through into `...` here (checked
  # before `x` is touched, since referencing a missing `x` errors first).
  if ("time" %in% names(list(...))) {
    stop(
      "The first argument is named `x` for S3 dispatch, not `time`. ",
      "Call positionally -- decompose_cindex(time, status, risk) -- ",
      "or use x = .",
      call. = FALSE
    )
  }
  time <- x
  validate_survival_inputs(time, status, risk)
  validate_weights_and_orientation(weights, higher_is_riskier)
  if (!is.numeric(n_boot) || length(n_boot) != 1L || n_boot < 0) {
    stop("`n_boot` must be a single non-negative number.", call. = FALSE)
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }
  n_boot <- as.integer(n_boot)

  # Orientation is applied once, before any counting, so C_ee, C_ec and
  # C_global stay mutually consistent and the identity survives.
  if (!higher_is_riskier) risk <- -risk

  pc <- pair_counts(time, status, risk, weights)
  dec <- decomp_from_pairs(pc)
  C_ee <- dec$C_ee
  C_ec <- dec$C_ec
  C_global <- dec$C_global

  boot <- if (n_boot > 0L) {
    bootstrap_decomp(time, status, risk, weights, n_boot)
  } else {
    NULL
  }
  n_boot_valid <- if (is.null(boot)) NA_integer_ else attr(boot, "n_valid")

  structure(
    list(
      C_ee = C_ee,
      C_ec = C_ec,
      C_global = C_global,
      gap = C_ec - C_ee,
      W_ee = pc$W_ee,
      W_ec = pc$W_ec,
      N_ee = pc$N_ee,
      N_ec = pc$N_ec,
      n = length(time),
      n_events = sum(status == 1),
      censoring_rate = mean(status == 0),
      weighting = weights$name,
      higher_is_riskier = higher_is_riskier,
      boot = boot,
      n_boot = n_boot,
      n_boot_valid = n_boot_valid,
      conf_level = conf_level
    ),
    class = "cindex_decomp"
  )
}

#' @rdname decompose_cindex
#' @export
decompose_cindex.formula <- function(x, data = parent.frame(),
                                     weights = weights_harrell(),
                                     higher_is_riskier = TRUE,
                                     n_boot = 0, conf_level = 0.95, ...) {
  mf <- stats::model.frame(x, data = data, na.action = stats::na.fail)
  resp <- stats::model.response(mf)
  if (!inherits(resp, "Surv")) {
    stop("The left-hand side of the formula must be a Surv() object.",
         call. = FALSE)
  }
  if (ncol(mf) != 2L) {
    stop("The right-hand side must name exactly one risk score variable.",
         call. = FALSE)
  }
  decompose_cindex.default(
    x = as.numeric(resp[, 1L]),
    status = as.numeric(resp[, 2L]),
    risk = as.numeric(mf[[2L]]),
    weights = weights,
    higher_is_riskier = higher_is_riskier,
    n_boot = n_boot,
    conf_level = conf_level,
    ...
  )
}
