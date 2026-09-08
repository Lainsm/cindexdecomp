#' Validate survival inputs
#'
#' @param time Numeric vector of observed times.
#' @param status Numeric vector of event indicators (1 = event, 0 = censored).
#' @param risk Numeric vector of risk scores.
#' @return `invisible(TRUE)` if valid; otherwise an error is thrown.
#' @keywords internal
#' @noRd
validate_survival_inputs <- function(time, status, risk) {
  n <- length(time)
  if (length(status) != n || length(risk) != n) {
    stop("`time`, `status` and `risk` must all have the same length.", call. = FALSE)
  }
  if (n == 0L) {
    stop("`time`, `status` and `risk` must not be empty.", call. = FALSE)
  }
  if (!is.numeric(time) || !is.numeric(status) || !is.numeric(risk)) {
    stop("`time`, `status` and `risk` must be numeric.", call. = FALSE)
  }
  if (anyNA(time) || anyNA(status) || anyNA(risk)) {
    stop("`time`, `status` and `risk` must not contain NA.", call. = FALSE)
  }
  if (any(!is.finite(time)) || any(!is.finite(risk))) {
    stop("`time` and `risk` must be finite.", call. = FALSE)
  }
  if (any(time < 0)) {
    stop("`time` must be non-negative.", call. = FALSE)
  }
  if (!all(status %in% c(0, 1))) {
    stop("`status` must be binary (0 = censored, 1 = event).", call. = FALSE)
  }
  if (sum(status == 1) < 2) {
    stop("At least 2 events are required to decompose the C-index.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate a weights object and the higher_is_riskier flag
#'
#' Shared by [decompose_cindex()] and [censoring_curve()] so the two entry
#' points can't drift apart on what counts as a valid `weights`/
#' `higher_is_riskier` pair or on the wording of the resulting error.
#'
#' @param weights A `cindex_weights` object.
#' @param higher_is_riskier Logical, length 1.
#' @return `invisible(TRUE)` if valid; otherwise an error is thrown.
#' @keywords internal
#' @noRd
validate_weights_and_orientation <- function(weights, higher_is_riskier) {
  if (!inherits(weights, "cindex_weights")) {
    stop("`weights` must be a `cindex_weights` object, e.g. weights_harrell().",
         call. = FALSE)
  }
  if (!is.logical(higher_is_riskier) || length(higher_is_riskier) != 1L) {
    stop("`higher_is_riskier` must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Kaplan-Meier estimate of the censoring distribution
#'
#' Returns a step function giving `G(t)`, the probability of remaining
#' uncensored beyond `t`. Used by IPCW weightings.
#'
#' @param time Numeric vector of observed times.
#' @param status Numeric vector of event indicators.
#' @return A function of one numeric argument.
#' @keywords internal
#' @noRd
censoring_km <- function(time, status) {
  fit <- survival::survfit(survival::Surv(time, 1 - status) ~ 1)
  ftime <- fit$time
  fsurv <- fit$surv
  function(t) {
    idx <- findInterval(t, ftime)
    out <- ifelse(idx == 0L, 1, fsurv[pmax(idx, 1L)])
    # guard against division by zero in IPCW weights
    pmax(out, .Machine$double.eps)
  }
}
