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
  if (!is.numeric(time) || !is.numeric(risk)) {
    stop("`time` and `risk` must be numeric.", call. = FALSE)
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
