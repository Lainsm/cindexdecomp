#' Count weighted concordant pairs by partner status
#'
#' The single place in the package where pairs are counted. Takes plain
#' vectors plus a weighting and returns weighted totals split by whether
#' the later member of each pair had an event or was censored.
#'
#' Tie rules (chosen to match [survival::concordance()] exactly):
#' * two events at the same time are NOT comparable;
#' * an event and a censoring at the same time ARE comparable, with the
#'   event treated as occurring first;
#' * equal risk scores receive half credit.
#'
#' @param time,status,risk Numeric vectors of equal length.
#' @param weights A `cindex_weights` object.
#' @return A list with `W_ee`, `S_ee`, `W_ec`, `S_ec`, `N_ee`, `N_ec`.
#'   `W_*` and `S_*` are weighted totals. `N_ee` and `N_ec` count
#'   *contributing* pairs -- those receiving non-zero weight. Under any
#'   strictly positive weighting (Harrell's included) that equals the total
#'   number of comparable pairs; under a truncating weighting such as
#'   `weights_truncated()` or `weights_uno(tau = )`, pairs beyond `tau` are
#'   excluded from both the weighted totals and these counts, because
#'   truncation excludes a pair from the estimator rather than
#'   down-weighting it.
#' @keywords internal
#' @noRd
pair_counts <- function(time, status, risk, weights) {
  G <- censoring_km(time, status)
  event_idx <- which(status == 1)

  W_ee <- 0; S_ee <- 0; N_ee <- 0
  W_ec <- 0; S_ec <- 0; N_ec <- 0

  for (i in event_idx) {
    # An event at time[i] precedes a censoring recorded at the same time,
    # so those pairs are comparable. Two events at the same time are not.
    # `status[i] == 1`, so the second clause can never select `i` itself.
    j <- which(time > time[i] | (time == time[i] & status == 0))
    if (length(j) == 0L) next

    w <- weights$fn(time[i], G)
    if (length(w) != 1L) {
      stop(
        "Weight function must return a single number; got length ",
        length(w), " at t = ", format(time[i]), ".",
        call. = FALSE
      )
    }
    if (!is.finite(w) || w < 0) {
      stop(
        "Weight function returned a negative or non-finite value at t = ",
        format(time[i]), ".",
        call. = FALSE
      )
    }
    if (w == 0) next

    conc <- (risk[i] > risk[j]) + 0.5 * (risk[i] == risk[j])
    is_ee <- status[j] == 1

    n_ee_i <- sum(is_ee)
    n_ec_i <- length(j) - n_ee_i

    N_ee <- N_ee + n_ee_i
    N_ec <- N_ec + n_ec_i
    W_ee <- W_ee + w * n_ee_i
    W_ec <- W_ec + w * n_ec_i
    S_ee <- S_ee + w * sum(conc[is_ee])
    S_ec <- S_ec + w * sum(conc[!is_ee])
  }

  list(
    W_ee = W_ee, S_ee = S_ee, N_ee = N_ee,
    W_ec = W_ec, S_ec = S_ec, N_ec = N_ec
  )
}

#' Derive C_ee, C_ec, C_global and the gap from a pair_counts() result
#'
#' Shared by [decompose_cindex()], [censoring_curve()] and
#' `bootstrap_decomp()` so the three don't each carry their own copy of the
#' zero-weight convention. Each quantity is NA'd independently: a side with
#' no comparable pairs makes only that side's `C_*` unusable, not the other
#' side or `C_global` (which reduces to the surviving side alone when one
#' side is empty).
#'
#' @param pc A list with `W_ee`, `S_ee`, `W_ec`, `S_ec`, as returned by
#'   `pair_counts()`.
#' @return A list with `C_ee`, `C_ec`, `C_global`, `gap`.
#' @keywords internal
#' @noRd
decomp_from_pairs <- function(pc) {
  C_ee <- if (pc$W_ee > 0) pc$S_ee / pc$W_ee else NA_real_
  C_ec <- if (pc$W_ec > 0) pc$S_ec / pc$W_ec else NA_real_
  total_w <- pc$W_ee + pc$W_ec
  C_global <- if (total_w > 0) (pc$S_ee + pc$S_ec) / total_w else NA_real_
  list(C_ee = C_ee, C_ec = C_ec, C_global = C_global, gap = C_ec - C_ee)
}
