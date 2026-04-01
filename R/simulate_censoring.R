#' Simulate Increasing Censoring on a Survival Dataset
#'
#' @description
#' Implements the "Time Machine" simulation: freezes risk scores from
#' a fitted model and applies artificial censoring thresholds iteratively.
#' At each threshold, events occurring after threshold t are masked as
#' censored. The C-Index is decomposed at each step to expose the
#' Censoring Illusion.
#'
#' @param time Numeric vector of survival times
#' @param status Numeric vector of event indicators (1 = event, 0 = censored)
#' @param risk Numeric vector of predicted risk scores
#' @param n_thresholds Integer. Number of censoring thresholds to evaluate.
#'   Default is 20.
#' @param min_pairs Integer. Minimum number of Event-Event pairs required
#'   to include a threshold. Default is 5000.
#'
#' @return A data frame with columns:
#' \itemize{
#'   \item \code{threshold} Applied censoring threshold
#'   \item \code{censoring} Resulting censoring rate
#'   \item \code{global_c} Global C-Index
#'   \item \code{ci_ee} Event-Event C-Index
#'   \item \code{ci_ec} Event-Censored C-Index
#'   \item \code{n_ee} Number of Event-Event pairs
#'   \item \code{n_ec} Number of Event-Censored pairs
#' }
#'
#' @examples
#' set.seed(42)
#' time   <- rexp(200, rate = 0.1)
#' status <- rbinom(200, 1, 0.6)
#' risk   <- rnorm(200)
#'
#' sim <- simulate_censoring(time, status, risk, min_pairs = 50)
#' head(sim)
#'
#' @importFrom survival concordance Surv
#' @export
simulate_censoring <- function(time, status, risk,
                                n_thresholds = 20,
                                min_pairs    = NULL) {

  event_idx <- which(status == 1)
  time_e    <- time[event_idx]
  status_e  <- status[event_idx]
  risk_e    <- risk[event_idx]

  if (is.null(min_pairs)) {
    min_pairs <- max(50, round(length(event_idx)^1.5 * 0.1))
    message("min_pairs automatically set to ", min_pairs)
  }

  if (length(event_idx) < 200) {
    warning("Small dataset detected (", length(event_idx),
            " events). Results may be unstable.")
  }

  thresholds <- quantile(time_e,
                         probs = seq(0.05, 0.70,
                                     length.out = n_thresholds))

  results_list <- list()

  for (t in thresholds) {
    sim_time   <- pmin(time_e, t)
    sim_status <- ifelse(time_e <= t, 1, 0)

    if (sum(sim_status) < 5 || sum(sim_status == 0) < 5) next

    censoring_rate <- 1 - mean(sim_status)

    global_c_raw <- survival::concordance(
      survival::Surv(sim_time, sim_status) ~ risk_e
    )$concordance
    global_c <- ifelse(global_c_raw < 0.5, 
                       1 - global_c_raw, global_c_raw)

    decomp <- decompose_cindex(sim_time, sim_status, risk_e)

    ci_ee <- ifelse(decomp$N_ee < min_pairs, NA, decomp$CI_ee)
    ci_ec <- ifelse(decomp$N_ec < min_pairs | decomp$CI_ec < 0.5,
                NA, decomp$CI_ec)

    results_list[[as.character(t)]] <- data.frame(
      threshold  = t,
      censoring  = censoring_rate,
      global_c   = global_c,
      ci_ee      = ci_ee,
      ci_ec      = ci_ec,
      n_ee       = decomp$N_ee,
      n_ec       = decomp$N_ec
    )
  }

  results <- dplyr::bind_rows(results_list)
  results <- results[order(results$censoring), ]   # sort low to high
  results <- results[results$censoring <= 0.80, ]  # hard cap at 80%
  
  return(results)
}