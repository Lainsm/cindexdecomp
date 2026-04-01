#' Decompose the C-Index into Event-Event and Event-Censored Components
#'
#' @description
#' Decomposes Harrell's Concordance Index into two components:
#' Event-Event (C_ee) pairs where both patients experience the event,
#' and Event-Censored (C_ec) pairs where only one patient experiences
#' the event. This exposes performance masking in highly censored
#' survival data.
#'
#' @param time Numeric vector of survival times
#' @param status Numeric vector of event indicators (1 = event, 0 = censored)
#' @param risk Numeric vector of predicted risk scores (higher = higher risk)
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{CI_ee} Event-Event concordance index
#'   \item \code{CI_ec} Event-Censored concordance index
#'   \item \code{N_ee} Number of Event-Event comparable pairs
#'   \item \code{N_ec} Number of Event-Censored comparable pairs
#' }
#'
#' @examples
#' # Simulate survival data
#' set.seed(42)
#' time   <- rexp(100, rate = 0.1)
#' status <- rbinom(100, 1, 0.6)
#' risk   <- rnorm(100)
#'
#' result <- decompose_cindex(time, status, risk)
#' print(result$CI_ee)
#' print(result$CI_ec)
#'
#' @export
decompose_cindex <- function(time, status, risk) {
  
  # Input validation
  if (length(time) != length(status) || length(time) != length(risk)) {
    stop("time, status and risk must all be the same length")
  }
  if (!all(status %in% c(0, 1))) {
    stop("status must be binary (0 = censored, 1 = event)")
  }
  if (sum(status == 1) < 2) {
    stop("At least 2 events required for decomposition")
  }
  
  idx <- which(status == 1)
  ee_comparable <- 0; ee_concordant <- 0
  ec_comparable <- 0; ec_concordant <- 0

  for (i in idx) {
    j_all <- which(time > time[i])
    if (length(j_all) == 0) next
    concordant <- (risk[i] > risk[j_all]) + 
                  0.5 * (risk[i] == risk[j_all])
    j_ee <- status[j_all] == 1 
    j_ec <- status[j_all] == 0
    ee_comparable <- ee_comparable + sum(j_ee)
    ee_concordant <- ee_concordant + sum(concordant[j_ee])
    ec_comparable <- ec_comparable + sum(j_ec)
    ec_concordant <- ec_concordant + sum(concordant[j_ec])
  }

  list(
    CI_ee = ifelse(ee_comparable > 0, 
                   ee_concordant / ee_comparable, NA),
    CI_ec = ifelse(ec_comparable > 0, 
                   ec_concordant / ec_comparable, NA),
    N_ee  = ee_comparable,
    N_ec  = ec_comparable
  )
}