# Shared fixtures. `ties = TRUE` rounds times to whole numbers, which is how
# real survival data arrives (days) and is the case that exercises the tie rules.
make_test_data <- function(n = 300, censor_rate = 0.4, ties = FALSE, seed = 1) {
  set.seed(seed)
  risk <- rnorm(n)
  event_time <- rexp(n, rate = exp(0.8 * risk) * 0.1)
  # solve roughly for a censoring rate
  cens_time <- rexp(n, rate = 0.1 * censor_rate / (1 - censor_rate))
  time <- pmin(event_time, cens_time)
  status <- as.numeric(event_time <= cens_time)
  if (ties) time <- round(time)
  # guarantee at least two events so validation never trips on the fixture
  if (sum(status) < 2) {
    status[order(time)[1:2]] <- 1
  }
  list(time = time, status = status, risk = risk)
}

# Fixture with a deliberate event/censoring tie at the same time.
make_tied_pair_data <- function() {
  list(
    time   = c(100, 100, 200, 300),
    status = c(  1,   0,   1,   1),
    risk   = c(0.9, 0.2, 0.5, 0.1)
  )
}
