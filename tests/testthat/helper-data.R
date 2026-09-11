# Shared fixtures. `ties = TRUE` rounds times to whole numbers, which is how
# real survival data arrives (days) and is the case that exercises the tie
# rules.
make_test_data <- function(n = 300, censor_rate = 0.4, ties = FALSE, seed = 1) {
  set.seed(seed)
  risk <- rnorm(n)
  event_time <- rexp(n, rate = exp(0.8 * risk) * 0.1)
  # solve roughly for a censoring rate
  cens_time <- rexp(n, rate = 0.1 * censor_rate / (1 - censor_rate))
  time <- pmin(event_time, cens_time)
  status <- as.numeric(event_time <= cens_time)
  # Day-scale rounding. Real survival data is recorded in whole days:
  # survival::lung has 149 unique times among 167 subjects, with a maximum
  # tie group of 3. Rounding the raw scale directly would instead pile 300
  # subjects onto ~31 distinct times with ~10% at time zero, which is not
  # a tie pattern any real dataset exhibits.
  if (ties) time <- round(time * 30) + 1
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
    status = c(1, 0, 1, 1),
    risk   = c(0.9, 0.2, 0.5, 0.1)
  )
}
