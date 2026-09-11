test_that("bootstrap_decomp returns the expected shape", {
  d <- make_test_data(150, seed = 1)
  set.seed(1)
  b <- bootstrap_decomp(
    d$time, d$status, d$risk, weights_harrell(),
    n_boot = 20
  )
  expect_true(is.matrix(b))
  expect_equal(nrow(b), 20)
  expect_equal(colnames(b), c("C_ee", "C_ec", "C_global", "gap"))
})

test_that("the bootstrap resamples SUBJECTS, not pairs", {
  # Regression guard. 150 subjects generate tens of thousands of pairs.
  # Treating pairs as independent understates SE(C_ee) by roughly 6x, so a
  # correct subject-level bootstrap must be far wider than the binomial
  # pair-based SE.
  d <- make_test_data(150, seed = 2)
  r <- decompose_cindex(d$time, d$status, d$risk, n_boot = 300)
  se_boot <- sd(r$boot[, "C_ee"], na.rm = TRUE)
  se_pairs <- sqrt(r$C_ee * (1 - r$C_ee) / r$N_ee)
  expect_gt(se_boot, 3 * se_pairs)
})

test_that("the gap is computed within each replicate, preserving correlation", {
  d <- make_test_data(150, seed = 3)
  r <- decompose_cindex(d$time, d$status, d$risk, n_boot = 300)

  # Consistency: the stored gap column matches its definition.
  expect_equal(
    r$boot[, "gap"],
    r$boot[, "C_ec"] - r$boot[, "C_ee"],
    tolerance = 1e-12
  )

  # Discrimination: C_ee and C_ec are computed from the SAME resampled
  # subjects, so they co-vary. An implementation that resampled them
  # independently would drive this correlation to ~0.
  expect_gt(cor(r$boot[, "C_ee"], r$boot[, "C_ec"], use = "complete.obs"), 0.1)

  # Consequence: because they co-vary positively, the paired bootstrap gives
  # a NARROWER gap SE than assuming independence would. Deriving the gap by
  # combining two marginal SEs overstates it (measured ~26% on n=300).
  sd_gap_paired <- stats::sd(r$boot[, "gap"], na.rm = TRUE)
  sd_gap_naive <- sqrt(stats::sd(r$boot[, "C_ee"], na.rm = TRUE)^2 +
    stats::sd(r$boot[, "C_ec"], na.rm = TRUE)^2)
  expect_lt(sd_gap_paired, sd_gap_naive)
})

test_that("the bootstrap is reproducible under set.seed", {
  d <- make_test_data(120, seed = 4)
  set.seed(99)
  a <- decompose_cindex(d$time, d$status, d$risk, n_boot = 40)
  set.seed(99)
  b <- decompose_cindex(d$time, d$status, d$risk, n_boot = 40)
  expect_equal(a$boot, b$boot)
})

test_that("n_boot = 0 stores no replicates", {
  d <- make_test_data(100, seed = 5)
  r <- decompose_cindex(d$time, d$status, d$risk)
  expect_null(r$boot)
})

test_that("confint errors clearly when there are no replicates", {
  d <- make_test_data(100, seed = 6)
  r <- decompose_cindex(d$time, d$status, d$risk)
  expect_error(confint(r), "n_boot")
})

test_that("confint returns intervals containing the point estimates", {
  d <- make_test_data(200, seed = 7)
  r <- decompose_cindex(d$time, d$status, d$risk, n_boot = 200)
  ci <- confint(r)
  expect_equal(rownames(ci), c("C_ee", "C_ec", "C_global", "gap"))
  expect_lte(ci["C_ee", 1], r$C_ee)
  expect_gte(ci["C_ee", 2], r$C_ee)
  expect_lte(ci["gap", 1], r$gap)
  expect_gte(ci["gap", 2], r$gap)
})

test_that("confint validates level the same way decompose_cindex validates conf_level", {
  # Regression: confint.cindex_decomp() used `level` raw. level = -0.2
  # produced probs c(0.6, 0.4) -- an inverted interval, lower bound above
  # upper -- and level = 0 produced a zero-width interval, both silently.
  # decompose_cindex() already rejects conf_level outside (0, 1); confint()
  # must match for `level`.
  d <- make_test_data(150, seed = 22)
  r <- decompose_cindex(d$time, d$status, d$risk, n_boot = 100)
  expect_error(confint(r, level = -0.2), "between 0 and 1")
  expect_error(confint(r, level = 0), "between 0 and 1")
  expect_error(confint(r, level = 1), "between 0 and 1")
  expect_error(confint(r, level = c(0.9, 0.95)), "between 0 and 1")
})

test_that("confint respects the level argument", {
  d <- make_test_data(200, seed = 8)
  r <- decompose_cindex(d$time, d$status, d$risk, n_boot = 300)
  wide <- confint(r, level = 0.99)
  narrow <- confint(r, level = 0.80)
  expect_gt(
    wide["C_ee", 2] - wide["C_ee", 1],
    narrow["C_ee", 2] - narrow["C_ee", 1]
  )
})

test_that("confint can select a subset of parameters", {
  d <- make_test_data(150, seed = 9)
  r <- decompose_cindex(d$time, d$status, d$risk, n_boot = 100)
  ci <- confint(r, parm = "gap")
  expect_equal(rownames(ci), "gap")
})

test_that("n_boot is validated", {
  d <- make_test_data(100, seed = 10)
  expect_error(
    decompose_cindex(d$time, d$status, d$risk, n_boot = -5),
    "non-negative"
  )
})

test_that("n_boot = Inf gives the clean validation message, not a coercion warning", {
  # Regression: n_boot = Inf passed the `n_boot < 0` check (Inf is not < 0),
  # then as.integer(Inf) silently produced NA with a "coercion" warning, and
  # `if (n_boot > 0L)` on that NA threw "missing value where TRUE/FALSE
  # needed" -- a confusing internal error instead of the package's own
  # message.
  d <- make_test_data(100, seed = 10)
  expect_error(
    decompose_cindex(d$time, d$status, d$risk, n_boot = Inf),
    "non-negative"
  )
})

test_that("confint warns when many bootstrap replicates are dropped", {
  # Few subjects and heavy censoring: most subject resamples will lack
  # enough event-event pairs, so W_ee <= 0 for a large share of replicates.
  d <- make_test_data(10, censor_rate = 0.9, seed = 11)
  r <- decompose_cindex(d$time, d$status, d$risk, n_boot = 200)
  expect_lt(r$n_boot_valid, 0.9 * r$n_boot)
  expect_warning(confint(r), "usable")
})

test_that("confint does not warn on a healthy fixture", {
  d <- make_test_data(150, seed = 21)
  r <- decompose_cindex(d$time, d$status, d$risk, n_boot = 200)
  expect_equal(r$n_boot_valid, r$n_boot)
  expect_no_warning(confint(r))
})

test_that("bootstrap_decomp keeps a usable C_ee when only C_ec's side is empty", {
  # Regression: bootstrap_decomp() used to discard the ENTIRE replicate (all
  # four columns NA'd together) whenever EITHER side had zero weight, even
  # though the point estimate computes C_ee and C_ec independently. Few
  # subjects and LOW censoring makes it common for a resample to miss the
  # handful of censored subjects entirely, giving W_ec == 0 while W_ee > 0.
  d <- make_test_data(10, censor_rate = 0.15, seed = 4)
  set.seed(1)
  b <- bootstrap_decomp(
    d$time, d$status, d$risk, weights_harrell(),
    n_boot = 300
  )
  n_ee_only <- sum(!is.na(b[, "C_ee"]) & is.na(b[, "C_ec"]))
  expect_gt(n_ee_only, 0)
})

test_that("confint computes bounds via the shared boot_percentile helper", {
  # Regression guard for Ruling R12: confint.cindex_decomp() used to
  # hand-duplicate the quantile formula instead of calling boot_percentile(),
  # the same helper compare_decompositions() and the dumbbell autoplot()
  # paths use. Mocking boot_percentile() to return a distinguishable sentinel
  # proves confint() actually routes through it, not just that the two
  # formulas happen to agree numerically.
  d <- make_test_data(150, seed = 30)
  r <- decompose_cindex(d$time, d$status, d$risk, n_boot = 100)
  testthat::local_mocked_bindings(
    boot_percentile = function(boot, col, level) c(-99, 99)
  )
  ci <- confint(r)
  expect_true(all(ci[, 1] == -99))
  expect_true(all(ci[, 2] == 99))
})

test_that("the formula method forwards n_boot and conf_level", {
  d <- make_test_data(120, seed = 12)
  df <- data.frame(time = d$time, status = d$status, risk = d$risk)
  r <- decompose_cindex(survival::Surv(time, status) ~ risk,
    data = df,
    n_boot = 25, conf_level = 0.9
  )
  expect_equal(nrow(r$boot), 25)
  expect_equal(r$n_boot, 25)
  expect_equal(r$conf_level, 0.9)
  expect_equal(rownames(confint(r)), c("C_ee", "C_ec", "C_global", "gap"))
})
