test_that("bootstrap_decomp returns the expected shape", {
  d <- make_test_data(150, seed = 1)
  set.seed(1)
  b <- bootstrap_decomp(d$time, d$status, d$risk, weights_harrell(), n_boot = 20)
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
  set.seed(99); a <- decompose_cindex(d$time, d$status, d$risk, n_boot = 40)
  set.seed(99); b <- decompose_cindex(d$time, d$status, d$risk, n_boot = 40)
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

test_that("confint respects the level argument", {
  d <- make_test_data(200, seed = 8)
  r <- decompose_cindex(d$time, d$status, d$risk, n_boot = 300)
  wide <- confint(r, level = 0.99)
  narrow <- confint(r, level = 0.80)
  expect_gt(wide["C_ee", 2] - wide["C_ee", 1],
            narrow["C_ee", 2] - narrow["C_ee", 1])
})

test_that("confint can select a subset of parameters", {
  d <- make_test_data(150, seed = 9)
  r <- decompose_cindex(d$time, d$status, d$risk, n_boot = 100)
  ci <- confint(r, parm = "gap")
  expect_equal(rownames(ci), "gap")
})

test_that("n_boot is validated", {
  d <- make_test_data(100, seed = 10)
  expect_error(decompose_cindex(d$time, d$status, d$risk, n_boot = -5),
               "non-negative")
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

test_that("the formula method forwards n_boot and conf_level", {
  d <- make_test_data(120, seed = 12)
  df <- data.frame(time = d$time, status = d$status, risk = d$risk)
  r <- decompose_cindex(survival::Surv(time, status) ~ risk, data = df,
                        n_boot = 25, conf_level = 0.9)
  expect_equal(nrow(r$boot), 25)
  expect_equal(r$n_boot, 25)
  expect_equal(r$conf_level, 0.9)
  expect_equal(rownames(confint(r)), c("C_ee", "C_ec", "C_global", "gap"))
})
