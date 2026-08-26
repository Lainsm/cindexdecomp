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
  r <- decompose_cindex(d$time, d$status, d$risk, n_boot = 100)
  expect_equal(
    r$boot[, "gap"],
    r$boot[, "C_ec"] - r$boot[, "C_ee"],
    tolerance = 1e-12
  )
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
