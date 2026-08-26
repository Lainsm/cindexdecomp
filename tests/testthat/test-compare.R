two_models <- function(seed = 1, n = 250) {
  d <- make_test_data(n, seed = seed)
  list(
    d = d,
    risks = list(good = d$risk, noisy = d$risk + rnorm(n, sd = 2))
  )
}

test_that("compare_decompositions produces the columns the dumbbell plot needs", {
  # The old plot_decomposition() required these exact columns and named a
  # monte_carlo_decompose() function that did not exist in the package.
  m <- two_models()
  cmp <- compare_decompositions(m$risks, m$d$time, m$d$status, n_boot = 30)
  expect_s3_class(cmp, "cindex_comparison")
  expect_true(all(c("model", "ci_ee", "ci_ec", "global_c",
                    "sd_ee", "sd_ec", "sd_global") %in% names(cmp$table)))
  expect_equal(nrow(cmp$table), 2)
})

test_that("model names are carried through from the list names", {
  m <- two_models()
  cmp <- compare_decompositions(m$risks, m$d$time, m$d$status, n_boot = 0)
  expect_equal(cmp$table$model, c("good", "noisy"))
})

test_that("a data frame of risks is accepted", {
  m <- two_models()
  cmp <- compare_decompositions(as.data.frame(m$risks), m$d$time,
                                m$d$status, n_boot = 0)
  expect_equal(nrow(cmp$table), 2)
})

test_that("unnamed risks are rejected", {
  m <- two_models()
  expect_error(
    compare_decompositions(unname(m$risks), m$d$time, m$d$status, n_boot = 0),
    "named"
  )
})

test_that("standard deviations are NA when n_boot is 0", {
  m <- two_models()
  cmp <- compare_decompositions(m$risks, m$d$time, m$d$status, n_boot = 0)
  expect_true(all(is.na(cmp$table$sd_ee)))
})

test_that("standard deviations are finite when n_boot > 0", {
  m <- two_models()
  cmp <- compare_decompositions(m$risks, m$d$time, m$d$status, n_boot = 40)
  expect_true(all(is.finite(cmp$table$sd_ee)))
  expect_true(all(cmp$table$sd_ee > 0))
})

test_that("per-model values match a direct decompose_cindex call", {
  m <- two_models()
  cmp <- compare_decompositions(m$risks, m$d$time, m$d$status, n_boot = 0)
  direct <- decompose_cindex(m$d$time, m$d$status, m$risks$good)
  expect_equal(cmp$table$ci_ee[1], direct$C_ee, tolerance = 1e-12)
  expect_equal(cmp$table$global_c[1], direct$C_global, tolerance = 1e-12)
})

test_that("fits are retained for downstream use", {
  m <- two_models()
  cmp <- compare_decompositions(m$risks, m$d$time, m$d$status, n_boot = 0)
  expect_named(cmp$fits, c("good", "noisy"))
  expect_s3_class(cmp$fits$good, "cindex_decomp")
})

test_that("a risk vector of the wrong length is rejected", {
  m <- two_models()
  bad <- list(good = m$d$risk, short = m$d$risk[1:10])
  expect_error(
    compare_decompositions(bad, m$d$time, m$d$status, n_boot = 0),
    "same length"
  )
})

test_that("a matrix of risks is accepted; column names become model names", {
  m <- two_models()
  mat <- cbind(good = m$risks$good, noisy = m$risks$noisy)
  cmp <- compare_decompositions(mat, m$d$time, m$d$status, n_boot = 0)
  expect_equal(cmp$table$model, c("good", "noisy"))
  expect_equal(nrow(cmp$table), 2)
})

test_that("a matrix without column names hits the existing named-input error", {
  m <- two_models()
  mat <- cbind(m$risks$good, m$risks$noisy)  # no dimnames
  expect_error(
    compare_decompositions(mat, m$d$time, m$d$status, n_boot = 0),
    "named"
  )
})

test_that("the table carries W_ee/W_ec so the identity is verifiable per model", {
  m <- two_models()
  cmp <- compare_decompositions(m$risks, m$d$time, m$d$status,
                                weights = weights_uno(), n_boot = 0)
  expect_true(all(c("w_ee", "w_ec") %in% names(cmp$table)))
  rhs <- (cmp$table$w_ee * cmp$table$ci_ee +
          cmp$table$w_ec * cmp$table$ci_ec) /
    (cmp$table$w_ee + cmp$table$w_ec)
  expect_equal(cmp$table$global_c, rhs, tolerance = 1e-12)
})

test_that("the table carries conf_level percentile CI bounds alongside sd_*", {
  m <- two_models()
  cmp <- compare_decompositions(m$risks, m$d$time, m$d$status,
                                n_boot = 200, conf_level = 0.9)
  expect_true(all(c("ci_ee_lo", "ci_ee_hi", "ci_ec_lo",
                    "ci_ec_hi") %in% names(cmp$table)))
  expect_true(all(is.finite(cmp$table$ci_ee_lo)))
  expect_true(all(cmp$table$ci_ee_lo <= cmp$table$ci_ee))
  expect_true(all(cmp$table$ci_ee_hi >= cmp$table$ci_ee))
  # sd_* is retained for compatibility alongside the new bounds.
  expect_true(all(is.finite(cmp$table$sd_ee)))

  # The bounds match the 90% CI computed directly from the same replicates.
  direct <- stats::quantile(cmp$fits$good$boot[, "C_ee"],
                            probs = c(0.05, 0.95), na.rm = TRUE,
                            names = FALSE)
  expect_equal(c(cmp$table$ci_ee_lo[1], cmp$table$ci_ee_hi[1]), direct,
              tolerance = 1e-12)
})

test_that("CI bounds are NA when n_boot = 0", {
  m <- two_models()
  cmp <- compare_decompositions(m$risks, m$d$time, m$d$status, n_boot = 0)
  expect_true(all(is.na(cmp$table$ci_ee_lo)))
  expect_true(all(is.na(cmp$table$ci_ee_hi)))
})

test_that("print lists every model", {
  m <- two_models()
  cmp <- compare_decompositions(m$risks, m$d$time, m$d$status, n_boot = 0)
  out <- capture.output(print(cmp))
  expect_true(any(grepl("good", out)))
  expect_true(any(grepl("noisy", out)))
})
