test_that("print shows the global C-index and the masking gap", {
  d <- make_test_data(200, seed = 1)
  r <- decompose_cindex(d$time, d$status, d$risk)
  out <- capture.output(print(r))
  expect_true(any(grepl("Global", out)))
  expect_true(any(grepl("Event-Event", out)))
  expect_true(any(grepl("Event-Censored", out)))
  expect_true(any(grepl("gap", out, ignore.case = TRUE)))
})

test_that("print names the weighting", {
  d <- make_test_data(150, seed = 2)
  r <- decompose_cindex(d$time, d$status, d$risk, weights = weights_uno())
  expect_true(any(grepl("Uno", capture.output(print(r)))))
})

test_that("print warns that pair counts are not precision when no CI exists", {
  # Constraint: N_ee must never be shown alone as though it were a sample
  # size. 16,471 pairs from 300 patients is 300 units of information.
  d <- make_test_data(200, seed = 3)
  r <- decompose_cindex(d$time, d$status, d$risk)
  out <- capture.output(print(r))
  expect_true(any(grepl("not precision|n_boot", out)))
})

test_that("print shows confidence intervals when they exist", {
  d <- make_test_data(150, seed = 4)
  r <- decompose_cindex(d$time, d$status, d$risk, n_boot = 50)
  out <- capture.output(print(r))
  expect_true(any(grepl("95%|CI", out)))
})

test_that("print flags near-chance event-event concordance", {
  d <- make_test_data(300, seed = 5)
  r <- decompose_cindex(d$time, d$status, d$risk)

  r$C_ee <- 0.51 # near chance: should get the "near chance" wording
  out_near <- capture.output(print(r))
  expect_true(any(grepl("near chance", out_near)))

  r$C_ee <- 0.10 # well below chance: a different, more accurate finding
  out_below <- capture.output(print(r))
  expect_true(any(grepl("systematically reversed", out_below)))
  expect_false(any(grepl("near chance", out_below)))
})

test_that("print returns its input invisibly", {
  d <- make_test_data(100, seed = 6)
  r <- decompose_cindex(d$time, d$status, d$risk)
  expect_invisible(print(r))
})

test_that("summary carries the bootstrap settings", {
  d <- make_test_data(150, seed = 7)
  r <- decompose_cindex(d$time, d$status, d$risk, n_boot = 50)
  s <- summary(r)
  expect_s3_class(s, "summary.cindex_decomp")
  expect_true(any(grepl("50", capture.output(print(s)))))
})

test_that("as.data.frame returns one tidy row", {
  d <- make_test_data(150, seed = 8)
  r <- decompose_cindex(d$time, d$status, d$risk)
  df <- as.data.frame(r)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 1)
  expect_true(all(c(
    "C_ee", "C_ec", "C_global", "gap", "N_ee", "N_ec",
    "weighting"
  ) %in% names(df)))
  expect_equal(df$C_ee, r$C_ee, tolerance = 1e-12)
})

test_that("print and summary agree on the confidence level string", {
  # Regression guard: the table header and the summary's Bootstrap line
  # used to format the same conf_level two different ways, so a level
  # like 0.975 rendered as "97.5 CI" in one place and "98% level" in the
  # other. Both must show the same string.
  d <- make_test_data(150, seed = 9)
  r <- decompose_cindex(d$time, d$status, d$risk,
    n_boot = 50,
    conf_level = 0.975
  )
  print_out <- capture.output(print(r))
  summary_out <- capture.output(print(summary(r)))
  expect_true(any(grepl("97.5%", print_out, fixed = TRUE)))
  expect_true(any(grepl("97.5%", summary_out, fixed = TRUE)))
})

test_that("print caveats a bootstrap with few usable replicates", {
  # Regression guard: a zero- or near-zero-width interval built from only
  # a handful of usable replicates must not print as though it were
  # precise. confint()'s own warning goes to stderr, so capture.output()
  # (and any report pipeline capturing printed text) would otherwise show
  # the interval with no caveat at all.
  d <- make_test_data(10, censor_rate = 0.9, seed = 11)
  r <- decompose_cindex(d$time, d$status, d$risk, n_boot = 200)
  expect_lt(r$n_boot_valid, 0.9 * r$n_boot)
  out <- suppressWarnings(capture.output(print(r)))
  expect_true(any(grepl("usable", out)))
})

test_that("print does not caveat a healthy bootstrap", {
  d <- make_test_data(150, seed = 21)
  r <- decompose_cindex(d$time, d$status, d$risk, n_boot = 200)
  expect_equal(r$n_boot_valid, r$n_boot)
  out <- capture.output(print(r))
  expect_false(any(grepl("usable", out)))
})
