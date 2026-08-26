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
  r$C_ee <- 0.51  # force the condition
  expect_true(any(grepl("chance", capture.output(print(r)))))
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
  expect_true(all(c("C_ee", "C_ec", "C_global", "gap", "N_ee", "N_ec",
                    "weighting") %in% names(df)))
  expect_equal(df$C_ee, r$C_ee, tolerance = 1e-12)
})
