skip_if_no_ggplot <- function() skip_if_not_installed("ggplot2")

test_that("autoplot dispatches on a decomposition", {
  skip_if_no_ggplot()
  d <- make_test_data(200, seed = 1)
  p <- ggplot2::autoplot(decompose_cindex(d$time, d$status, d$risk))
  expect_s3_class(p, "ggplot")
})

test_that("autoplot dispatches on a censoring curve", {
  skip_if_no_ggplot()
  d <- make_test_data(250, seed = 2)
  p <- ggplot2::autoplot(censoring_curve(d$time, d$status, d$risk,
                                         n_thresholds = 6))
  expect_s3_class(p, "ggplot")
})

test_that("autoplot dispatches on a comparison", {
  skip_if_no_ggplot()
  d <- make_test_data(200, seed = 3)
  cmp <- compare_decompositions(
    list(a = d$risk, b = d$risk + rnorm(200)),
    d$time, d$status, n_boot = 20
  )
  expect_s3_class(ggplot2::autoplot(cmp), "ggplot")
})

test_that("no plot hardcodes axis limits", {
  # Regression: scale_x_continuous(limits = c(0.48, 0.74)) silently dropped
  # any model outside that window, including worse-than-chance models.
  skip_if_no_ggplot()
  d <- make_test_data(200, seed = 4)
  cmp <- compare_decompositions(list(a = d$risk), d$time, d$status, n_boot = 0)
  p <- ggplot2::autoplot(cmp)
  xs <- p$scales$get_scales("x")
  expect_true(is.null(xs) || is.null(xs$limits))
})

test_that("a worse-than-chance model still appears in the plot data", {
  skip_if_no_ggplot()
  d <- make_test_data(250, seed = 5)
  cmp <- compare_decompositions(list(backwards = -d$risk), d$time,
                                d$status, n_boot = 0)
  expect_lt(cmp$table$global_c, 0.5)
  p <- ggplot2::autoplot(cmp)
  built <- ggplot2::ggplot_build(p)
  xs <- unlist(lapply(built$data, function(l) l$x))
  expect_true(any(xs < 0.5, na.rm = TRUE))
})

test_that("no plot sets a hardcoded font family", {
  # Regression: family = "Arial" appeared 14 times; Arial is absent from
  # CRAN's Linux check machines and emits warnings on every plot call.
  skip_if_no_ggplot()
  d <- make_test_data(150, seed = 6)
  ps <- list(
    ggplot2::autoplot(decompose_cindex(d$time, d$status, d$risk)),
    ggplot2::autoplot(censoring_curve(d$time, d$status, d$risk,
                                      n_thresholds = 4))
  )
  for (p in ps) {
    expect_false(identical(p$theme$text$family, "Arial"))
  }
})

test_that("theme_cindex defaults to light and offers dark", {
  skip_if_no_ggplot()
  light <- theme_cindex()
  dark <- theme_cindex(dark = TRUE)
  expect_s3_class(light, "theme")
  expect_equal(light$plot.background$fill, "#FFFFFF")
  expect_equal(dark$plot.background$fill, "#282A36")
})

test_that("theme_cindex sets no font family", {
  skip_if_no_ggplot()
  expect_false(identical(theme_cindex()$text$family, "Arial"))
})

test_that("the curve plot's subtitle names what the ribbon actually spans", {
  # Regression: the ribbon is geom_ribbon(ymin = C_ee, ymax = C_global),
  # but the subtitle called it "the masking gap" -- a different quantity
  # (C_ec - C_ee) that is drawn nowhere on this plot. Measured: ribbon
  # spanned 0.054-0.431 while the true gap spanned 0.064-0.437.
  skip_if_no_ggplot()
  d <- make_test_data(250, seed = 8)
  cv <- censoring_curve(d$time, d$status, d$risk, n_thresholds = 6)
  p <- ggplot2::autoplot(cv)
  subtitle <- p$labels$subtitle
  expect_false(grepl("masking gap", subtitle))
  expect_true(grepl("C_ee", subtitle, fixed = TRUE))
})

test_that("the dumbbell plot's error bars are labelled with the level they draw", {
  # Regression: error bars were xmin/xmax = ci_ee -/+ sd_ee with nothing in
  # the legend, subtitle or docs stating they were +/-1 SD (~68% coverage),
  # readable as a 95% interval. Bars now draw the object's own conf_level
  # percentile interval and the subtitle names the level.
  skip_if_no_ggplot()
  d <- make_test_data(200, seed = 9)
  r <- decompose_cindex(d$time, d$status, d$risk, n_boot = 150,
                        conf_level = 0.9)
  p <- ggplot2::autoplot(r)
  subtitle <- p$labels$subtitle
  expect_true(grepl("90% CI", subtitle, fixed = TRUE))
  expect_false(grepl("SD", subtitle, fixed = TRUE))

  built <- ggplot2::ggplot_build(p)
  errorbar_layers <- vapply(
    p$layers, function(l) inherits(l$geom, "GeomErrorbar"), logical(1)
  )
  expect_true(any(errorbar_layers))
})

test_that("the dumbbell plot omits error bars, and any CI claim, when n_boot = 0", {
  skip_if_no_ggplot()
  d <- make_test_data(150, seed = 10)
  r <- decompose_cindex(d$time, d$status, d$risk, n_boot = 0)
  p <- ggplot2::autoplot(r)
  expect_false(grepl("CI", p$labels$subtitle, fixed = TRUE))
  errorbar_layers <- vapply(
    p$layers, function(l) inherits(l$geom, "GeomErrorbar"), logical(1)
  )
  expect_false(any(errorbar_layers))
})

test_that("the curve plot marks low-precision points without dropping them", {
  skip_if_no_ggplot()
  d <- make_test_data(250, seed = 7)
  cv <- censoring_curve(d$time, d$status, d$risk, n_thresholds = 6,
                        min_pairs = 10^9)
  expect_true(all(cv$data$low_precision))
  p <- ggplot2::autoplot(cv)
  built <- ggplot2::ggplot_build(p)
  n_points <- sum(vapply(built$data, function(l) nrow(l), numeric(1)))
  expect_gt(n_points, 0)
})
