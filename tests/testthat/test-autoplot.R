skip_if_no_ggplot <- function() skip_if_not_installed("ggplot2")

test_that("autoplot dispatches on a decomposition", {
  skip_if_no_ggplot()
  d <- make_test_data(200, seed = 1)
  p <- ggplot2::autoplot(decompose_cindex(d$time, d$status, d$risk))
  expect_s3_class(p, "ggplot")
})

test_that("a single decomposition's plot doesn't label its row with the literal placeholder", {
  # Regression: autoplot.cindex_decomp() built a synthetic one-row table
  # with model = "model" (an internal placeholder, not a real name), so a
  # single decomposition rendered with the literal word "model" as its
  # y-axis label.
  skip_if_no_ggplot()
  d <- make_test_data(200, seed = 1)
  p <- ggplot2::autoplot(decompose_cindex(d$time, d$status, d$risk))
  expect_true(inherits(p$theme$axis.text.y, "element_blank"))
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

test_that("the dumbbell plot draws error bars per-model, not all-or-nothing", {
  # Regression: has_ci was `all(is.finite(...))` across every row of the
  # comparison table, so one model with an unstable (NA) bootstrap CI
  # suppressed error bars for EVERY model, not just the affected one.
  skip_if_no_ggplot()
  tab <- data.frame(
    model = c("stable", "unstable"),
    ci_ee = c(0.60, 0.55), ci_ec = c(0.70, 0.50), global_c = c(0.65, 0.52),
    ci_ee_lo = c(0.55, NA), ci_ee_hi = c(0.65, NA),
    ci_ec_lo = c(0.65, NA), ci_ec_hi = c(0.75, NA),
    stringsAsFactors = FALSE
  )
  p <- dumbbell_plot(tab, "Harrell", dark = FALSE, conf_level = 0.95)
  built <- ggplot2::ggplot_build(p)
  is_errorbar <- vapply(p$layers, function(l) inherits(l$geom, "GeomErrorbar"),
                        logical(1))
  expect_true(any(is_errorbar))
  for (dd in built$data[is_errorbar]) expect_equal(nrow(dd), 1)
  expect_true(grepl("CI", p$labels$subtitle, fixed = TRUE))
})

test_that("the curve plot's ribbon never inverts when the gap is negative", {
  # Regression: geom_ribbon(ymin = C_ee, ymax = C_global) assumed
  # C_ee <= C_global, but C_global is a weighted average of C_ee and C_ec,
  # so a threshold with a negative gap (C_ec < C_ee) makes C_global < C_ee
  # and the ribbon draw inverted.
  skip_if_no_ggplot()
  cv <- structure(
    list(
      data = data.frame(
        threshold = c(10, 20), censoring = c(0.3, 0.6),
        C_ee = c(0.70, 0.60), C_ec = c(0.65, 0.75),
        C_global = c(0.68, 0.68),  # row 1: C_ee > C_global, negative gap
        low_precision = c(FALSE, FALSE)
      ),
      weighting = "Harrell"
    ),
    class = "cindex_curve"
  )
  p <- ggplot2::autoplot(cv)
  built <- ggplot2::ggplot_build(p)
  is_ribbon <- vapply(p$layers, function(l) inherits(l$geom, "GeomRibbon"),
                      logical(1))
  ribbon_data <- built$data[[which(is_ribbon)]]
  expect_true(all(ribbon_data$ymin <= ribbon_data$ymax))
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
