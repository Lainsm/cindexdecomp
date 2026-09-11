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
