# ggplot2 is a hard Import, but the plotting tests are skipped rather than
# failed where it is somehow unavailable. Lives in a helper-*.R file so
# every test file can reach it -- testthat only auto-sources these.
skip_if_no_ggplot <- function() skip_if_not_installed("ggplot2")
