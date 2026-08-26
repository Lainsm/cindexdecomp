#' Plot a concordance decomposition
#'
#' @description
#' `autoplot()` methods for the three result classes. Because the object
#' carries its own class, the right plot is selected automatically and a
#' result can never be paired with the wrong chart.
#'
#' The dumbbell plots (`cindex_decomp`, `cindex_comparison`) draw error
#' bars from the object's own bootstrap replicates as a percentile interval
#' at its `conf_level` (e.g. 95%), not a symmetric +/-1 SD band; the level
#' actually drawn is named in the subtitle. Bars are omitted when no
#' bootstrap was run (`n_boot = 0`).
#'
#' @param object A `cindex_decomp`, `cindex_curve` or `cindex_comparison`
#'   object.
#' @param dark Logical. Passed to [theme_cindex()].
#' @param ... Ignored.
#' @return A `ggplot` object.
#' @name autoplot.cindexdecomp
NULL

dumbbell_plot <- function(tab, weighting, dark, conf_level = NULL) {
  p <- cindex_palette(dark)
  tab$model <- factor(tab$model, levels = rev(unique(tab$model)))
  ci_cols <- c("ci_ee_lo", "ci_ee_hi", "ci_ec_lo", "ci_ec_hi")
  has_ci <- all(ci_cols %in% names(tab)) &&
    all(vapply(tab[ci_cols], function(col) all(is.finite(col)), logical(1)))

  gg <- ggplot2::ggplot(tab) +
    ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed",
                        colour = p$fg, linewidth = 0.6) +
    ggplot2::geom_segment(
      ggplot2::aes(x = .data$ci_ee, xend = .data$ci_ec,
                   y = .data$model, yend = .data$model),
      colour = p$grid, linewidth = 2
    )

  if (has_ci) {
    gg <- gg +
      ggplot2::geom_errorbar(
        ggplot2::aes(xmin = .data$ci_ee_lo, xmax = .data$ci_ee_hi,
                     y = .data$model),
        orientation = "y", width = 0.12, colour = p$ee, linewidth = 0.7
      ) +
      ggplot2::geom_errorbar(
        ggplot2::aes(xmin = .data$ci_ec_lo, xmax = .data$ci_ec_hi,
                     y = .data$model),
        orientation = "y", width = 0.12, colour = p$ec, linewidth = 0.7
      )
  }

  ci_txt <- if (has_ci && !is.null(conf_level) && is.finite(conf_level)) {
    sprintf("error bars: %s%% CI", format(100 * conf_level, trim = TRUE))
  } else {
    NULL
  }
  subtitle <- paste(
    c(paste0("Weighting: ", weighting), ci_txt,
      "dashed line marks chance (0.50)"),
    collapse = " \u00b7 "
  )

  gg +
    ggplot2::geom_point(
      ggplot2::aes(x = .data$ci_ee, y = .data$model,
                   colour = "Event-Event"), size = 4.5
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data$ci_ec, y = .data$model,
                   colour = "Event-Censored"), size = 4.5
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data$global_c, y = .data$model,
                   colour = "Global"), size = 3.5, shape = 18
    ) +
    ggplot2::scale_colour_manual(
      name = NULL,
      values = c("Event-Event" = p$ee, "Event-Censored" = p$ec,
                 "Global" = p$global),
      breaks = c("Event-Event", "Global", "Event-Censored")
    ) +
    ggplot2::labs(
      x = "Concordance index", y = NULL,
      title = "C-index decomposition",
      subtitle = subtitle
    ) +
    theme_cindex(dark = dark)
}

#' @rdname autoplot.cindexdecomp
#' @export
autoplot.cindex_comparison <- function(object, dark = FALSE, ...) {
  dumbbell_plot(object$table, object$weighting, dark,
               conf_level = object$conf_level)
}

#' @rdname autoplot.cindexdecomp
#' @export
autoplot.cindex_decomp <- function(object, dark = FALSE, ...) {
  ee_ci <- boot_percentile(object$boot, "C_ee", object$conf_level)
  ec_ci <- boot_percentile(object$boot, "C_ec", object$conf_level)
  tab <- data.frame(
    model = "model",
    ci_ee = object$C_ee,
    ci_ec = object$C_ec,
    global_c = object$C_global,
    ci_ee_lo = ee_ci[1], ci_ee_hi = ee_ci[2],
    ci_ec_lo = ec_ci[1], ci_ec_hi = ec_ci[2],
    stringsAsFactors = FALSE
  )
  dumbbell_plot(tab, object$weighting, dark, conf_level = object$conf_level)
}

#' @rdname autoplot.cindexdecomp
#' @export
autoplot.cindex_curve <- function(object, dark = FALSE, ...) {
  p <- cindex_palette(dark)
  dat <- object$data
  dat$precision <- ifelse(dat$low_precision, "low", "ok")

  ggplot2::ggplot(dat, ggplot2::aes(x = .data$censoring)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$C_ee, ymax = .data$C_global),
      fill = p$ee, alpha = 0.10, na.rm = TRUE
    ) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed",
                        colour = p$fg, linewidth = 0.6) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$C_ec, colour = "Event-Censored"),
      linewidth = 1, na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$C_global, colour = "Global"),
      linewidth = 1, na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$C_ee, colour = "Event-Event"),
      linewidth = 1, na.rm = TRUE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = .data$C_ec, colour = "Event-Censored",
                   alpha = .data$precision), size = 2.4, na.rm = TRUE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = .data$C_global, colour = "Global",
                   alpha = .data$precision), size = 2.4, na.rm = TRUE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = .data$C_ee, colour = "Event-Event",
                   alpha = .data$precision), size = 2.4, na.rm = TRUE
    ) +
    # Low-precision points are drawn faintly, never removed.
    ggplot2::scale_alpha_manual(
      values = c(ok = 1, low = 0.35), guide = "none"
    ) +
    ggplot2::scale_colour_manual(
      name = NULL,
      values = c("Event-Event" = p$ee, "Event-Censored" = p$ec,
                 "Global" = p$global),
      breaks = c("Event-Event", "Global", "Event-Censored")
    ) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(
      x = "Cohort censoring rate", y = "Concordance index",
      title = "Concordance under increasing censoring",
      subtitle = paste0("Weighting: ", object$weighting,
                        " \u00b7 shaded band spans C_ee to the global C-index")
    ) +
    theme_cindex(dark = dark)
}
