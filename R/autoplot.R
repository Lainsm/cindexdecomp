#' Plot a concordance decomposition
#'
#' @description
#' `autoplot()` methods for the three result classes. Because the object
#' carries its own class, the right plot is selected automatically and a
#' result can never be paired with the wrong chart.
#'
#' @param object A `cindex_decomp`, `cindex_curve` or `cindex_comparison`
#'   object.
#' @param dark Logical. Passed to [theme_cindex()].
#' @param ... Ignored.
#' @return A `ggplot` object.
#' @name autoplot.cindexdecomp
NULL

dumbbell_plot <- function(tab, weighting, dark) {
  p <- cindex_palette(dark)
  tab$model <- factor(tab$model, levels = rev(unique(tab$model)))
  has_sd <- all(is.finite(tab$sd_ee))

  gg <- ggplot2::ggplot(tab) +
    ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed",
                        colour = p$fg, linewidth = 0.6) +
    ggplot2::geom_segment(
      ggplot2::aes(x = .data$ci_ee, xend = .data$ci_ec,
                   y = .data$model, yend = .data$model),
      colour = p$grid, linewidth = 2
    )

  if (has_sd) {
    gg <- gg +
      ggplot2::geom_errorbar(
        ggplot2::aes(xmin = .data$ci_ee - .data$sd_ee,
                     xmax = .data$ci_ee + .data$sd_ee,
                     y = .data$model),
        orientation = "y", width = 0.12, colour = p$ee, linewidth = 0.7
      ) +
      ggplot2::geom_errorbar(
        ggplot2::aes(xmin = .data$ci_ec - .data$sd_ec,
                     xmax = .data$ci_ec + .data$sd_ec,
                     y = .data$model),
        orientation = "y", width = 0.12, colour = p$ec, linewidth = 0.7
      )
  }

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
      subtitle = paste0("Weighting: ", weighting,
                        " · dashed line marks chance (0.50)")
    ) +
    theme_cindex(dark = dark)
}

#' @rdname autoplot.cindexdecomp
#' @export
autoplot.cindex_comparison <- function(object, dark = FALSE, ...) {
  dumbbell_plot(object$table, object$weighting, dark)
}

#' @rdname autoplot.cindexdecomp
#' @export
autoplot.cindex_decomp <- function(object, dark = FALSE, ...) {
  sd_of <- function(col) {
    if (is.null(object$boot)) NA_real_ else stats::sd(object$boot[, col],
                                                      na.rm = TRUE)
  }
  tab <- data.frame(
    model = "model",
    ci_ee = object$C_ee,
    ci_ec = object$C_ec,
    global_c = object$C_global,
    sd_ee = sd_of("C_ee"),
    sd_ec = sd_of("C_ec"),
    sd_global = sd_of("C_global"),
    stringsAsFactors = FALSE
  )
  dumbbell_plot(tab, object$weighting, dark)
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
                        " · shaded band is the masking gap")
    ) +
    theme_cindex(dark = dark)
}
