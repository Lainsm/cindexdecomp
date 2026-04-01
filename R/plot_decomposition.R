#' Plot C-Index Decomposition as a Dumbbell Chart
#'
#' @description
#' Visualises the decomposed C-Index across multiple models as a
#' dumbbell plot, showing Event-Event (C_ee), Event-Censored (C_ec)
#' and Global C-Index side by side with uncertainty bands.
#'
#' @param results_df A data frame with columns: model, ci_ee, ci_ec,
#'   global_c, sd_ee, sd_ec, sd_global. Typically output from
#'   \code{monte_carlo_decompose()}.
#' @param col_ee Character. Colour for Event-Event points.
#'   Default is IBM magenta "#DC267F".
#' @param col_ec Character. Colour for Event-Censored points.
#'   Default is IBM blue "#648FFF".
#' @param background Character. Plot background colour.
#'   Default is Dracula "#282a36".
#'
#' @return A ggplot2 object
#'
#' @examples
#' # Using example results
#' results <- data.frame(
#'   model    = c("Cox PH", "XGBoost"),
#'   ci_ee    = c(0.585, 0.576),
#'   ci_ec    = c(0.664, 0.659),
#'   global_c = c(0.641, 0.635),
#'   sd_ee    = c(0.022, 0.024),
#'   sd_ec    = c(0.022, 0.024),
#'   sd_global = c(0.018, 0.019)
#' )
#' plot_decomposition(results)
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot aes geom_segment geom_vline geom_rect
#'   geom_point geom_text annotate scale_colour_manual
#'   scale_x_continuous labs theme_classic theme element_rect
#'   element_text element_line element_blank margin unit
#' @importFrom dplyr select mutate
#' @export
plot_decomposition <- function(results_df,
                                col_ee     = "#DC267F",
                                col_ec     = "#648FFF",
                                background = "#282a36") {

  # Ensure model is a factor
  fig_data <- results_df
  fig_data$model <- factor(fig_data$model, 
                          levels = rev(unique(fig_data$model)))

  ggplot2::ggplot(fig_data) +

    # Background track
    ggplot2::geom_segment(
      ggplot2::aes(x = 0.48, xend = 0.73,
                   y = model, yend = model),
      colour = "#3a3c4e", linewidth = 10) +

    # Connecting line
    ggplot2::geom_segment(
      ggplot2::aes(x = ci_ee, xend = ci_ec,
                   y = model, yend = model),
      colour = "#6272a4", linewidth = 2) +

    # Baseline
    ggplot2::geom_vline(xintercept = 0.5,
                        linetype  = "dashed",
                        colour    = "#f8f8f2",
                        linewidth = 1) +

    # Uncertainty bands
    ggplot2::geom_rect(
      ggplot2::aes(xmin = ci_ee - sd_ee,
                   xmax = ci_ee + sd_ee,
                   ymin = as.numeric(model) - 0.18,
                   ymax = as.numeric(model) + 0.18),
      fill = col_ee, alpha = 0.4) +

    ggplot2::geom_rect(
      ggplot2::aes(xmin = ci_ec - sd_ec,
                   xmax = ci_ec + sd_ec,
                   ymin = as.numeric(model) - 0.18,
                   ymax = as.numeric(model) + 0.18),
      fill = col_ec, alpha = 0.4) +

    ggplot2::geom_rect(
      ggplot2::aes(xmin = global_c - sd_global,
                   xmax = global_c + sd_global,
                   ymin = as.numeric(model) - 0.18,
                   ymax = as.numeric(model) + 0.18),
      fill = "#f8f8f2", alpha = 0.15) +

    # Points
    ggplot2::geom_point(
      ggplot2::aes(x = ci_ee, y = model,
                   colour = "Event-Event (C_ee)"), size = 7) +
    ggplot2::geom_point(
      ggplot2::aes(x = ci_ec, y = model,
                   colour = "Event-Censored (C_ec)"), size = 7) +
    ggplot2::geom_point(
      ggplot2::aes(x = global_c, y = model,
                   colour = "Global C-Index"),
      size = 6, shape = 18) +

    # Value labels
    ggplot2::geom_text(
      ggplot2::aes(x = ci_ee, y = model,
                   label = sprintf("%.3f", ci_ee)),
      vjust = -1.8, size = 5.5, fontface = "bold",
      colour = "#f8f8f2", family = "Arial") +
    ggplot2::geom_text(
      ggplot2::aes(x = ci_ec, y = model,
                   label = sprintf("%.3f", ci_ec)),
      vjust = -1.8, size = 5.5, fontface = "bold",
      colour = "#f8f8f2", family = "Arial") +
    ggplot2::geom_text(
      ggplot2::aes(x = global_c, y = model,
                   label = sprintf("%.3f", global_c)),
      vjust = 2.5, size = 5.5, fontface = "bold",
      colour = "#f8f8f2", family = "Arial") +

    # Scales
    ggplot2::scale_colour_manual(
      name   = NULL,
      values = c(
        "Event-Event (C_ee)"    = col_ee,
        "Event-Censored (C_ec)" = col_ec,
        "Global C-Index"        = "#f8f8f2"
      ),
      breaks = c("Event-Event (C_ee)",
                 "Global C-Index",
                 "Event-Censored (C_ec)")
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0.48, 0.74),
      breaks = seq(0.5, 0.7, by = 0.05)
    ) +
    ggplot2::labs(x = "Concordance Index", y = NULL) +

    # Theme
    ggplot2::theme_classic(base_size = 24) +
    ggplot2::theme(
      plot.background  = ggplot2::element_rect(fill = background,
                                                colour = NA),
      panel.background = ggplot2::element_rect(fill = background,
                                                colour = NA),
      text             = ggplot2::element_text(family   = "Arial",
                                               colour   = "#f8f8f2"),
      axis.text        = ggplot2::element_text(colour   = "#f8f8f2"),
      axis.text.y      = ggplot2::element_text(size     = 22,
                                               face     = "bold",
                                               colour   = "#f8f8f2",
                                               family   = "Arial"),
      axis.title.x     = ggplot2::element_text(size     = 24,
                                               face     = "bold",
                                               colour   = "#f8f8f2",
                                               margin   = ggplot2::margin(t = 15),
                                               family   = "Arial"),
      axis.text.x      = ggplot2::element_blank(),
      axis.ticks.x     = ggplot2::element_blank(),
      axis.line        = ggplot2::element_line(colour   = "#f8f8f2"),
      axis.ticks       = ggplot2::element_line(colour   = "#f8f8f2"),
      axis.line.y      = ggplot2::element_blank(),
      axis.ticks.y     = ggplot2::element_blank(),
      legend.position  = "right",
      legend.text      = ggplot2::element_text(size     = 22,
                                               face     = "bold",
                                               colour   = "#f8f8f2",
                                               family   = "Arial"),
      legend.background = ggplot2::element_rect(fill    = background,
                                                colour  = NA),
      legend.key        = ggplot2::element_rect(fill    = background,
                                                colour  = NA),
      panel.grid.major.x = ggplot2::element_line(colour = "#3a3c4e",
                                                  linewidth = 0.4),
      panel.grid.minor   = ggplot2::element_blank(),
      plot.margin        = ggplot2::margin(20, 20, 20, 20)
    )
}