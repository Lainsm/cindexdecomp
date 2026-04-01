#' Plot the Censoring Illusion Simulation
#'
#' @description
#' Visualises the result of \code{simulate_censoring()} as a line chart
#' showing how C_ee, C_ec and Global C-Index change as artificial
#' censoring increases. The shaded region between Global C-Index and
#' C_ee represents the "Censoring Illusion Zone".
#'
#' @param sim_df A data frame output from \code{simulate_censoring()}
#' @param col_ee Character. Colour for C_ee line.
#'   Default is IBM magenta "#DC267F".
#' @param col_ec Character. Colour for C_ec line.
#'   Default is IBM blue "#648FFF".
#' @param background Character. Plot background colour.
#'   Default is Dracula "#282a36".
#'
#' @return A ggplot2 object
#'
#' @examples
#' set.seed(42)
#' time   <- rexp(200, rate = 0.1)
#' status <- rbinom(200, 1, 0.6)
#' risk   <- rnorm(200)
#'
#' sim <- simulate_censoring(time, status, risk)
#' plot_simulation(sim)
#'
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_point
#'   geom_hline annotate scale_colour_manual scale_x_continuous
#'   scale_y_continuous labs theme_classic theme element_rect
#'   element_text element_line element_blank margin unit
#' @importFrom scales percent_format
#' @export
plot_simulation <- function(sim_df,
                             col_ee     = "#DC267F",
                             col_ec     = "#648FFF",
                             background = "#282a36") {

  # Get last points for end labels
  last_point <- sim_df[which.max(sim_df$censoring), ]

  ggplot2::ggplot(sim_df,
                  ggplot2::aes(x = censoring)) +

    # Illusion Zone ribbon — between global_c and ci_ee
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ci_ee, ymax = global_c),
      fill = "#ffffff", alpha = 0.08,
      na.rm = TRUE) +

    # Random chance line
    ggplot2::geom_hline(yintercept = 0.50,
                        linetype   = "dashed",
                        colour     = "#f8f8f2",
                        linewidth  = 0.8) +

    # Illusion Zone label
    ggplot2::annotate("text",
                      x        = 0.45,
                      y        = 0.57,
                      label    = "Illusion Zone",
                      colour   = "#ffffff",
                      hjust    = 0.5,
                      size     = 5,
                      fontface = "italic",
                      alpha    = 0.6,
                      family   = "Arial") +

    # Random chance label
    ggplot2::annotate("text",
                      x        = 0.02,
                      y        = 0.503,
                      label    = "Random Chance (0.50)",
                      colour   = "#f8f8f2",
                      hjust    = 0,
                      size     = 5,
                      fontface = "italic",
                      family   = "Arial") +

    # Lines
    ggplot2::geom_line(
      ggplot2::aes(y = ci_ec,
                   colour = "Event-Censored (C_ec)"),
      linewidth = 1.2, na.rm = TRUE) +
    ggplot2::geom_point(
      ggplot2::aes(y = ci_ec,
                   colour = "Event-Censored (C_ec)"),
      size = 3, na.rm = TRUE) +

    ggplot2::geom_line(
      ggplot2::aes(y = ci_ee,
                   colour = "Event-Event (C_ee)"),
      linewidth = 1.2, na.rm = TRUE) +
    ggplot2::geom_point(
      ggplot2::aes(y = ci_ee,
                   colour = "Event-Event (C_ee)"),
      size = 3, na.rm = TRUE) +

    ggplot2::geom_line(
      ggplot2::aes(y = global_c,
                   colour = "Global C-Index"),
      linewidth = 1.2, na.rm = TRUE) +
    ggplot2::geom_point(
      ggplot2::aes(y = global_c,
                   colour = "Global C-Index"),
      size = 3, na.rm = TRUE) +

    # End of line labels
    ggplot2::annotate("text",
                      x        = last_point$censoring + 0.02,
                      y        = last_point$ci_ec,
                      label    = "Event-Censored (C_ec)",
                      colour   = col_ec,
                      hjust    = 0, size = 5,
                      fontface = "bold",
                      family   = "Arial") +
    ggplot2::annotate("text",
                      x        = last_point$censoring + 0.02,
                      y        = last_point$global_c,
                      label    = "Global C-Index",
                      colour   = "#f8f8f2",
                      hjust    = 0, size = 5,
                      fontface = "bold",
                      family   = "Arial") +
    ggplot2::annotate("text",
                      x        = last_point$censoring + 0.02,
                      y        = last_point$ci_ee + 0.008,
                      label    = "Event-Event (C_ee)",
                      colour   = col_ee,
                      hjust    = 0, size = 5,
                      fontface = "bold",
                      family   = "Arial") +

    # Scales
    ggplot2::scale_colour_manual(
      name   = NULL,
      values = c(
        "Event-Censored (C_ec)" = col_ec,
        "Event-Event (C_ee)"    = col_ee,
        "Global C-Index"        = "#f8f8f2"
      )
    ) +
    ggplot2::scale_x_continuous(
      labels  = scales::percent_format(accuracy = 1),
      name    = "Artificial Censoring Rate",
      expand  = ggplot2::expansion(mult = c(0.02, 0.18))
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(0.40, 0.80, by = 0.02),
      name   = "Concordance Index"
    ) +

    # Theme
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::theme(
      plot.background  = ggplot2::element_rect(fill   = background,
                                               colour = NA),
      panel.background = ggplot2::element_rect(fill   = background,
                                               colour = NA),
      text             = ggplot2::element_text(colour = "#f8f8f2",
                                               family = "Arial"),
      axis.text        = ggplot2::element_text(colour = "#f8f8f2",
                                               size   = 22),
      axis.title       = ggplot2::element_text(colour = "#f8f8f2",
                                               face   = "bold",
                                               size   = 24),
      axis.line        = ggplot2::element_line(colour = "#f8f8f2"),
      axis.ticks       = ggplot2::element_line(colour = "#f8f8f2"),
      panel.grid.major = ggplot2::element_line(colour    = "#3a3c4e",
                                               linewidth = 0.4),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "none",
      axis.title.x     = ggplot2::element_text(face   = "bold",
                                               margin = ggplot2::margin(t = 15)),
      axis.title.y     = ggplot2::element_text(face   = "bold",
                                               margin = ggplot2::margin(r = 15)),
      plot.margin      = ggplot2::margin(20, 120, 20, 40)
    )
}