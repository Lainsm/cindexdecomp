#' Colour palette for cindexdecomp plots
#'
#' Event-Event and Event-Censored keep the IBM colourblind-safe magenta and
#' blue in both variants; only the neutrals change.
#'
#' @param dark Logical. `FALSE` (default) returns the light palette.
#' @return A named list of hex colours.
#' @keywords internal
#' @noRd
cindex_palette <- function(dark = FALSE) {
  list(
    ee     = "#DC267F",
    ec     = "#648FFF",
    global = if (dark) "#F8F8F2" else "#2B2B2B",
    fg     = if (dark) "#F8F8F2" else "#1A1A1A",
    bg     = if (dark) "#282A36" else "#FFFFFF",
    grid   = if (dark) "#3A3C4E" else "#E6E6E6",
    track  = if (dark) "#3A3C4E" else "#EDEDED"
  )
}

#' Plot theme for cindexdecomp
#'
#' @description
#' Light by default, because journal figures print on white paper. Pass
#' `dark = TRUE` for the high-contrast presentation variant.
#'
#' No font family is set anywhere: the system default is used, so plots
#' render without warnings on machines that lack any particular typeface.
#'
#' @param dark Logical. `FALSE` (default) for the light theme.
#' @param base_size Base font size in points.
#' @return A ggplot2 theme object.
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) + geom_point() + theme_cindex()
#' @export
theme_cindex <- function(dark = FALSE, base_size = 14) {
  p <- cindex_palette(dark)
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.background   = ggplot2::element_rect(fill = p$bg, colour = NA),
      panel.background  = ggplot2::element_rect(fill = p$bg, colour = NA),
      legend.background = ggplot2::element_rect(fill = p$bg, colour = NA),
      legend.key        = ggplot2::element_rect(fill = p$bg, colour = NA),
      text              = ggplot2::element_text(colour = p$fg),
      axis.text         = ggplot2::element_text(colour = p$fg),
      axis.line         = ggplot2::element_line(colour = p$fg),
      axis.ticks        = ggplot2::element_line(colour = p$fg),
      panel.grid.major  = ggplot2::element_line(colour = p$grid,
                                                linewidth = 0.3),
      panel.grid.minor  = ggplot2::element_blank(),
      plot.title        = ggplot2::element_text(face = "bold"),
      plot.subtitle     = ggplot2::element_text(colour = p$fg),
      plot.margin       = ggplot2::margin(12, 16, 12, 12)
    )
}
