# Plot theme for cindexdecomp

Light by default, because journal figures print on white paper. Pass
`dark = TRUE` for the high-contrast presentation variant.

No font family is set anywhere: the system default is used, so plots
render without warnings on machines that lack any particular typeface.

## Usage

``` r
theme_cindex(dark = FALSE, base_size = 14)
```

## Arguments

- dark:

  Logical. `FALSE` (default) for the light theme.

- base_size:

  Base font size in points.

## Value

A ggplot2 theme object.

## Examples

``` r
library(ggplot2)
ggplot(mtcars, aes(wt, mpg)) +
  geom_point() +
  theme_cindex()
```
