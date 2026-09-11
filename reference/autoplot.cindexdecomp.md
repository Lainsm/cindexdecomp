# Plot a concordance decomposition

[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
methods for the three result classes. Because the object carries its own
class, the right plot is selected automatically and a result can never
be paired with the wrong chart.

The dumbbell plots (`cindex_decomp`, `cindex_comparison`) draw error
bars from the object's own bootstrap replicates as a percentile interval
at its `conf_level` (e.g. 95%), not a symmetric +/-1 SD band; the level
actually drawn is named in the subtitle. Bars are omitted when no
bootstrap was run (`n_boot = 0`).

## Usage

``` r
# S3 method for class 'cindex_comparison'
autoplot(object, dark = FALSE, ...)

# S3 method for class 'cindex_decomp'
autoplot(object, dark = FALSE, ...)

# S3 method for class 'cindex_curve'
autoplot(object, dark = FALSE, ...)
```

## Arguments

- object:

  A `cindex_decomp`, `cindex_curve` or `cindex_comparison` object.

- dark:

  Logical. Passed to
  [`theme_cindex()`](https://lainsm.github.io/cindexdecomp/reference/theme_cindex.md).

- ...:

  Ignored.

## Value

A `ggplot` object.
