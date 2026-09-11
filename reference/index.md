# Package index

## Decomposition

The entry point, and the bootstrap intervals that say whether the
masking gap is distinguishable from zero.

- [`decompose_cindex()`](https://Lainsm.github.io/cindexdecomp/reference/decompose_cindex.md)
  [`print(`*`<cindex_decomp>`*`)`](https://Lainsm.github.io/cindexdecomp/reference/decompose_cindex.md)
  [`summary(`*`<cindex_decomp>`*`)`](https://Lainsm.github.io/cindexdecomp/reference/decompose_cindex.md)
  [`as.data.frame(`*`<cindex_decomp>`*`)`](https://Lainsm.github.io/cindexdecomp/reference/decompose_cindex.md)
  : Decompose a concordance index into Event-Event and Event-Censored
  parts
- [`confint(`*`<cindex_decomp>`*`)`](https://Lainsm.github.io/cindexdecomp/reference/confint.cindex_decomp.md)
  : Bootstrap confidence intervals for a concordance decomposition

## Weightings

A concordance index is fully determined by its pair weighting. These
constructors supply it.

- [`weights_harrell()`](https://Lainsm.github.io/cindexdecomp/reference/cindex_weights.md)
  [`weights_uno()`](https://Lainsm.github.io/cindexdecomp/reference/cindex_weights.md)
  [`weights_truncated()`](https://Lainsm.github.io/cindexdecomp/reference/cindex_weights.md)
  [`weights_custom()`](https://Lainsm.github.io/cindexdecomp/reference/cindex_weights.md)
  : Pair weightings for concordance decomposition

## Censoring and comparison

Trace the decomposition as censoring rises, or run several models over
one cohort.

- [`censoring_curve()`](https://Lainsm.github.io/cindexdecomp/reference/censoring_curve.md)
  [`print(`*`<cindex_curve>`*`)`](https://Lainsm.github.io/cindexdecomp/reference/censoring_curve.md)
  [`as.data.frame(`*`<cindex_curve>`*`)`](https://Lainsm.github.io/cindexdecomp/reference/censoring_curve.md)
  : Trace the decomposition as administrative censoring increases
- [`compare_decompositions()`](https://Lainsm.github.io/cindexdecomp/reference/compare_decompositions.md)
  [`print(`*`<cindex_comparison>`*`)`](https://Lainsm.github.io/cindexdecomp/reference/compare_decompositions.md)
  [`as.data.frame(`*`<cindex_comparison>`*`)`](https://Lainsm.github.io/cindexdecomp/reference/compare_decompositions.md)
  : Compare the decomposition across several models

## Plotting

- [`autoplot(`*`<cindex_comparison>`*`)`](https://Lainsm.github.io/cindexdecomp/reference/autoplot.cindexdecomp.md)
  [`autoplot(`*`<cindex_decomp>`*`)`](https://Lainsm.github.io/cindexdecomp/reference/autoplot.cindexdecomp.md)
  [`autoplot(`*`<cindex_curve>`*`)`](https://Lainsm.github.io/cindexdecomp/reference/autoplot.cindexdecomp.md)
  : Plot a concordance decomposition
- [`theme_cindex()`](https://Lainsm.github.io/cindexdecomp/reference/theme_cindex.md)
  : Plot theme for cindexdecomp

## Package

- [`cindexdecomp`](https://Lainsm.github.io/cindexdecomp/reference/cindexdecomp-package.md)
  [`cindexdecomp-package`](https://Lainsm.github.io/cindexdecomp/reference/cindexdecomp-package.md)
  : cindexdecomp: C-Index Decomposition for Survival Models
- [`reexports`](https://Lainsm.github.io/cindexdecomp/reference/reexports.md)
  [`autoplot`](https://Lainsm.github.io/cindexdecomp/reference/reexports.md)
  : Objects exported from other packages
