# Bootstrap confidence intervals for a concordance decomposition

Bootstrap confidence intervals for a concordance decomposition

## Usage

``` r
# S3 method for class 'cindex_decomp'
confint(object, parm = c("C_ee", "C_ec", "C_global", "gap"), level = NULL, ...)
```

## Arguments

- object:

  A `cindex_decomp` object created with `n_boot > 0`.

- parm:

  Character vector selecting quantities. Defaults to all of `"C_ee"`,
  `"C_ec"`, `"C_global"`, `"gap"`.

- level:

  Confidence level. Defaults to the object's `conf_level`.

- ...:

  Ignored.

## Value

A matrix of percentile bounds, one row per quantity.

## Examples

``` r
set.seed(1)
time <- rexp(150, rate = 0.1)
status <- rbinom(150, 1, 0.6)
risk <- rnorm(150)
fit <- decompose_cindex(time, status, risk, n_boot = 50)
confint(fit)
#>                2.5 %    97.5 %
#> C_ee      0.40790408 0.5491088
#> C_ec      0.39844130 0.6196709
#> C_global  0.42648720 0.5416201
#> gap      -0.09831969 0.1368331
```
