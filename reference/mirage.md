# mirage: MIxture model based Rare variant Analysis on GEnes

Gene-level rare variant association test using a Bayesian mixture model.

## Usage

``` r
mirage(
  data,
  n1,
  n2,
  gamma = 3,
  sigma = 2,
  eta.init = 0.1,
  delta.init = 0.1,
  estimate.delta = TRUE,
  estimate.eta = TRUE,
  fixed_eta = NULL,
  max.iter = 10000L,
  tol = 1e-05,
  verbose = TRUE
)
```

## Arguments

- data:

  Data frame with columns: ID, Gene, No.case, No.contr, category. A
  4-column input (without ID) is accepted with a warning.

- n1:

  Sample size in cases.

- n2:

  Sample size in controls.

- gamma:

  Prior mean of the relative-risk parameter. Scalar or length-K vector.

- sigma:

  Prior rate of the relative-risk parameter. Scalar or length-K vector.
  The Gamma prior has shape gamma \* sigma and scale 1 / sigma.

- eta.init:

  Initial value for proportion of risk variants.

- delta.init:

  Initial value for proportion of risk genes.

- estimate.delta:

  Logical. Whether to estimate delta.

- estimate.eta:

  Logical. Whether to estimate eta.

- fixed_eta:

  Numeric vector of fixed eta values (required when estimate.eta =
  FALSE). Must have length K.

- max.iter:

  Maximum EM iterations.

- tol:

  Convergence tolerance.

- verbose:

  Logical. Print progress messages.

## Value

An object of class `mirage_result`.

## Examples

``` r
toy <- data.frame(
  ID = c("v1", "v2", "v3", "v4"),
  Gene = c("A", "A", "B", "B"),
  No.case = c(3L, 1L, 0L, 2L),
  No.contr = c(1L, 0L, 2L, 0L),
  category = c(1L, 2L, 1L, 2L)
)
mirage(toy, n1 = 500, n2 = 500, verbose = FALSE)
#> 
#> ── MIRAGE Gene-Level Analysis Result ───────────────────────────────────────────
#> Number of genes analyzed: 2
#> Number of variant categories: 2
#> Estimated proportion of risk genes (delta): 1
#>  category      eta.est eta.pvalue
#>         1 0.0001946449  1.0000000
#>         2 1.0000000000  0.1323021
#> 
#> ── Top 5 Genes by Posterior Probability ──
#> 
#>  Gene       BF post.prob
#>     B 2.137230         1
#>     A 1.452247         1
```
