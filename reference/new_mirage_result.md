# Construct a mirage_result object

Construct a mirage_result object

## Usage

``` r
new_mirage_result(
  delta_est,
  eta_est,
  bf_pp_gene,
  bf_all,
  original_categories,
  call
)
```

## Arguments

- delta_est:

  Data frame of delta estimates.

- eta_est:

  Data frame of eta estimates.

- bf_pp_gene:

  Data frame of gene-level Bayes Factors and Posterior Probabilities.

- bf_all:

  List of variant-level Bayes Factors per gene.

- original_categories:

  The original category labels before internal re-indexing.

- call:

  The matched call.

## Value

An object of class `mirage_result`.
