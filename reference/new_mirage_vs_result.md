# Construct a mirage_vs_result object

Construct a mirage_vs_result object

## Usage

``` r
new_mirage_vs_result(eta_est, full_info, post_prob, original_categories, call)
```

## Arguments

- eta_est:

  Data frame of eta estimates.

- full_info:

  Data frame of variant-level information and Bayes Factors.

- post_prob:

  Data frame of variant-level Posterior Probabilities.

- original_categories:

  The original category labels.

- call:

  The matched call.

## Value

An object of class `mirage_vs_result`.
