# Calculate Log Bayes Factors for Variants

Vectorized C++ helper for computing variant-level log Bayes Factors.
Uses adaptive Simpson's rule in log-space to prevent numerical
over/underflow.

## Usage

``` r
calc_log_bf(var_case, var_contr, gamma, sigma, N1, N0)
```

## Arguments

- var_case:

  Integer vector of variant counts in cases.

- var_contr:

  Integer vector of variant counts in controls.

- gamma:

  Numeric vector of shape hyperparameters.

- sigma:

  Numeric vector of scale hyperparameters.

- N1:

  Total sample size in cases.

- N0:

  Total sample size in controls.

## Value

Numeric vector of log Bayes Factors.
