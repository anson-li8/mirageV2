# Stable Log of (1 - eta) + eta \* exp(log_bf) — Vectorized

Stable Log of (1 - eta) + eta \* exp(log_bf) — Vectorized

## Usage

``` r
log_mixture_bf(log_bf, eta)
```

## Arguments

- log_bf:

  Numeric vector of log Bayes factors.

- eta:

  Numeric scalar or vector in \\\\\\\[0, 1\]\\.

## Value

Numeric vector of log((1-eta) + eta \* exp(log_bf)).
