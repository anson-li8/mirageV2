# MIRAGE EM Algorithm

MIRAGE EM Algorithm

## Usage

``` r
em_mirage(
  type = c("gene", "vs"),
  log_var_bf,
  category,
  gene_index = NULL,
  n_genes = NULL,
  n_categories,
  delta_init = 0.1,
  eta_init = rep(0.1, n_categories),
  estimate_delta = TRUE,
  estimate_eta = TRUE,
  fixed_eta = NULL,
  max_iter = 10000L,
  tol = 1e-05,
  verbose = TRUE
)
```

## Arguments

- type:

  Either "gene" or "vs".

- log_var_bf:

  Numeric vector of log Bayes factors for each variant.

- category:

  Integer vector of category indices (1:K) for each variant.

- gene_index:

  Integer vector of gene indices (1:G) for each variant.

- n_genes:

  Number of genes. Required when type = "gene".

- n_categories:

  Number of variant categories (K).

- delta_init:

  Initial value for delta.

- eta_init:

  Numeric vector of initial eta values (length K).

- estimate_delta:

  Logical.

- estimate_eta:

  Logical.

- fixed_eta:

  Numeric vector of fixed eta values.

- max_iter:

  Maximum number of EM iterations.

- tol:

  Convergence tolerance.

- verbose:

  Logical.

## Value

A list containing EM results.
