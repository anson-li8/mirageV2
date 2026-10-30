
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mirageV2

<!-- badges: start -->

[![R-CMD-check](https://github.com/anson-li8/mirageV2/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/anson-li8/mirageV2/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/anson-li8/mirageV2/graph/badge.svg)](https://app.codecov.io/gh/anson-li8/mirageV2)
<!-- badges: end -->

`mirageV2` is a high-performance rewrite of the
[`mirage`](https://github.com/xinhe-lab/mirage) package for rare variant
(RV) association testing. It fits a Bayesian mixture model in which each
variant in a gene is treated as either a risk or non-risk variant, with
prior probabilities that depend on functional annotations such as
conservation scores and protein-impact predictions.

## Main improvements over `mirage`

| Area | `mirage` | `mirageV2` |
|----|----|----|
| **Numerical stability** | Probability-space EM. Overflow is caught by swapping Inf for 3e300 with ifelse(bb==Inf, 3e300, bb), and underflow leaves denormalized values like 4.9e-324. | Log-space EM using log_sum_exp. Overflow and underflow do not occur. |
| **BF computation** | R `integrate()` per variant in a loop | C++ composite Simpson’s rule (vectorized) |
| **EM engine** | Duplicated R `vapply`/nested loops | One C++ EM implementation |
| **Input validation** | Columns guessed from names | Validated up front; errors use cli |
| **Output** | Unnamed lists | S3 classes with `print`, `summary`, `tidy` methods |
| **Testing** | None | 56 `testthat` tests covering edge cases |
| **Performance** | ~86.69s on `mirage_toy` (282 genes) | **~2.27s** (**38.2x speedup**) |
| **Numerical equivalence** | n/a | Max posterior-probability difference 4.6e-5. delta is 0.925347 vs 0.925379. |

Runtime numbers come from a single run of both packages on the
mirage_toy dataset (1000 variants across 282 genes) on Windows 11 with R
4.6.0. Timings depend on the machine, and the original mirage package
has to be installed for a side-by-side comparison. See the [benchmark
vignette](https://anson-li8.github.io/mirageV2/articles/benchmark.html)
for detailed timing and equivalence results.

## Installation

``` r
# install.packages("devtools")
devtools::install_github("anson-li8/mirageV2")
```

## Quick start

``` r
library(mirageV2)

# Load toy data
toy <- read.csv(system.file("extdata", "mirage_toy.csv", package = "mirageV2"),
                stringsAsFactors = FALSE)

# Gene-level analysis
res <- mirage(toy, n1 = 4315, n2 = 4315)
print(res)
summary(res)

# Tidy output via a broom-style tidy() method
tidy(res)

# Variant-set analysis
vs_toy <- read.csv(system.file("extdata", "mirage_vs_toy.csv", package = "mirageV2"),
                   stringsAsFactors = FALSE)
res_vs <- mirage_vs(vs_toy, n1 = 4315, n2 = 4315)
print(res_vs)
```

## Input format

### Gene-level (`mirage`)

A data frame with 5 columns (or 4 without `ID`):

| Column | Description |
|----|----|
| `ID` | Variant identifier |
| `Gene` | Gene name |
| `No.case` | Allele count in cases |
| `No.contr` | Allele count in controls |
| `category` | Variant annotation group (or `group.index` for backward compatibility) |

### Variant-set level (`mirage_vs`)

Same as above without the `Gene` column.

## Method

`mirageV2` implements the same Bayesian mixture model as the original
`mirage`. For a gene $i$ with variants $j = 1, \ldots, m_i$:

- Each variant is modeled as risk ($Z_{ij} = 1$) or non-risk
  ($Z_{ij} = 0$) with $P(Z_{ij} = 1 \mid U_i = 1) = \eta_k$, where $k$
  is the variant’s annotation category.
- Each gene is a risk gene ($U_i = 1$) with probability $\delta$.
- Variant-level Bayes factors are computed by integrating over a
  Gamma-distributed relative risk parameter.
- Parameters $(\delta, \eta_1, \ldots, \eta_K)$ are estimated via EM
  (Expectation-Maximization) in log-space.
- Gene-level posterior probabilities and likelihood ratio test p-values
  are reported.

See the [original `mirage`
documentation](https://xinhe-lab.github.io/mirage/) for full method
details, and the [method
tutorial](https://anson-li8.github.io/mirageV2/articles/methodtutorial.html)
for a usage walkthrough.

## Session information

``` r
sessionInfo()
#> R version 4.6.0 (2026-04-24 ucrt)
#> Platform: x86_64-w64-mingw32/x64
#> Running under: Windows 11 x64 (build 29639)
#> 
#> Matrix products: default
#>   LAPACK version 3.12.1
#> 
#> locale:
#> [1] LC_COLLATE=Spanish_Latin America.utf8 
#> [2] LC_CTYPE=Spanish_Latin America.utf8   
#> [3] LC_MONETARY=Spanish_Latin America.utf8
#> [4] LC_NUMERIC=C                          
#> [5] LC_TIME=Spanish_Latin America.utf8    
#> 
#> time zone: America/Chicago
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] compiler_4.6.0    fastmap_1.2.0     cli_3.6.6         tools_4.6.0      
#>  [5] htmltools_0.5.9   otel_0.2.0        rstudioapi_0.19.0 yaml_2.3.12      
#>  [9] rmarkdown_2.32    knitr_1.52        xfun_0.62         digest_0.6.39    
#> [13] rlang_1.3.0       evaluate_1.0.5
```
