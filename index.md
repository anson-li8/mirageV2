# mirageV2

`mirageV2` is a refactor and performance improvement of the
[`mirage`](https://github.com/xinhe-lab/mirage) package for rare variant
(RV) association testing using a Bayesian mixture model. It models
variants in a gene as a mixture of risk and non-risk variants, with
prior probabilities informed by functional annotations such as
conservation scores and protein impact predictions.

## Key improvements over `mirage`

| Area | `mirage` | `mirageV2` |
|----|----|----|
| **Numerical stability** | Probability-space EM with `ifelse(bb==Inf, 3e300, bb)` overflow patches | Log-space EM with `log_sum_exp` throughout |
| **BF computation** | R [`integrate()`](https://rdrr.io/r/stats/integrate.html) per variant in a loop | Vectorized C++ composite Simpson’s rule via Rcpp |
| **EM engine** | Duplicated `vapply`/nested loops in [`mirage()`](reference/mirage.md) and [`mirage_vs()`](reference/mirage_vs.md) | Single unified engine, vectorized with `rowsum`/`table` |
| **Input validation** | Implicit column guessing | Strict S3 validation with `cli` error messages |
| **Output** | Loose unnamed lists | S3 classes (`mirage_result`, `mirage_vs_result`) with `print`, `summary`, `tidy` methods |
| **Testing** | None | 56 `testthat` tests covering edge cases |
| **Performance** | ~160s on `mirage_toy` (282 genes) | **~10s** (15x speedup) |

## Installation

``` r

# Install devtools if needed
install.packages("devtools")

# Install mirageV2
devtools::install_github("anson-li8/mirageV2")
```

## Quick start

``` r

library(mirageV2)

# Load bundled toy data
toy <- read.csv(system.file("extdata", "mirage_toy.csv", package = "mirageV2"),
                stringsAsFactors = FALSE)

# Gene-level analysis
res <- mirage(toy, n1 = 4315, n2 = 4315)
print(res)
summary(res)

# Tidy output (broom-compatible)
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
`mirage`. For a gene $`i`$ with variants $`j = 1, \ldots, m_i`$:

- Each variant is modeled as risk ($`Z_{ij} = 1`$) or non-risk
  ($`Z_{ij} = 0`$) with $`P(Z_{ij} = 1 \mid U_i = 1) = \eta_k`$, where
  $`k`$ is the variant’s annotation category.
- Each gene is a risk gene ($`U_i = 1`$) with probability $`\delta`$.
- Variant-level Bayes factors are computed by integrating over a
  Gamma-distributed relative risk parameter.
- Parameters $`(\delta, \eta_1, \ldots, \eta_K)`$ are estimated via EM
  in log-space.
- Gene-level posterior probabilities and likelihood ratio test p-values
  are reported.

See the [original `mirage`
documentation](https://xinhe-lab.github.io/mirage/) for full method
details.

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
#>  [9] rmarkdown_2.31    knitr_1.51        xfun_0.60         digest_0.6.39    
#> [13] rlang_1.3.0       evaluate_1.0.5
```
