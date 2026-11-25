# Getting Started with mirageV2

## Installation

The package is installed from GitHub:

``` r

devtools::install_github("anson-li8/mirageV2")
```

## Input format

Both functions take a data frame of variant-level summary statistics.
The required columns are:

### Gene-level (`mirage`)

| Column | Description |
|----|----|
| `ID` | Variant identifier |
| `Gene` | Gene name |
| `No.case` | Allele count in cases |
| `No.contr` | Allele count in controls |
| `category` | Variant annotation group; `group.index` is accepted as an alias |

### Variant-set level (`mirage_vs`)

Same columns, without the `Gene` column.

## Gene-level analysis

``` r

library(mirageV2)
library(generics)
#> 
#> Attaching package: 'generics'
#> The following objects are masked from 'package:base':
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union

toy <- read.csv(
  system.file("extdata", "mirage_toy.csv", package = "mirageV2"),
  stringsAsFactors = FALSE
)
head(toy)
#>                             ID   Gene No.case No.contr group.index
#> 1    6:161006077-161006077_C_T    LPA     298      272           1
#> 2      1:11205024-11205024_C_T   MTOR       1        0           2
#> 3      1:11272449-11272449_G_A   MTOR       0        1           2
#> 4 1:11290989-11290994_CTGACT_C   MTOR       1        0           2
#> 5      1:11594465-11594465_C_T PTCHD2       0        1           2
#> 6   1:19441304-19441307_CCTT_C   UBR4       2        0           2
```

Fit the model with 4315 cases and 4315 controls:

``` r

res <- mirage(toy, n1 = 4315, n2 = 4315)
#> ℹ Computing variant-level Bayes factors
#> ℹ Running EM algorithm over 282 genes
#> ✔ EM converged in 5303 iterations
#> ℹ Computing LRT statistics and p-values
```

[`mirage()`](https://anson-li8.github.io/mirageV2/reference/mirage.md)
returns a `mirage_result` object:

``` r

print(res)
#> 
#> ── MIRAGE Gene-Level Analysis Result ───────────────────────────────────────────
#> Number of genes analyzed: 282
#> Number of variant categories: 5
#> Estimated proportion of risk genes (delta): 0.9894
#>  category    eta.est eta.pvalue
#>         1 0.00000000  0.9999998
#>         2 0.09617393  0.4198750
#>         3 0.00000000  0.9999998
#>         4 0.00000000  0.9999998
#>         5 0.07594910  0.2715666
#> 
#> ── Top 5 Genes by Posterior Probability ──
#> 
#>    Gene       BF post.prob
#>   MACF1 2.097148 0.9949148
#>    UBR4 1.278714 0.9916871
#>    MIB1 1.234471 0.9913918
#>    SPEN 1.163957 0.9908750
#>  MYCBP2 1.136189 0.9906541
```

The main elements are:

- `res$delta.est`: estimated proportion of risk genes, with a p-value.
- `res$eta.est`: estimated proportion of risk variants per category,
  with p-values.
- `res$BF.PP.gene`: Bayes factors and posterior probabilities per gene.
- `res$BF.all`: variant-level Bayes factors grouped by gene.

``` r

summary(res)
#> ── MIRAGE Gene-Level Summary ───────────────────────────────────────────────────
#> Genes analyzed: 282
#> Variant categories: 5
#> Delta (prop. risk genes): 0.9894 (p = 0.932)
#> 
#> ── Eta estimates (prop. risk variants per category) ──
#> 
#>      eta.est eta.pvalue
#> 1 0.00000000  0.9999998
#> 2 0.09617393  0.4198750
#> 3 0.00000000  0.9999998
#> 4 0.00000000  0.9999998
#> 5 0.07594910  0.2715666
#> ── Top 10 genes by posterior probability ──
#>     Gene       BF post.prob
#>    MACF1 2.097148 0.9949148
#>     UBR4 1.278714 0.9916871
#>     MIB1 1.234471 0.9913918
#>     SPEN 1.163957 0.9908750
#>   MYCBP2 1.136189 0.9906541
#>     BAI2 1.122968 0.9905451
#>  RPS6KA1 1.119002 0.9905119
#>    DIP2C 1.109399 0.9904306
#>     DHX8 1.109399 0.9904306
#>     CHD5 1.092099 0.9902804
```

[`tidy()`](https://generics.r-lib.org/reference/tidy.html) returns the
gene-level results as a data frame:

``` r

head(tidy(res))
#>     gene        bf post.prob
#> 1    LPA 0.9230083 0.9885203
#> 2   MTOR 1.0389620 0.9897884
#> 3 PTCHD2 0.8397172 0.9873960
#> 4   UBR4 1.2787141 0.9916871
#> 5  CSMD2 0.8074922 0.9868996
#> 6   NCDN 0.9978913 0.9893726
```

## Variant-set analysis

``` r

vs_toy <- read.csv(
  system.file("extdata", "mirage_vs_toy.csv", package = "mirageV2"),
  stringsAsFactors = FALSE
)
head(vs_toy)
#>   No.case No.contr group.index
#> 1       1        0           1
#> 2       3        1           1
#> 3       2        4           1
#> 4       1        3           1
#> 5       1        0           1
#> 6       2        0           1
```

Fit the variant-set model on the variant-set toy data:

``` r

res_vs <- mirage_vs(vs_toy, n1 = 4315, n2 = 4315)
#> Warning: Input has 3 columns. Assuming columns are "No.case", "No.contr", "category" and
#> generating "ID".
#> ℹ Computing variant-level Bayes factors
#> ℹ Running EM algorithm over 213 variants
#> ✔ EM converged in 47 iterations
#> ℹ Computing LRT statistics and p-values
print(res_vs)
#> 
#> ── MIRAGE Variant-Set Analysis Result ──────────────────────────────────────────
#> Number of variants analyzed: 213
#> Number of variant categories: 1
#> 
#> ── Category-specific risk proportions (eta) ──
#> 
#>  category   eta.est eta.pvalue
#>         1 0.1063957 0.05768073
```

## Fixed eta

If you have external annotation data with risk proportions you trust,
set `estimate.eta = FALSE` and pass them with `fixed_eta`:

``` r

res_fixed <- mirage(toy, n1 = 4315, n2 = 4315,
                    estimate.eta = FALSE,
                    fixed_eta = c(0.60, 0.31, 0.16, 0.40, 0.20))
#> ℹ Computing variant-level Bayes factors
#> ℹ Running EM algorithm over 282 genes
#> ✔ EM converged in 66 iterations
#> ℹ Computing LRT statistics and p-values
head(res_fixed$BF.PP.gene)
#>     Gene          BF    post.prob
#> 1    LPA 0.009867872 9.003802e-07
#> 2   MTOR 0.984854125 8.985364e-05
#> 3 PTCHD2 0.585846233 5.345191e-05
#> 4   UBR4 1.335285048 1.218215e-04
#> 5  CSMD2 0.077842145 7.102553e-06
#> 6   NCDN 1.023760656 9.340297e-05
```

## Key parameters

| Parameter    | Default | Description                               |
|--------------|---------|-------------------------------------------|
| `gamma`      | 3       | Prior mean of the relative-risk parameter |
| `sigma`      | 2       | Prior rate. Gamma shape is gamma \* sigma |
| `eta.init`   | 0.1     | Initial risk variant proportion           |
| `delta.init` | 0.1     | Initial risk gene proportion              |
| `max.iter`   | 10000   | Maximum EM iterations                     |
| `tol`        | 1e-5    | Convergence tolerance                     |

## Session info

``` r

sessionInfo()
#> R version 4.6.2 (2026-10-29)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 26.04.1 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.32.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] generics_0.1.4 mirageV2_0.1.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
#>  [5] xfun_0.62         cachem_1.1.0      knitr_1.52        htmltools_0.5.9  
#>  [9] rmarkdown_2.32    lifecycle_1.0.5   cli_3.6.6         sass_0.5.0       
#> [13] pkgdown_2.2.1     textshaping_1.1.0 jquerylib_0.1.4   systemfonts_1.3.4
#> [17] compiler_4.6.2    tools_4.6.2       ragg_1.5.2        bslib_0.12.0     
#> [21] evaluate_1.0.5    Rcpp_1.1.2        yaml_2.3.12       otel_0.3.0       
#> [25] jsonlite_2.0.0    rlang_1.3.0       fs_2.1.0
```
