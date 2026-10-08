# Getting Started with mirageV2

## Installation

``` r

devtools::install_github("anson-li8/mirageV2")
```

## Input format

`mirageV2` requires a data frame with variant-level summary statistics.

### Gene-level (`mirage`)

| Column     | Description                                 |
|------------|---------------------------------------------|
| `ID`       | Variant identifier                          |
| `Gene`     | Gene name                                   |
| `No.case`  | Allele count in cases                       |
| `No.contr` | Allele count in controls                    |
| `category` | Variant annotation group (or `group.index`) |

### Variant-set level (`mirage_vs`)

Same as above without the `Gene` column.

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

Run the analysis:

``` r

res <- mirage(toy, n1 = 4315, n2 = 4315)
#> ℹ Computing variant-level Bayes factors
#> ✔ Computing variant-level Bayes factors [329ms]
#> 
#> ℹ Running EM algorithm (282 genes)
#> ✔ EM converged in 2947 iterations
#> ℹ Running EM algorithm (282 genes)✔ Running EM algorithm (282 genes) [348ms]
#> 
#> ℹ Computing LRT statistics and p-values
#> ✔ Computing LRT statistics and p-values [149ms]
```

The result is a `mirage_result` object:

``` r

print(res)
#> 
#> ── MIRAGE Gene-Level Analysis Result ───────────────────────────────────────────
#> Number of genes analyzed: 282
#> Number of variant categories: 5
#> Estimated proportion of risk genes (delta): 0.9254
#> 
#> ── Top 5 Genes by Posterior Probability ──
#> 
#>    Gene       BF post.prob
#>   MACF1 2.135502 0.9636131
#>    UBR4 1.289898 0.9411628
#>    MIB1 1.251474 0.9394656
#>    SPEN 1.168992 0.9354701
#>  MYCBP2 1.146437 0.9342839
```

Key outputs:

- `res$delta.est`: Estimated proportion of risk genes with p-value.
- `res$eta.est`: Estimated proportion of risk variants per category with
  p-values.
- `res$BF.PP.gene`: Bayes factors and posterior probabilities per gene.
- `res$BF.all`: Variant-level Bayes factors grouped by gene.

``` r

summary(res)
#> ── MIRAGE Gene-Level Summary ───────────────────────────────────────────────────
#> Genes analyzed: 282
#> Variant categories: 5
#> Delta (prop. risk genes): 0.9254 (p = 0.935)
#> 
#> ── Eta estimates (prop. risk variants per category) ──
#> 
#>      eta.est eta.pvalue
#> 1 0.00000000  0.9999999
#> 2 0.10309231  0.4206781
#> 3 0.00000000  0.9999999
#> 4 0.00000000  0.9999999
#> 5 0.07871355  0.2782095
#> ── Top 10 genes by posterior probability ──
#>     Gene       BF post.prob
#>    MACF1 2.135502 0.9636131
#>     UBR4 1.289898 0.9411628
#>     MIB1 1.251474 0.9394656
#>     SPEN 1.168992 0.9354701
#>   MYCBP2 1.146437 0.9342839
#>     BAI2 1.127566 0.9332575
#>  RPS6KA1 1.123288 0.9330204
#>    DIP2C 1.117268 0.9326838
#>     DHX8 1.117268 0.9326838
#>    LTBP3 1.095385 0.9314312
```

Tidy output:

``` r

head(tidy(res))
#>     gene        bf post.prob
#> 1    LPA 0.9164302 0.9191244
#> 2   MTOR 1.0405788 0.9280794
#> 3 PTCHD2 0.8327233 0.9117123
#> 4   UBR4 1.2898984 0.9411628
#> 5  CSMD2 0.8010658 0.9085424
#> 6   NCDN 1.0006240 0.9254218
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

``` r

res_vs <- mirage_vs(vs_toy, n1 = 4315, n2 = 4315)
#> Warning: Input has 3 columns. Assuming columns are "No.case", "No.contr", "category" and
#> generating "ID".
#> ℹ Computing variant-level Bayes factors
#> ✔ Computing variant-level Bayes factors [74ms]
#> 
#> ℹ Running EM algorithm (213 variants)
#> ✔ EM converged in 47 iterations
#> ℹ Running EM algorithm (213 variants)✔ Running EM algorithm (213 variants) [19ms]
#> 
#> ℹ Computing LRT statistics and p-values
#> ✔ Computing LRT statistics and p-values [15ms]
print(res_vs)
#> 
#> ── MIRAGE Variant-Set Analysis Result ──────────────────────────────────────────
#> Number of variants analyzed: 213
#> Number of variant categories: 1
#> 
#> ── Category-specific risk proportions (eta) ──
#>    eta.est eta.pvalue
#>  0.1063957 0.05768073
```

## Fixed eta

When external annotation data provides reliable risk proportions, fix
`eta` instead of estimating it:

``` r

res_fixed <- mirage(toy, n1 = 4315, n2 = 4315,
                    estimate.eta = FALSE,
                    fixed_eta = c(0.60, 0.31, 0.16, 0.40, 0.20))
#> ℹ Computing variant-level Bayes factors
#> ✔ Computing variant-level Bayes factors [314ms]
#> 
#> ℹ Running EM algorithm (282 genes)
#> ✔ EM converged in 2 iterations
#> ℹ Running EM algorithm (282 genes)✔ Running EM algorithm (282 genes) [18ms]
#> 
#> ℹ Computing LRT statistics and p-values
#> ✔ Computing LRT statistics and p-values [135ms]
head(res_fixed$BF.PP.gene)
#>     Gene          BF    post.prob
#> 1    LPA 0.009867872 0.0009560602
#> 2   MTOR 0.984854125 0.0871831704
#> 3 PTCHD2 0.585846233 0.0537603275
#> 4   UBR4 1.335285048 0.1146481371
#> 5  CSMD2 0.077842145 0.0074924820
#> 6   NCDN 1.023760656 0.0903162663
```

## Key parameters

| Parameter    | Default | Description                       |
|--------------|---------|-----------------------------------|
| `gamma`      | 3       | Gamma prior shape for effect size |
| `sigma`      | 2       | Gamma prior scale for effect size |
| `eta.init`   | 0.1     | Initial risk variant proportion   |
| `delta.init` | 0.1     | Initial risk gene proportion      |
| `max.iter`   | 10000   | Maximum EM iterations             |
| `tol`        | 1e-5    | Convergence tolerance             |

## Session info

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.5 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
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
#>  [5] xfun_0.61         cachem_1.1.0      knitr_1.52        htmltools_0.5.9  
#>  [9] rmarkdown_2.32    lifecycle_1.0.5   cli_3.6.6         sass_0.4.10      
#> [13] pkgdown_2.2.1     textshaping_1.0.5 jquerylib_0.1.4   systemfonts_1.3.2
#> [17] compiler_4.6.1    tools_4.6.1       ragg_1.5.2        bslib_0.12.0     
#> [21] evaluate_1.0.5    Rcpp_1.1.2        yaml_2.3.12       otel_0.2.0       
#> [25] jsonlite_2.0.0    rlang_1.3.0       fs_2.1.0
```
