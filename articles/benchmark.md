# Performance Benchmark: mirageV2 vs mirage

## Overview

`mirageV2` is a performance-oriented rewrite of the original
[`mirage`](https://github.com/xinhe-lab/mirage) package. This vignette
benchmarks both packages on the same dataset and compares results for
numerical equivalence.

### Architecture differences

| Component | `mirage` | `mirageV2` |
|----|----|----|
| Bayes factor integration | R [`integrate()`](https://rdrr.io/r/stats/integrate.html) per variant in a loop | C++ composite Simpson’s (vectorized) |
| EM algorithm | R `vapply`/nested loops, probability space | C++ EM, log-space with `log_sum_exp` |
| Numerical overflow handling | `ifelse(bb==Inf, 3e300, bb)` patch | Log-space arithmetic (no overflow possible) |

## Setup

``` r

library(mirageV2)

toy <- read.csv(
  system.file("extdata", "mirage_toy.csv", package = "mirageV2"),
  stringsAsFactors = FALSE
)
cat("Dataset:", nrow(toy), "variants across",
    length(unique(toy$Gene)), "genes\n")
#> Dataset: 1000 variants across 282 genes
```

Check if the original `mirage` package is available:

``` r

has_original <- requireNamespace("mirage", quietly = TRUE)
cat("Original mirage installed:", has_original, "\n")
#> Original mirage installed: TRUE
```

## Runtime comparison

``` r

# Time mirageV2
t_v2 <- system.time(
  res_v2 <- mirage(toy, n1 = 4315, n2 = 4315, verbose = FALSE)
)

# Time original mirage
t_v1 <- system.time(
  res_v1 <- mirage::mirage(toy, n1 = 4315, n2 = 4315, verbose = FALSE)
)

cat("\n--- Timing Results ---\n")
#> 
#> --- Timing Results ---
cat(sprintf("mirage  (original): %.2fs\n", t_v1["elapsed"]))
#> mirage  (original): 79.88s
cat(sprintf("mirageV2:           %.2fs\n", t_v2["elapsed"]))
#> mirageV2:           0.71s
cat(sprintf("Speedup:            %.1fx\n", t_v1["elapsed"] / t_v2["elapsed"]))
#> Speedup:            113.1x
```

``` r

# Original mirage not installed; benchmark mirageV2 alone
t_v2 <- system.time(
  res_v2 <- mirage(toy, n1 = 4315, n2 = 4315, verbose = FALSE)
)
cat(sprintf("mirageV2: %.2fs\n", t_v2["elapsed"]))
cat("Install the original package for a direct comparison:\n")
cat("  devtools::install_github('xinhe-lab/mirage')\n")
```

## Result equivalence

``` r

# Compare posterior probabilities
v2_pp <- res_v2$BF.PP.gene$post.prob
v1_pp <- res_v1$BF.PP.gene[, 3]  # post.prob column

cat("Max absolute difference in posterior probabilities:",
    max(abs(v2_pp - v1_pp)), "\n")
#> Max absolute difference in posterior probabilities: 4.594879e-05

# Compare delta estimates
cat(sprintf("mirage  delta: %.6f\n", res_v1$delta.est$delta.est))
#> mirage  delta: 0.925347
cat(sprintf("mirageV2 delta: %.6f\n", res_v2$delta.est$delta.est))
#> mirageV2 delta: 0.925379

# Compare eta estimates
cat("\nEta estimates (mirage vs mirageV2):\n")
#> 
#> Eta estimates (mirage vs mirageV2):
print(data.frame(
  category = res_v2$eta.est |> rownames(),
  mirage = res_v1$eta.est$eta.est,
  mirageV2 = res_v2$eta.est$eta.est
))
#>   category        mirage   mirageV2
#> 1        1  0.000000e+00 0.00000000
#> 2        2  1.030961e-01 0.10309231
#> 3        3 4.940656e-324 0.00000000
#> 4        4 1.976263e-323 0.00000000
#> 5        5  7.871498e-02 0.07871355
```

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
#> [1] mirageV2_0.1.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] vctrs_0.7.3       crayon_1.5.3      cli_3.6.6         knitr_1.52       
#>  [5] rlang_1.3.0       xfun_0.61         otel_0.2.0        generics_0.1.4   
#>  [9] textshaping_1.0.5 jsonlite_2.0.0    prettyunits_1.2.0 htmltools_0.5.9  
#> [13] mirage_0.1.0.0    ragg_1.5.2        sass_0.4.10       hms_1.1.4        
#> [17] rmarkdown_2.32    evaluate_1.0.5    jquerylib_0.1.4   fastmap_1.2.0    
#> [21] progress_1.2.3    yaml_2.3.12       lifecycle_1.0.5   compiler_4.6.1   
#> [25] fs_2.1.0          pkgconfig_2.0.3   Rcpp_1.1.2        systemfonts_1.3.2
#> [29] digest_0.6.39     R6_2.6.1          bslib_0.12.0      tools_4.6.1      
#> [33] pkgdown_2.2.1     cachem_1.1.0      desc_1.4.3
```
