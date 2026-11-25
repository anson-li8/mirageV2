# Changelog

## mirageV2 0.1.0

- First release of a rewrite of the MIRAGE rare-variant association
  tests (Han et al., 2025).
- [`mirage()`](https://anson-li8.github.io/mirageV2/reference/mirage.md)
  fits the gene-level mixture model.
  [`mirage_vs()`](https://anson-li8.github.io/mirageV2/reference/mirage_vs.md)
  fits the variant-set model.
- Variant Bayes factors are computed with a C++ composite Simpson’s
  rule. The EM algorithm runs in log-space, and convergence tracks both
  delta and eta.
- Reported Bayes factors stay finite, and gene posteriors are taken from
  the log-space EM fit.
- Results are S3 objects with
  [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html), and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
  methods. Inputs are validated before the fit.
