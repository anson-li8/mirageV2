# mirageV2 0.1.0

* First release of a rewrite of the MIRAGE rare-variant association tests (Han et al., 2025).
* `mirage()` fits the gene-level mixture model. `mirage_vs()` fits the variant-set model.
* Variant Bayes factors are computed with a C++ composite Simpson's rule. The EM algorithm runs in log-space, and convergence tracks both delta and eta.
* Reported Bayes factors stay finite, and gene posteriors are taken from the log-space EM fit.
* Results are S3 objects with `print()`, `summary()`, `tidy()`, and `as.data.frame()` methods. Inputs are validated before the fit.
