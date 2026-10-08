#' Summary method for mirage_result
#'
#' @param object A \code{mirage_result} object.
#' @param n_top Number of top genes to display. Default 10.
#' @param ... Ignored.
#' @export
summary.mirage_result <- function(object, n_top = 10L, ...) {
  cli::cli_h1("MIRAGE Gene-Level Summary")
  cli::cli_text("Genes analyzed: {nrow(object$BF.PP.gene)}")
  cli::cli_text("Variant categories: {length(object$original_categories)}")
  cli::cli_text("Delta (prop. risk genes): {round(object$delta.est$delta.est, 4)} (p = {signif(object$delta.est$delta.pvalue, 3)})")
  cli::cli_h2("Eta estimates (prop. risk variants per category)")
  print(object$eta.est)
  cli::cli_h2("Top {n_top} genes by posterior probability")
  top <- head(object$BF.PP.gene[order(object$BF.PP.gene$post.prob, decreasing = TRUE), ], n_top)
  print(top, row.names = FALSE)
  invisible(object)
}

#' Summary method for mirage_vs_result
#'
#' @param object A \code{mirage_vs_result} object.
#' @param n_top Number of top variants to display. Default 10.
#' @param ... Ignored.
#' @export
summary.mirage_vs_result <- function(object, n_top = 10L, ...) {
  cli::cli_h1("MIRAGE Variant-Set Summary")
  cli::cli_text("Variants analyzed: {nrow(object$full.info)}")
  cli::cli_text("Variant categories: {length(object$original_categories)}")
  cli::cli_h2("Eta estimates")
  print(object$eta.est)
  cli::cli_h2("Top {n_top} variants by posterior probability")
  top <- head(object$post.prob[order(object$post.prob$post.prob, decreasing = TRUE), ], n_top)
  print(top, row.names = FALSE)
  invisible(object)
}

#' Tidy method for mirage_result
#'
#' Returns gene-level results as a tidy data frame.
#'
#' @param x A \code{mirage_result} object.
#' @param ... Ignored.
#' @return A data frame with columns: gene, bf, post.prob, p.value.
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.mirage_result <- function(x, ...) {
  df <- x$BF.PP.gene
  names(df) <- c("gene", "bf", "post.prob")
  df
}

#' Tidy method for mirage_vs_result
#'
#' Returns variant-level results as a tidy data frame.
#'
#' @param x A \code{mirage_vs_result} object.
#' @param ... Ignored.
#' @return A data frame with columns: variant, bf, post.prob, category.
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.mirage_vs_result <- function(x, ...) {
  df <- data.frame(
    variant = x$post.prob$variant,
    bf = x$post.prob$BF,
    post.prob = x$post.prob$post.prob,
    category = x$full.info$category,
    stringsAsFactors = FALSE
  )
  df
}

#' Coerce mirage_result to data frame
#'
#' @param x A \code{mirage_result} object.
#' @param ... Ignored.
#' @return A data frame of gene-level results.
#' @export
as.data.frame.mirage_result <- function(x, ...) {
  x$BF.PP.gene
}

#' Coerce mirage_vs_result to data frame
#'
#' @param x A \code{mirage_vs_result} object.
#' @param ... Ignored.
#' @return A data frame of variant-level results.
#' @export
as.data.frame.mirage_vs_result <- function(x, ...) {
  x$full.info
}
