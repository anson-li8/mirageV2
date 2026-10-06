#' Validate and Clean MIRAGE Input Data
#'
#' @param data A data frame containing variant counts.
#' @param type Either "gene" or "vs".
#' @return A cleaned data frame with S3 class `mirage_data`.
#' @keywords internal
validate_mirage_data <- function(data, type = c("gene", "vs")) {
  type <- match.arg(type)

  if (!is.data.frame(data)) {
    cli::cli_abort("{.arg data} must be a data frame.")
  }

  # Accept 'group.index' as alias for 'category' (original mirage convention)
  if ("group.index" %in% names(data) && !("category" %in% names(data))) {
    names(data)[names(data) == "group.index"] <- "category"
  }

  # Handle legacy 4-column input for gene-level (Gene, No.case, No.contr, category)
  if (type == "gene" && ncol(data) == 4 && !("ID" %in% names(data))) {
    cli::cli_warn("Input has 4 columns. Assuming columns are {.val Gene}, {.val No.case}, {.val No.contr}, {.val category} and generating {.val ID}.")
    data <- cbind(ID = as.character(seq_len(nrow(data))), data)
  }

  # Handle legacy 3-column input for vs-level
  if (type == "vs" && ncol(data) == 3 && !("ID" %in% names(data))) {
    cli::cli_warn("Input has 3 columns. Assuming columns are {.val No.case}, {.val No.contr}, {.val category} and generating {.val ID}.")
    data <- cbind(ID = as.character(seq_len(nrow(data))), data)
  }

  req_cols <- if (type == "gene") {
    c("ID", "Gene", "No.case", "No.contr", "category")
  } else {
    c("ID", "No.case", "No.contr", "category")
  }

  if (!all(req_cols %in% names(data))) {
    cli::cli_abort("{.arg data} must contain columns: {.val {req_cols}}.")
  }

  # Subset to required columns in exact order
  data <- data[, req_cols, drop = FALSE]


  if (any(is.na(data))) {
    cli::cli_abort(paste0(
      "Input {.arg data} contains missing values ({.val NA}). ",
      "Please remove or impute them before running MIRAGE."
    ))
  }

  # Type checks and coercion
  data$ID <- as.character(data$ID)
  if (type == "gene") data$Gene <- as.character(data$Gene)

  if (!is.numeric(data$No.case) || any(data$No.case < 0)) {
    cli::cli_abort("{.field No.case} must be non-negative numeric values.")
  }
  if (!is.numeric(data$No.contr) || any(data$No.contr < 0)) {
    cli::cli_abort("{.field No.contr} must be non-negative numeric values.")
  }

  data$No.case <- as.integer(round(data$No.case))
  data$No.contr <- as.integer(round(data$No.contr))

  # Category validation & re-indexing to 1:K contiguous integers
  if (is.factor(data$category)) {
    unique_cats <- levels(data$category)
    data$category <- as.integer(data$category)
  } else if (is.numeric(data$category) || is.character(data$category)) {
    unique_cats <- unique(data$category)
    cat_map <- setNames(seq_along(unique_cats), as.character(unique_cats))
    data$category <- as.integer(cat_map[as.character(data$category)])
  } else {
    cli::cli_abort("{.field category} must be numeric, character, or a factor.")
  }

  attr(data, "original_categories") <- unique_cats
  class(data) <- c("mirage_data", "data.frame")

  return(data)
}

#' Validate MIRAGE EM Algorithm Parameters
#' @keywords internal
validate_mirage_params <- function(n1, n2, gamma, sigma, eta.init, delta.init,
                                   estimate.delta, estimate.eta, fixed.eta,
                                   max.iter, tol, num_categories, type = "gene") {
  if (!is.numeric(n1) || length(n1) != 1 || n1 <= 0) {
    cli::cli_abort("{.arg n1} must be a single positive number.")
  }
  if (!is.numeric(n2) || length(n2) != 1 || n2 <= 0) {
    cli::cli_abort("{.arg n2} must be a single positive number.")
  }

  if (!is.numeric(gamma) || any(gamma <= 0)) {
    cli::cli_abort("{.arg gamma} must be positive numeric values.")
  }
  if (!is.numeric(sigma) || any(sigma <= 0)) {
    cli::cli_abort("{.arg sigma} must be positive numeric values.")
  }
  if (!(length(gamma) %in% c(1, num_categories))) {
    cli::cli_abort("{.arg gamma} must have length 1 or match the number of categories ({num_categories}).")
  }
  if (!(length(sigma) %in% c(1, num_categories))) {
    cli::cli_abort("{.arg sigma} must have length 1 or match the number of categories ({num_categories}).")
  }

  if (!is.numeric(eta.init) || length(eta.init) != 1 || eta.init < 0 || eta.init > 1) {
    cli::cli_abort("{.arg eta.init} must be a single numeric value between 0 and 1.")
  }

  if (type == "gene") {
    if (!is.numeric(delta.init) || length(delta.init) != 1 || delta.init < 0 || delta.init > 1) {
      cli::cli_abort("{.arg delta.init} must be a single numeric value between 0 and 1.")
    }
    if (!is.logical(estimate.delta) || length(estimate.delta) != 1) {
      cli::cli_abort("{.arg estimate.delta} must be TRUE or FALSE.")
    }
  }

  if (!is.logical(estimate.eta) || length(estimate.eta) != 1) {
    cli::cli_abort("{.arg estimate.eta} must be TRUE or FALSE.")
  }

  if (!estimate.eta) {
    if (is.null(fixed.eta)) {
      cli::cli_abort("{.arg fixed.eta} must be provided when {.arg estimate.eta} is FALSE.")
    }
    if (!is.numeric(fixed.eta) || length(fixed.eta) != num_categories) {
      cli::cli_abort("{.arg fixed.eta} must be a numeric vector of length {num_categories}.")
    }
    if (any(fixed.eta < 0 | fixed.eta > 1)) {
      cli::cli_abort("All values in {.arg fixed.eta} must be between 0 and 1.")
    }
  }

  if (!is.numeric(max.iter) || length(max.iter) != 1 || max.iter < 1) {
    cli::cli_abort("{.arg max.iter} must be a single positive integer.")
  }

  if (!is.numeric(tol) || length(tol) != 1 || tol <= 0) {
    cli::cli_abort("{.arg tol} must be a single positive number.")
  }

  invisible(TRUE)
}
