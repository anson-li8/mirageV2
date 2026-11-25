# Helper function: load bundled test data
load_toy_gene <- function() {
  read.csv(system.file("extdata", "mirage_toy.csv", package = "mirageV2"),
           stringsAsFactors = FALSE)
}

load_toy_vs <- function() {
  read.csv(system.file("extdata", "mirage_vs_toy.csv", package = "mirageV2"),
           stringsAsFactors = FALSE)
}

# Helper function: minimal synthetic data for fast tests
make_tiny_gene_data <- function() {
  data.frame(
    ID = paste0("v", 1:6),
    Gene = c("A", "A", "B", "B", "C", "C"),
    No.case = c(3L, 1L, 0L, 2L, 1L, 5L),
    No.contr = c(1L, 0L, 2L, 0L, 4L, 1L),
    category = c(1L, 2L, 1L, 2L, 1L, 2L),
    stringsAsFactors = FALSE
  )
}

make_tiny_vs_data <- function() {
  data.frame(
    ID = paste0("v", 1:6),
    No.case = c(3L, 1L, 0L, 2L, 1L, 5L),
    No.contr = c(1L, 0L, 2L, 0L, 4L, 1L),
    category = c(1L, 2L, 1L, 2L, 1L, 2L),
    stringsAsFactors = FALSE
  )
}
