# Tidy method for mirage_result

Returns gene-level results as a tidy data frame.

## Usage

``` r
# S3 method for class 'mirage_result'
tidy(x, ...)
```

## Arguments

- x:

  A `mirage_result` object.

- ...:

  Ignored.

## Value

A data frame with columns: gene, bf, post.prob, p.value.
