# Validate and Clean MIRAGE Input Data

Validate and Clean MIRAGE Input Data

## Usage

``` r
validate_mirage_data(data, type = c("gene", "vs"))
```

## Arguments

- data:

  A data frame containing variant counts.

- type:

  Either "gene" or "vs".

## Value

A cleaned data frame with S3 class `mirage_data`.
