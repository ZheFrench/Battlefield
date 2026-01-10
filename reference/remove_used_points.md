# Remove rows whose rounded (x, y) coordinates have already been used

Filters \`df\` to drop any rows whose rounded coordinate key (as
produced by \`.xy_key(x, y)\`) appears in \`used_df\`.

## Usage

``` r
remove_used_points(df, used_df)
```

## Arguments

- df:

  A data frame containing at least columns \`x\` and \`y\`.

- used_df:

  A data frame containing at least columns \`x\` and \`y\` defining the
  set of already-used points to exclude.

## Value

A filtered data frame with the same columns as \`df\`, excluding rows
whose rounded \`(x, y)\` matches any row in \`used_df\`.

## Details

Matching is performed on rounded coordinates, not exact floating-point
values. Internally, keys are computed as \`paste0(round(x), "\_",
round(y))\`.

This function uses \`dplyr::filter()\` (and the \`\|\>\` pipe), so
\`dplyr\` must be installed and a pipe operator must be available.

## Examples

``` r
df <- data.frame(x = c(0.1, 1.2, 1.6, 2.0),
                 y = c(0.2, 3.4, 3.4, 9.5),
                 id = 1:4)
used_df <- data.frame(x = c(1.49, 2.1),
                      y = c(3.49, 9.49))

# Removes rows rounding to (1,3) and (2,9)
remove_used_points(df, used_df)
#>     x   y id
#> 1 0.1 0.2  1
#> 2 1.6 3.4  3
#> 3 2.0 9.5  4
```
