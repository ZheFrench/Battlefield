# Find the closest spot to a target point

Returns the row of \`df\` corresponding to the spot nearest to the
target coordinates \`(tx, ty)\` using squared Euclidean distance.

## Usage

``` r
closest_spot(df, tx, ty)
```

## Arguments

- df:

  A data frame containing at least numeric columns \`x\` and \`y\`
  representing spot coordinates.

- tx:

  Numeric scalar; target x-coordinate.

- ty:

  Numeric scalar; target y-coordinate.

## Value

A one-row data frame (same columns as \`df\`) corresponding to the
closest spot. The result is always a data frame (\`drop = FALSE\`).

## Details

The function minimizes \`(x - tx)^2 + (y - ty)^2\`. Squared distances
are used for efficiency; taking the square root is unnecessary for
argmin.

If multiple spots are tied for the minimum distance, the first
occurrence is returned (as per \`which.min()\`).

## Examples

``` r
df <- data.frame(x = c(0, 2, 5), y = c(0, 2, 1), id = c("a", "b", "c"))
closest_spot(df, tx = 1.9, ty = 2.1)
#>   x y id
#> 2 2 2  b
```
