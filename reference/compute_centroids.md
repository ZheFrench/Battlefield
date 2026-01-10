# Compute cluster centroids (mean x/y per cluster)

Computes the centroid of each cluster by taking the mean of the \`x\`
and \`y\` coordinates for each \`cluster\` level.

## Usage

``` r
compute_centroids(df)
```

## Arguments

- df:

  A \`data.frame\` containing at least the columns \`cluster\`, \`x\`,
  and \`y\`. \`cluster\` can be any type accepted by \`aggregate()\`
  grouping (e.g., integer, factor, character). \`x\` and \`y\` must be
  numeric.

## Value

A \`data.frame\` with one row per cluster and columns:

- cluster:

  Cluster identifier.

- x:

  Mean x-coordinate for the cluster.

- y:

  Mean y-coordinate for the cluster.

## Details

This is a thin wrapper around
[`aggregate`](https://rdrr.io/r/stats/aggregate.html) using
`FUN = mean`. Missing values are handled according to \`mean()\`'s
default behavior (i.e., \`NA\`s will propagate unless you pre-handle
them).

## Examples

``` r
df <- data.frame(
  cluster = c(1, 1, 2, 2),
  x = c(0, 2, 10, 12),
  y = c(1, 3, 5, 7)
)
compute_centroids(df)
#>   cluster  x y
#> 1       1  1 2
#> 2       2 11 6
```
