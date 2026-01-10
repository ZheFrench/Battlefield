# Create and visualize layers for multiple clusters

Batch create layer classifications for multiple clusters and optionally
generate individual plots for each cluster.

## Usage

``` r
create_all_layers(
  df,
  clusters = NULL,
  k = 6,
  max_dist = NULL,
  intermediate_quantile = 0.5,
  coord_cols = c("x", "y"),
  cluster_col = "cluster"
)
```

## Arguments

- df:

  A data.frame containing at least coordinate and cluster columns.

- clusters:

  Optional vector of cluster labels to process. If \`NULL\`, all unique
  clusters are used.

- k:

  Integer. Number of nearest neighbors to consider. Default is 6.

- max_dist:

  Optional numeric. Maximum Euclidean distance for neighbors.

- intermediate_quantile:

  Numeric between 0 and 1. The quantile of distances to border used to
  define the intermediate layer threshold. Default is 0.5 (median).

- coord_cols:

  Character vector of length 2 giving the coordinate column names.
  Default is \`c("x","y")\`.

- cluster_col:

  Character. Name of the column containing cluster labels. Default is
  \`"cluster"\`.

## Value

A data.frame combining all processed clusters with the \`layer\` column
added.

## Examples

``` r
data("visiumHD_16um_simulated_spe", package = "Battlefield")
spe <- visiumHD_16um_simulated_spe
df <- data.frame(
  spot_id = colnames(spe),
  x = spatialCoords(spe)[, 1],
  y = spatialCoords(spe)[, 2],
  cluster = colData(spe)$cluster
)
#> Error in spatialCoords(spe): could not find function "spatialCoords"
all_layers <- create_all_layers(df, k = 6, intermediate_quantile = 0.5)
#> Error in create_all_layers(df, k = 6, intermediate_quantile = 0.5): all(coord_cols %in% colnames(df)) is not TRUE
head(all_layers)
#> Error: object 'all_layers' not found
```
