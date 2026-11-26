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
  x = SpatialExperiment::spatialCoords(spe)[, 1],
  y = SpatialExperiment::spatialCoords(spe)[, 2],
  cluster = SummarizedExperiment::colData(spe)$cluster
)
all_layers <- create_all_layers(df, k = 6, intermediate_quantile = 0.5)
#> Cluster 2: No 'core' layer, all inner spots touch the border.
head(all_layers)
#>     spot_id x    y cluster layer
#> 1   bin16_1 0    0       3  core
#> 2  bin16_41 0  400       3  core
#> 3  bin16_81 0  800       3  core
#> 4 bin16_121 0 1200       3  core
#> 5 bin16_161 0 1600       3  core
#> 6 bin16_201 0 2000       3  core
```
