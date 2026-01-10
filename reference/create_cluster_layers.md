# Create layer classification for spots in a cluster

This function classifies spots within a target cluster into three
layers: core, intermediate, and border. The border layer consists of
spots at the cluster interface (touching other clusters). The core layer
consists of spots far from any border. The intermediate layer bridges
the two, with depth defined by proximity to the border.

## Usage

``` r
create_cluster_layers(
  df,
  target_cluster,
  k = 6,
  max_dist = NULL,
  intermediate_quantile = 0.5,
  coord_cols = c("x", "y"),
  cluster_col = "cluster"
)
```

## Arguments

- df:

  A data.frame containing at least the coordinate columns and a cluster
  label column.

- target_cluster:

  Cluster label for which to create layers.

- k:

  Integer. Number of nearest neighbors to consider (excluding self) when
  identifying border spots. Default is 6.

- max_dist:

  Optional numeric. Maximum Euclidean distance for neighbors. If
  \`NULL\`, no distance filtering is applied.

- intermediate_quantile:

  Numeric between 0 and 1. The quantile of distances to border used to
  define the intermediate layer threshold. Default is 0.5 (median). \#'
  Spots with distance to border \<= this quantile threshold are
  classified as intermediate. Use lower values (e.g., 0.33) for narrower
  intermediate layers, higher values (e.g., 0.75) for wider intermediate
  layers.

- coord_cols:

  Character vector of length 2 giving the coordinate column names.
  Default is \`c("x","y")\`.

- cluster_col:

  Character. Name of the column containing cluster labels. Default is
  \`"cluster"\`.

## Value

A data.frame based on \`df\` filtered to contain only spots from
\`target_cluster\`, with an additional column: - \`layer\`: character,
one of "border", "intermediate", or "core"

## Details

Layer definitions: - \*\*Border\*\*: Spots in the target cluster that
touch at least one spot from a different cluster (kNN-based). -
\*\*Intermediate\*\*: Non-border spots whose distance to the nearest
border spot is within the \`intermediate_quantile\` threshold of all
non-border spots. - \*\*Core\*\*: All remaining non-border spots that
are far from the border.

The \`intermediate_quantile\` parameter controls the depth of the
intermediate layer: - 0.33: Narrow intermediate layer (only very close
to border) - 0.50: Moderate intermediate layer (median distance
threshold) - 0.75: Wide intermediate layer (most non-border spots
included)

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
# Default: moderate depth
layers <- create_cluster_layers(df, target_cluster = 1, k = 6)
#> Error in create_cluster_layers(df, target_cluster = 1, k = 6): all(coord_cols %in% colnames(df)) is not TRUE
table(layers$layer)
#> Error: object 'layers' not found

# Narrow intermediate layer
layers_narrow <- create_cluster_layers(df, target_cluster = 1, k = 6,
intermediate_quantile = 0.33)
#> Error in create_cluster_layers(df, target_cluster = 1, k = 6, intermediate_quantile = 0.33): all(coord_cols %in% colnames(df)) is not TRUE
table(layers_narrow$layer)
#> Error: object 'layers_narrow' not found
```
