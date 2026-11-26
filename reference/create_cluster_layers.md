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
  cluster_col = "cluster",
  center = NULL,
  core_radius = NULL
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

- center:

  Optional spot used as the center of the core: its identifier, or a
  one-row data.frame with its exact coordinates (e.g., \`A\` copied from
  \[view_spots()\]). It must be an existing non-border spot of
  \`target_cluster\` that does not touch the border layer either. If
  \`NULL\` (default), layers are defined by the distance to the border.

- core_radius:

  Optional numeric, only used with \`center\`. Maximum distance to
  \`center\` of core spots, in coordinate units. If \`NULL\` (default),
  it is set from \`intermediate_quantile\`.

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

Spots touching a border spot (kNN-based, in either direction) are always
intermediate, so that core and border are never in contact. Narrow
clusters may thus have no core.

If \`center\` is provided, border spots are unchanged, but core spots
are the non-border spots closest to \`center\`: those within
\`core_radius\` or, by default, the closest \`1 -
intermediate_quantile\` proportion.

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
  x = SpatialExperiment::spatialCoords(spe)[, 1],
  y = SpatialExperiment::spatialCoords(spe)[, 2],
  cluster = SummarizedExperiment::colData(spe)$cluster
)
# Default: moderate depth
layers <- create_cluster_layers(df, target_cluster = 1, k = 6)
table(layers$layer)
#> 
#>       border         core intermediate 
#>           34           20           27 

# Narrow intermediate layer
layers_narrow <- create_cluster_layers(df, target_cluster = 1, k = 6,
intermediate_quantile = 0.33)
table(layers_narrow$layer)
#> 
#>       border         core intermediate 
#>           34           24           23 

# Core centered on a chosen spot: here the core spot closest to the
# cluster centroid
inner <- layers[layers$layer == "core", ]
center <- closest_spot(inner, mean(inner$x), mean(inner$y))$spot_id
layers_center <- create_cluster_layers(df, target_cluster = 1, k = 6,
center = center)
table(layers_center$layer)
#> 
#>       border         core intermediate 
#>           34           20           27 
```
