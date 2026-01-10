# Add layer classifications to SpatialExperiment colData

This function takes layer classification results (border, intermediate,
core) and adds them as annotations in the colData of a SpatialExperiment
object.

## Usage

``` r
add_layers_to_spe(spe, layer = NULL, erase = FALSE)
```

## Arguments

- spe:

  A SpatialExperiment object from which the layer dataframes were
  derived.

- layer:

  A data.frame of layer classifications. Can be output from either
  \[create_cluster_layers()\] or \[create_all_layers()\]. Expected
  columns: spot_id, layer.

- erase:

  Logical. If \`TRUE\`, erase pre-existing layer columns before adding
  new ones. If \`FALSE\` (default), issue a warning if columns already
  exist and skip adding them.

## Value

The SpatialExperiment object with updated colData.

## Details

The \`layer\` parameter must be provided and contain at least one row.

The following column is added: - \`layer\`: character, one of "border",
"intermediate", "core", or NA

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
# Or for a single cluster
cluster_1_layers <- create_cluster_layers(df, target_cluster = 1, k = 6)
#> Error in create_cluster_layers(df, target_cluster = 1, k = 6): all(coord_cols %in% colnames(df)) is not TRUE
spe <- add_layers_to_spe(spe, layer = cluster_1_layers)
#> Error: object 'cluster_1_layers' not found
```
