# Add border or core spot selections to SpatialExperiment colData

This function takes spot selection results (either border or core spots)
and adds them as annotations in the colData of a SpatialExperiment
object.

## Usage

``` r
add_borders_to_spe(spe, border = NULL, core = NULL, erase = FALSE)
```

## Arguments

- spe:

  A SpatialExperiment object from which the selection dataframes were
  derived.

- border:

  Optional data.frame of border spots. Can be output from either
  \[select_border_spots()\] or \[build_all_borders()\]. Expected
  columns: spot_id, interface, mode.

- core:

  Optional data.frame of core spots. Can be output from either
  \[select_core_spots()\] or \[build_all_cores()\]. Expected columns:
  spot_id, interface, mode.

- erase:

  Logical. If \`TRUE\`, erase pre-existing battlefield columns before
  adding new ones. If \`FALSE\` (default), issue a warning if columns
  already exist and skip adding them.

## Value

The SpatialExperiment object with updated colData.

## Details

At least \`border\` must be provided. \`core\` is optional. If \`core\`
is provided without \`border\`, an error is raised.

The following unified columns are added: - \`is_border\`: logical, TRUE
if spot is a border spot, FALSE or NA otherwise - \`is_core\`: logical,
TRUE if spot is a core spot, FALSE or NA otherwise - \`interface\`: the
target interface cluster - \`border_mode\`: character, "inner", "outer",
or NA (from the \`mode\` column in input data)

## Examples

``` r
data("visiumHD_16um_simulated_spe", package = "Battlefield")
spe <- visiumHD_16um_simulated_spe
#> Loading required namespace: SpatialExperiment
df <- data.frame(
  spot_id = colnames(spe),
  x = spatialCoords(spe)[, 1],
  y = spatialCoords(spe)[, 2],
  cluster = colData(spe)$cluster
)
#> Error in spatialCoords(spe): could not find function "spatialCoords"
# Using build_all_borders and build_all_cores
all_borders <- build_all_borders(df, k = 6)
#> Error in df[[cluster_col]]: object of type 'closure' is not subsettable
all_cores <- build_all_cores(df, all_borders, region = "inner")
#> Error in build_all_cores(df, all_borders, region = "inner"): unused argument (region = "inner")

# Or using individual select functions
border_3_to_4 <- select_border_spots(df, cluster = 3, interface = 4, k = 6)
#> Error in select_border_spots(df, cluster = 3, interface = 4, k = 6): all(coord_cols %in% colnames(df)) is not TRUE
core_3_to_4 <- select_core_spots(df, all_borders, cluster = 3, interface = 4)
#> Error in select_core_spots(df, all_borders, cluster = 3, interface = 4): all(coord_cols %in% colnames(df)) is not TRUE
```
