# Get inlaid spots within a source cluster

This function retrieves all spots from a source cluster along with their
inlaid/annotation values. Unlike \[get_neighborhood_spots()\], this
returns spots \*inside\* the source cluster itself, not around it.

## Usage

``` r
get_inlaid_spots(
  df,
  cluster,
  inlaid_col = "cluster",
  cluster_col = "cluster",
  coords = c("x", "y")
)
```

## Arguments

- df:

  A data.frame containing at least the columns \`x\`, \`y\`, cluster
  identifier column, inlaid column, and \`spot_id\`. - \`x\`, \`y\`:
  numeric coordinates of spots - cluster column (name specified by
  \`cluster_col\`): cluster assignment for each spot - inlaid column
  (name specified by \`inlaid_col\`): inlaid/annotation for each spot -
  \`spot_id\`: unique identifier for each spot

- cluster:

  The cluster label to retrieve inlaid spots for.

- inlaid_col:

  Character. Name of the column containing inlaid/annotation values.
  Default is \`"cluster"\`.

- cluster_col:

  Character. Name of the column containing cluster assignments. Default
  is \`"cluster"\`.

- coords:

  Character vector of length 2 giving the coordinate column names.
  Default is \`c("x", "y")\`.

## Value

A data.frame with columns:

- spot_id:

  Unique spot identifier.

- x:

  X coordinate.

- y:

  Y coordinate.

- inlaid:

  Inlaid/annotation value.

- cluster:

  The source cluster being analyzed.

- is_inlaid:

  Logical, always TRUE for returned rows.

Rows are sorted by inlaid and spot_id.

## Examples

``` r
data("visium_simulated_spe", package = "Battlefield")
spe <- visium_simulated_spe
df <- data.frame(
  spot_id = colnames(spe),
  x = spatialCoords(spe)[, 1],
  y = spatialCoords(spe)[, 2],
  cluster = colData(spe)$cluster,
inlaid = sample(paste0("type_", 1:3), length(colnames(spe)), replace = TRUE)
)
#> Error in spatialCoords(spe): could not find function "spatialCoords"
# Get all inlaid spots within cluster 1
inlaid_spots <- get_inlaid_spots(df, cluster = 1, inlaid_col = "inlaid")
#> Error in get_inlaid_spots(df, cluster = 1, inlaid_col = "inlaid"): is.data.frame(df) is not TRUE
head(inlaid_spots)
#> Error: object 'inlaid_spots' not found
```
