# Get neighborhood spots around a target cluster or point

This function identifies all spots in the neighborhood of a target
cluster or a specific point using k-nearest neighbors search. It returns
all neighbor spots that were identified, marking which ones are
neighbors.

## Usage

``` r
get_neighborhood_spots(
  df,
  cluster = NULL,
  spot_id = NULL,
  k = 100,
  max_dist = NULL,
  inlaid_col = NULL,
  coords = c("x", "y"),
  cluster_col = "cluster"
)
```

## Arguments

- df:

  A data.frame containing at least the columns \`x\`, \`y\`, cluster
  identifier column, and \`spot_id\`. - \`x\`, \`y\`: numeric
  coordinates of spots - cluster column (name specified by
  \`cluster_col\`): cluster assignment for each spot - \`spot_id\`:
  unique identifier for each spot

- cluster:

  The cluster label for which to analyze the neighborhood. Either
  \`cluster\` or \`spot_id\` must be provided, but not both.

- spot_id:

  The spot to find neighbors for: its identifier, or a one-row
  data.frame with its exact coordinates (e.g., \`A\` copied from
  \[view_spots()\]). Coordinates must be those of an existing spot.
  Either \`cluster\` or \`spot_id\` must be provided, but not both.

- k:

  Integer. Number of nearest neighbors to consider. Default is 100.
  Larger values capture more distant neighbors. If k exceeds the number
  of spots, it will be capped at nrow(df) - 1.

- max_dist:

  Numeric. Maximum distance from the target to consider. If \`NULL\`
  (default), all neighbors up to k are included without distance
  filtering.

- inlaid_col:

  Character. Name of the column containing inlaid/annotation values to
  return in the \`neighborhood\` column. If \`NULL\` (default), uses the
  \`cluster_col\` value. This allows counting neighborhood composition
  by any column (e.g., cell type, tissue type).

- coords:

  Character vector of length 2 giving the coordinate column names.
  Default is \`c("x", "y")\`.

- cluster_col:

  Character. Name of the column containing cluster assignments. Default
  is \`"cluster"\`.

## Value

A data.frame with columns:

- spot_id:

  Unique spot identifier.

- x:

  X coordinate.

- y:

  Y coordinate.

- neighborhood:

  Value from \`inlaid_col\` (or cluster if inlaid_col is NULL) of the
  neighbor spot.

- cluster:

  The source cluster or spot_id being analyzed.

- is_neighborhood:

  Logical, always TRUE for returned rows.

Rows are sorted by neighborhood and spot_id.

## Details

The function can work in two modes: - \*\*Cluster mode\*\*: If
\`cluster\` is provided, computes k-nearest neighbors for all spots in
that cluster. Returns all neighbors found (excluding source cluster
spots). - \*\*Point mode\*\*: If \`spot_id\` is provided, finds the
k-nearest neighbors to specific spot.

## Examples

``` r
data("visium_simulated_spe", package = "Battlefield")
spe <- visium_simulated_spe
df <- data.frame(
  spot_id = colnames(spe),
  x = SpatialExperiment::spatialCoords(spe)[, 1],
  y = SpatialExperiment::spatialCoords(spe)[, 2],
  cluster = SummarizedExperiment::colData(spe)$cluster
)
# Get all neighbor spots for cluster 1
neighbors_cluster <- get_neighborhood_spots(df, cluster = 1, k = 50)
head(neighbors_cluster)
#>           spot_id      x        y cluster neighborhood is_neighborhood
#> SPOT0465 SPOT0465 1347.5 523.9454       1            2            TRUE
#> SPOT0466 SPOT0466 1402.5 523.9454       1            2            TRUE
#> SPOT0506 SPOT0506 1375.0 571.5768       1            2            TRUE
#> SPOT0545 SPOT0545 1347.5 619.2082       1            2            TRUE
#> SPOT0546 SPOT0546 1402.5 619.2082       1            2            TRUE
#> SPOT0219 SPOT0219 1017.5 238.1570       1            3            TRUE

# Get neighbors for a specific point
first_spot <- df$spot_id[1]
neighbors_point <- get_neighborhood_spots(df, spot_id = first_spot, k = 10)
head(neighbors_point)
#>           spot_id     x        y  cluster neighborhood is_neighborhood
#> SPOT0002 SPOT0002  55.0  0.00000 SPOT0001            3            TRUE
#> SPOT0003 SPOT0003 110.0  0.00000 SPOT0001            3            TRUE
#> SPOT0041 SPOT0041  27.5 47.63140 SPOT0001            3            TRUE
#> SPOT0042 SPOT0042  82.5 47.63140 SPOT0001            3            TRUE
#> SPOT0043 SPOT0043 137.5 47.63140 SPOT0001            3            TRUE
#> SPOT0081 SPOT0081   0.0 95.26279 SPOT0001            3            TRUE
```
