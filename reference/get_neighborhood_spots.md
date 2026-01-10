# Get neighborhood spots around a target cluster or point

This function identifies all spots in the neighborhood of a target
cluster or a specific point using k-nearest neighbors search. It returns
all neighbor spots that were identified, marking which ones are
neighbors.

## Usage

``` r
get_neighborhood_spots(
  df,
  source = NULL,
  spot_id = NULL,
  k = 100,
  max_dist = NULL,
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

- source:

  The cluster label for which to analyze the neighborhood. Either
  \`source\` or \`spot_id\` must be provided, but not both.

- spot_id:

  Character. The spot identifier to find neighbors for. Either
  \`source\` or \`spot_id\` must be provided, but not both.

- k:

  Integer. Number of nearest neighbors to consider. Default is 100.
  Larger values capture more distant neighbors. If k exceeds the number
  of spots, it will be capped at nrow(df) - 1.

- max_dist:

  Numeric. Maximum distance from the target to consider. If \`NULL\`
  (default), all neighbors up to k are included without distance
  filtering.

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

- cluster:

  Cluster label of the neighbor spot.

- source:

  The source cluster or spot_id being analyzed.

- is_neighbourhood:

  Logical, always TRUE for returned rows.

Rows are sorted by cluster and spot_id.

## Details

The function can work in two modes: - \*\*Cluster mode\*\*: If
\`source\` is provided, computes k-nearest neighbors for all spots in
that cluster. Returns all neighbors found (excluding source cluster
spots). - \*\*Point mode\*\*: If \`spot_id\` is provided, finds the
k-nearest neighbors to specific spot.

## Examples

``` r
data("visium_simulated_spe", package = "Battlefield")
spe <- visium_simulated_spe
df <- data.frame(
  spot_id = colnames(spe),
  x = spatialCoords(spe)[, 1],
  y = spatialCoords(spe)[, 2],
  cluster = colData(spe)$cluster
)
#> Error in spatialCoords(spe): could not find function "spatialCoords"
# Get all neighbor spots for cluster 1
neighbors_cluster <- get_neighborhood_spots(df, source = 1, k = 50)
#> Error in get_neighborhood_spots(df, source = 1, k = 50): is.data.frame(df) is not TRUE
head(neighbors_cluster)
#> Error: object 'neighbors_cluster' not found

# Get neighbors for a specific point
first_spot <- df$spot_id[1]
#> Error in df$spot_id: object of type 'closure' is not subsettable
neighbors_point <- get_neighborhood_spots(df, spot_id = first_spot, k = 10)
#> Error in get_neighborhood_spots(df, spot_id = first_spot, k = 10): is.data.frame(df) is not TRUE
head(neighbors_point)
#> Error: object 'neighbors_point' not found
```
