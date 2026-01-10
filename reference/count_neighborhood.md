# Count annotated spot types in the neighborhood of a cluster or point

This function computes summary statistics of cluster types in the
neighborhood of a target cluster or a specific point. It uses
\[get_neighborhood_spots()\] internally to identify neighbors, then
counts how many belong to each cluster type.

## Usage

``` r
count_neighborhood(
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
  unique identifier for each spot (optional, required for detailed
  analysis)

- source:

  The cluster label for which to analyze the neighborhood. Either
  \`source\` or \`spot_id\` must be provided, but not both.

- spot_id:

  Character. The spot identifier to find neighbors for. Either
  \`source\` or \`spot_id\` must be provided, but not both.

- k:

  Integer. Number of nearest neighbors to consider. Default is 100.
  Larger values capture more distant neighbors.

- max_dist:

  Numeric. Maximum distance from the target to consider. If \`NULL\`
  (default), all neighbors up to k are included without distance \#'
  filtering.

- coords:

  Character vector of length 2 giving the coordinate column names.
  Default is \`c("x", "y")\`.

- cluster_col:

  Character. Name of the column containing cluster assignments. Default
  is \`"cluster"\`.

## Value

A data.frame with columns:

- cluster:

  Cluster label (from neighboring spots).

- count:

  Number of times this cluster appears in the neighborhood.

- proportion:

  Proportion of this cluster relative to all neighbors.

Rows are sorted by count in descending order.

## Details

This function wraps \[get_neighborhood_spots()\] and summarizes the
results by counting spots from each neighboring cluster. In cluster
mode, the source cluster itself is excluded from the counts. In point
mode, all neighboring clusters are counted.

The result shows the composition of the neighborhood: how many spots of
each cluster type surround the target.

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
# Count cluster types in neighborhood of cluster 1
neighbor_counts <- count_neighborhood(df, source = 1, k = 50)
#> Error in get_neighborhood_spots(df, source, spot_id, k, max_dist, coords,     cluster_col): is.data.frame(df) is not TRUE
neighbor_counts
#> Error: object 'neighbor_counts' not found

# Count cluster types around a specific point
first_spot <- df$spot_id[1]
#> Error in df$spot_id: object of type 'closure' is not subsettable
point_neighbors <- count_neighborhood(df, spot_id = first_spot, k = 10)
#> Error in get_neighborhood_spots(df, source, spot_id, k, max_dist, coords,     cluster_col): is.data.frame(df) is not TRUE
point_neighbors
#> Error: object 'point_neighbors' not found
```
