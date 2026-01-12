# Count neighborhood composition for all clusters

This function computes neighborhood statistics for all clusters in a
dataset in a single call. It returns a dataframe with counts of
neighboring cluster types for each source cluster.

## Usage

``` r
count_all_neighborhoods(
  df,
  clusters = NULL,
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
  unique identifier for each spot (optional, required for detailed
  analysis)

- clusters:

  Optional vector of cluster labels to process. If \`NULL\`, all unique
  clusters in \`df\` are used.

- k:

  Integer. Number of nearest neighbors to consider. Default is 100.
  Larger values capture more distant neighbors.

- max_dist:

  Numeric. Maximum distance from the target to consider. If \`NULL\`
  (default), all neighbors up to k are included without distance
  filtering.

- inlaid_col:

  Character. Name of the column containing inlaid/annotation values to
  count in neighborhoods. If \`NULL\` (default), uses the
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

- cluster:

  The source cluster being analyzed.

- neighborhood:

  The neighbor cluster type or inlaid value.

- count:

  Number of neighbor spots of this cluster type.

- proportion:

  Proportion of this cluster among all neighbors of the cluster.

Rows are grouped by source cluster and sorted by count within each group
(descending).

## Details

This is a batch processing function that applies
\[count_neighborhood()\] to all clusters and combines the results into a
single dataframe. Each row represents a neighbor cluster found around a
source cluster, with counts and proportions.

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
# Get neighborhood statistics for all clusters
all_neighbors <- count_all_neighborhoods(df, k = 50)
#> Error in count_all_neighborhoods(df, k = 50): is.data.frame(df) is not TRUE
head(all_neighbors)
#> Error: object 'all_neighbors' not found
```
