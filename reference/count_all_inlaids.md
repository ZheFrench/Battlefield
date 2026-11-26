# Count inlaid composition for all clusters

This function counts inlaid/annotation composition within each cluster
in a single call. It returns a dataframe with counts of inlaid types for
each source cluster.

## Usage

``` r
count_all_inlaids(
  df,
  clusters = NULL,
  inlaid_col = "cluster",
  cluster_col = "cluster"
)
```

## Arguments

- df:

  A data.frame containing at least the columns \`x\`, \`y\`, cluster
  identifier column, and \`spot_id\`. - \`x\`, \`y\`: numeric
  coordinates of spots - cluster column (name specified by
  \`cluster_col\`): cluster assignment for \#' each spot - inlaid column
  (name specified by \`inlaid_col\`): inlaid/annotation for each spot -
  \`spot_id\`: unique identifier for each spot

- clusters:

  Optional vector of cluster labels to process. If \`NULL\`, all unique
  clusters in \`df\` are used.

- inlaid_col:

  Character. Name of the column containing inlaid/annotation values.
  Default is \`"cluster"\`.

- cluster_col:

  Character. Name of the column containing cluster assignments. Default
  is \`"cluster"\`.

## Value

A data.frame with columns:

- cluster:

  The source cluster being analyzed.

- inlaid:

  The inlaid/annotation type.

- count:

  Number of spots with this inlaid type in the cluster.

- proportion:

  Proportion of this inlaid type among all spots in cluster.

Rows are grouped by source cluster and sorted by count within each group
(descending).

## Examples

``` r
data("visium_simulated_spe", package = "Battlefield")
spe <- visium_simulated_spe
df <- data.frame(
  spot_id = colnames(spe),
  x = SpatialExperiment::spatialCoords(spe)[, 1],
  y = SpatialExperiment::spatialCoords(spe)[, 2],
  cluster = SummarizedExperiment::colData(spe)$cluster,
inlaid = sample(paste0("type_", 1:3), length(colnames(spe)), replace = TRUE)
)
# Get inlaid statistics for all clusters
all_inlaids <- count_all_inlaids(df, inlaid_col = "inlaid")
head(all_inlaids)
#>   cluster inlaid count proportion
#> 1       1 type_1    30  0.3703704
#> 2       1 type_2    26  0.3209877
#> 3       1 type_3    25  0.3086420
#> 4       2 type_1    19  0.4042553
#> 5       2 type_3    15  0.3191489
#> 6       2 type_2    13  0.2765957
```
