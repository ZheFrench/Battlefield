# Count inlaid composition within a source cluster

This function counts the composition of an inlaid/annotation column
within a target source cluster. Unlike \[count_neighborhood()\], this
counts types \*inside\* the source cluster itself, not around it.

## Usage

``` r
count_inlaid(df, cluster, inlaid_col = "cluster", cluster_col = "cluster")
```

## Arguments

- df:

  A data.frame containing at least the columns \`x\`, \`y\`, cluster
  identifier column, and \`spot_id\`. - \`x\`, \`y\`: numeric
  coordinates of spots - cluster column (name specified by
  \`cluster_col\`): cluster assignment for each spot - inlaid column
  (name specified by \`inlaid_col\`): inlaid/annotation for each spot -
  \`spot_id\`: unique identifier for each spot

- cluster:

  The cluster label to analyze inlaid composition for.

- inlaid_col:

  Character. Name of the column containing inlaid/annotation values.
  Default is \`"cluster"\` (to analyze cluster composition within
  another grouping).

- cluster_col:

  Character. Name of the column containing cluster assignments. Default
  is \`"cluster"\`.

## Value

A data.frame with columns:

- inlaid:

  Inlaid/annotation value.

- count:

  Number of spots with this inlaid type in the cluster.

- proportion:

  Proportion of this inlaid type among all spots in cluster.

Rows are sorted by count in descending order.

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
# Count inlaid types within cluster 1
inlaid_counts <- count_inlaid(df, cluster = 1, inlaid_col = "inlaid")
inlaid_counts
#>   cluster inlaid count proportion
#> 1       1 type_3    29  0.3580247
#> 2       1 type_1    28  0.3456790
#> 3       1 type_2    24  0.2962963
```
