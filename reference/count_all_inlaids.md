# Count inlaid composition for all clusters

This function counts inlaid/annotation composition within each cluster
in a single call. It returns a dataframe with counts of inlaid types for
each source cluster.

## Usage

``` r
count_all_inlaids(
  df,
  sources = NULL,
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

- sources:

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

- source:

  The source cluster being analyzed.

- inlaid:

  The inlaid/annotation type.

- count:

  Number of spots with this inlaid type in the source.

- proportion:

  Proportion of this inlaid type among all spots in source.

Rows are grouped by source cluster and sorted by count within each group
(descending).

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
# Get inlaid statistics for all clusters
all_inlaids <- count_all_inlaids(df, inlaid_col = "inlaid")
#> Error in count_all_inlaids(df, inlaid_col = "inlaid"): is.data.frame(df) is not TRUE
head(all_inlaids)
#> Error: object 'all_inlaids' not found
```
