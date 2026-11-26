# Add split cluster labels to SpatialExperiment colData

This function adds the labels of a split cluster (e.g., \`C2_1\`,
\`C2_2\`) to the colData of a SpatialExperiment object. Labels are
matched by \`spot_id\`, whatever the order of the rows. The original
cluster column is not changed.

## Usage

``` r
add_split_to_spe(spe, split, split_col = "cluster_split", erase = FALSE)
```

## Arguments

- spe:

  A SpatialExperiment object from which the split data.frame was
  derived.

- split:

  A data.frame with one row per spot and columns \`spot_id\` and
  \`split_col\`, typically the output of \[split_clusters()\].

- split_col:

  Character. Name of the column of \`split\` holding the labels, also
  used as the name of the new colData column. Default is
  \`"cluster_split"\`.

- erase:

  Logical. If \`TRUE\`, replace a pre-existing column of that name. If
  \`FALSE\` (default), the object is returned unchanged with a message
  if the column already exists.

## Value

The SpatialExperiment object with updated colData.

## Details

Spots of \`spe\` that are absent from \`split\` get \`NA\`. It is an
error if a spot of \`split\` is absent from \`spe\`, or if identifiers
are duplicated.

Spot identifiers of \`spe\` are taken from the \`spot_id\` column of
\`colData\` if present, otherwise from \`colnames(spe)\`.

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
split <- split_clusters(df, "1", min_size = 5)
spe <- add_split_to_spe(spe, split)
table(spe$cluster, spe$cluster_split)
#>    
#>       1   2   3   4   5
#>   1  81   0   0   0   0
#>   2   0  47   0   0   0
#>   3   0   0 428   0   0
#>   4   0   0   0 411   0
#>   5   0   0   0   0  33
```
