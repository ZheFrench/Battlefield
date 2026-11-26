# Remove the orphan spots of one cluster

Removes the spots of the cluster that belong to a contiguous piece of
fewer than \`min_size\` spots. Spots of other clusters are always kept.

## Usage

``` r
clean_orphans(
  df,
  clusters,
  min_size = 20,
  max_dist = NULL,
  coord_cols = c("x", "y"),
  cluster_col = "cluster"
)
```

## Arguments

- df:

  A data.frame with coordinates and a cluster column.

- clusters:

  The cluster to split (a single label).

- min_size:

  Pieces with fewer spots are removed.

- max_dist:

  Maximal distance between touching spots. \`NULL\` uses 1.5 times the
  spot spacing (direct grid neighbors).

- coord_cols, cluster_col:

  Names of the coordinate and cluster columns.

## Value

\`df\` without the orphan spots. If every piece is smaller than
\`min_size\`, \`df\` is returned unchanged with a warning.

## Details

To clean border spots, use it on the output of
\[create_cluster_layers()\]: layers are then computed with all spots,
and labels are not changed.

## See also

\[split_clusters()\]

## Examples

``` r
df <- expand.grid(x = 1:20, y = 1:10)
df$spot_id <- paste0("s", seq_len(nrow(df)))
df$cluster <- "B"
df$cluster[df$x <= 4 | df$x >= 18] <- "A"   # two pieces
df$cluster[df$x == 10 & df$y == 5] <- "A"   # one orphan spot
cleaned <- clean_orphans(df, "A", min_size = 5)
nrow(df) - nrow(cleaned)
#> [1] 1

# Clean layers: computed with all spots, then orphans are removed
layers <- create_cluster_layers(df, "A", k = 4)
layers <- clean_orphans(layers, "A", min_size = 5)
table(layers$layer)
#> 
#>       border         core intermediate 
#>           20           10           40 
```
