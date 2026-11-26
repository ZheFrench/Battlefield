# Select border spots from cluster A that touch cluster B (and flag junctions)

This function identifies \*\*directed interface spots\*\*: spots
belonging to \`cluster\` (A) that have at least one neighbor in
\`interface\` (B) within a k-nearest-neighbors neighborhood (optionally
constrained by a distance cutoff).

## Usage

``` r
select_border_spots(
  df,
  cluster,
  interface,
  mode = "inner",
  k = 7,
  max_dist = NULL,
  coord_cols = c("x", "y"),
  cluster_col = "cluster"
)
```

## Arguments

- df:

  A data.frame containing at least the coordinate columns and a cluster
  label column.

- cluster:

  Cluster label for the source cluster (A). Spots in this cluster are
  tested for adjacency to \`interface\`.

- interface:

  Cluster label for the target cluster (B).

- mode:

  Character. One of "inner", "outer", or "both". - "inner": returns
  border spots from cluster → interface (default). - "outer": returns
  border spots from interface → cluster (swaps cluster/interface). -
  "both": returns border spots from both directions, row-bound together.
  Default is "inner".

- k:

  Integer. Number of nearest neighbors to consider (excluding self).
  Internally uses \`k + 1\` to include self then removes it. Default is
  7.

- max_dist:

  Optional numeric. Maximum Euclidean distance for a neighbor to be
  considered (distance cutoff). If \`NULL\`, no distance filtering is
  applied.

- coord_cols:

  Character vector of length 2 giving the coordinate column names.
  Default is \`c("x","y")\`.

- cluster_col:

  Character. Name of the column containing cluster labels. Default is
  \`"cluster"\`.

## Value

A data.frame subset of \`df\` containing border spots from \`cluster\`
to \`interface\` with columns: spot_id, x, y, cluster, interface,
directed_pair, is_border, is_border_multiple, other_adjacent_borders.

## Details

In addition, it flags spots that also touch \*\*other clusters\*\*
(i.e., junction / multi-interface spots) and reports which other
clusters are touched.

Steps:

1.  Builds a kNN neighborhood for each spot in \`cluster\`.

2.  Optionally removes neighbors farther than \`max_dist\`.

3.  Marks a spot as a border spot if it has at least one neighbor in
    \`interface\`.

4.  Flags a spot as multi-interface if it also has a neighbor belonging
    to a cluster different from \`cluster\` and \`interface\`.

\*\*Note on mode column\*\*: The \`mode\` column in the returned
data.frame will always contain either "inner" or "outer", even if the
\`mode\` parameter is set to "both". When \`mode="both"\`, the function
internally calls the selection algorithm twice (once for
cluster→interface and once for interface→cluster), then row-binds the
results, so each row is labeled with its actual direction.

The returned data.frame contains only spots from \`cluster\` that touch
\`interface\`, with columns:

- spot_id:

  Spot identifier.

- x:

  X coordinate.

- y:

  Y coordinate.

- cluster:

  Cluster label (same as \`cluster\`).

- interface:

  The value of \`interface\`.

- directed_pair:

  Interface label combining cluster and \`interface\`, e.g. \`"A-B"\`.

- undirected_pair:

  Undirected pair label (both directions), e.g. \`"A-B / B-A"\`.

- is_border:

  Always \`TRUE\` for returned rows (kept for clarity).

- is_border_multiple:

  \`TRUE\` if the spot also touches other clusters.

- other_adjacent_borders:

  Comma-separated list of other clusters touched, or \`NA\`.

- mode:

  The mode used for selection: "inner", "outer", or "both".

## Examples

``` r
# Minimal example (requires RANN and dplyr)
set.seed(1)
df_ex <- data.frame(
  x = rnorm(200),
  y = rnorm(200),
  cluster = sample(c("A","B","C"), 200, replace = TRUE)
)
res <- select_border_spots(df_ex, cluster = "A", interface = "B", k = 7)
head(res)
#>             x           y directed_pair undirected_pair cluster interface
#> 2   0.1836433  1.68887329           A-B       A-B / B-A       A         B
#> 3  -0.8356286  1.58658843           A-B       A-B / B-A       A         B
#> 5   0.3295078 -2.28523554           A-B       A-B / B-A       A         B
#> 9   0.5757814 -0.01339952           A-B       A-B / B-A       A         B
#> 10 -0.3053884  0.51010842           A-B       A-B / B-A       A         B
#> 11  1.5117812 -0.16437583           A-B       A-B / B-A       A         B
#>    is_border is_border_multiple other_adjacent_borders  mode
#> 2       TRUE               TRUE                      C inner
#> 3       TRUE               TRUE                      C inner
#> 5       TRUE               TRUE                      C inner
#> 9       TRUE               TRUE                      C inner
#> 10      TRUE               TRUE                      C inner
#> 11      TRUE               TRUE                      C inner
```
