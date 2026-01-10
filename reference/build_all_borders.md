# Build border spots for all oriented cluster pairs

This function iterates over all \*\*ordered\*\* cluster pairs (A -\> B)
and returns a single data.frame containing the border spots for each
interface, as computed by a border-detection function (e.g.
\[select_border_spots()\]).

## Usage

``` r
build_all_borders(
  df,
  k = 6,
  max_dist = NULL,
  mode = "both",
  pairs = NULL,
  coord_cols = c("x", "y"),
  cluster_col = "cluster"
)
```

## Arguments

- df:

  A data.frame containing at least coordinate columns and a cluster
  label column.

- k:

  Integer. Number of nearest neighbors to consider (excluding self) when
  computing borders. Default is 6.

- max_dist:

  Optional numeric. Maximum Euclidean distance for neighbors to be
  considered. If \`NULL\`, no distance filtering is applied.

- mode:

  Character. One of "inner", "outer", or "both". Controls which border
  direction(s) to select for each pair. Default is "both".

- pairs:

  Optional data.frame of oriented cluster pairs, typically produced by
  \[directed_cluster_interface_pairs()\]. Must contain \`cluster\` and
  \`interface\` columns. If \`NULL\`, it is computed from
  \`df\[\[cluster_col\]\]\`.

- coord_cols:

  Character vector of length 2 giving the coordinate column names.
  Default is \`c("x","y")\`.

- cluster_col:

  Character. Name of the column containing cluster labels. Default is
  \`"cluster"\`.

## Value

A data.frame produced by row-binding the result of border selection for
each oriented pair. Typically contains the original columns of \`df\`
plus interface annotation columns from the border selector.

\*\*Note on mode column\*\*: The \`mode\` column in the returned
data.frame will always contain either "inner" or "outer", never "both".
Even when the \`mode\` parameter is set to "both", each border spot row
is labeled with its actual direction (cluster→interface is "inner",
interface→cluster is "outer").

## Details

If \`pairs_df\` is not provided, it is generated from the cluster labels
using \[directed_cluster_interface_pairs()\].

## Examples

``` r
# Example with synthetic data
set.seed(1)
df_ex <- data.frame(
  x = rnorm(200),
  y = rnorm(200),
  cluster = sample(c("A","B","C"), 200, replace = TRUE)
)
all_borders <- build_all_borders(df_ex, k = 6)
#> Error in dplyr::select(dplyr::mutate(out, interface = i_val, directed_pair = paste0(c_val,     "-", i_val), undirected_pair = undirected_pair_val, is_border = TRUE,     is_border_multiple = touches_other[touches_to], other_adjacent_borders = other_clusters[touches_to],     mode = ""), spot_id, x, y, directed_pair, undirected_pair,     cluster, interface, is_border, is_border_multiple, other_adjacent_borders,     mode): Can't select columns that don't exist.
#> ✖ Column `spot_id` doesn't exist.
head(all_borders)
#> Error: object 'all_borders' not found
```
