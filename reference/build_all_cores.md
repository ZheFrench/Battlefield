# Build core spots for all oriented cluster pairs

This function iterates over all \*\*ordered\*\* cluster pairs (A -\> B)
and returns a single data.frame containing the core (control) spots for
each pair, as computed by \[select_core_spots()\].

## Usage

``` r
build_all_cores(
  df,
  border_df,
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

- border_df:

  A data.frame of border spots from \[build_all_borders()\].

- mode:

  Character. One of "inner", "outer", or "both". Controls which mode
  values from border_df to use when counting border spots. Default is
  "both".

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

A data.frame produced by row-binding the result of core spot selection
for each oriented pair. Typically contains the original columns of
\`df\` plus annotation columns (is_core, interface, mode) from the core
spot selector.

\*\*Note on mode column\*\*: The \`mode\` column in the returned
data.frame reflects which border modes were used for each core spot
selection: "inner", "outer", or "both". This differs from
\`build_all_borders()\` which only returns "inner" or "outer".

## Details

If \`pairs\` is not provided, it is generated from the cluster labels
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
all_cores <- build_all_cores(df_ex, all_borders, mode = "both")
#> Error: object 'all_borders' not found
head(all_cores)
#> Error: object 'all_cores' not found
```
