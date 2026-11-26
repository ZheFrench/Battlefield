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
head(all_borders)
#>            x           y directed_pair undirected_pair cluster interface
#> 1  0.1836433  1.68887329           A-C       A-C / C-A       A         C
#> 2 -0.8356286  1.58658843           A-C       A-C / C-A       A         C
#> 3  0.3295078 -2.28523554           A-C       A-C / C-A       A         C
#> 4  0.4874291  0.66706617           A-C       A-C / C-A       A         C
#> 5  0.5757814 -0.01339952           A-C       A-C / C-A       A         C
#> 6 -0.3053884  0.51010842           A-C       A-C / C-A       A         C
#>   is_border is_border_multiple other_adjacent_borders  mode
#> 1      TRUE               TRUE                      B inner
#> 2      TRUE              FALSE                   <NA> inner
#> 3      TRUE               TRUE                      B inner
#> 4      TRUE              FALSE                   <NA> inner
#> 5      TRUE               TRUE                      B inner
#> 6      TRUE               TRUE                      B inner
```
