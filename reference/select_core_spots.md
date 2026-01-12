# Select core (non-interface) spots for a directed pair

This function selects core (non-interface) spots from a cluster that
match the count of border spots for a specific directed pair. The
\`mode\` value is read from the \`border_df\` dataframe for this
cluster/interface pair. If multiple mode values exist in border_df for
the pair, a warning is printed and the first value is used.

## Usage

``` r
select_core_spots(
  df,
  border_df,
  cluster,
  interface,
  mode = "both",
  coord_cols = c("x", "y"),
  cluster_col = "cluster"
)
```

## Arguments

- df:

  A data.frame containing at least the coordinate columns and a cluster
  label column.

- border_df:

  A data.frame of border spots from \[build_all_borders()\], containing
  columns \`spot_id\`, \`cluster\`, \`interface\`, \`directed_pair\`,
  etc.

- cluster:

  Character. Label of the cluster from which to select inner spots.

- interface:

  Character. Label of the target cluster for the directed pair.

- mode:

  Character. One of "inner", "outer", or "both". Controls which mode
  values from border_df to use when counting border spots. When
  mode="both", both "inner" and "outer" mode rows from border_df are
  used. Default is "both".

- coord_cols:

  Character vector of length 2 giving the coordinate column names.
  Default is \`c("x","y")\`.

- cluster_col:

  Character. Name of the column containing cluster labels. Default is
  \`"cluster"\`.

## Value

A data.frame of sampled core spots from \`cluster\`, with all columns
from \`df\` plus \`interface\` annotation.

## Details

For example with \`cluster="1"\` and \`interface="2"\`: - If mode in df
is "inner": counts only border spots from 1→2 - If mode in df is
"outer": counts only border spots from 2→1 - If mode in df is "both":
counts border spots from both 1→2 AND 2→1

This returns N random inner spots from the cluster, where N equals the
border count.

Steps:

1.  Counts directed border spots based on \`mode\` parameter.

2.  Gets all spots in \`cluster\` that are not at any border.

3.  Randomly samples the same number of core spots as the border count.

4.  Returns these sampled core spots with \`interface\` annotation.

If not enough core spots exist to match the border count, all available
core spots are returned with a warning.

## Examples

``` r
# Create synthetic spatial transcriptomics data with 3 clusters
df_ex <- data.frame(
  spot_id = paste0("spot_", seq_len(300)),
  x = rnorm(300),
  y = rnorm(300),
  cluster = sample(c("A", "B", "C"), 300, replace = TRUE)
)

# Step 1: Identify all border spots between clusters
all_borders <- build_all_borders(df_ex, k = 6)

# Step 2: Select core spots for the directed pair A→B
# matching the number of border spots found
cores_inner <- select_core_spots(
  df_ex,
  border_df = all_borders,
  cluster = "A",
  interface = "B",
  mode = "inner"
)
#> Not enough core spots for pair(cluster=A, interface=B, mode=inner). Requested:85, Available:0. Returning all available.
#> Error in `$<-.data.frame`(`*tmp*`, "interface", value = "B"): replacement has 1 row, data has 0
head(cores_inner)
#> Error: object 'cores_inner' not found

# Step 3: Select core spots considering both directions (A→B AND B→A)
cores_both <- select_core_spots(
  df_ex,
  border_df = all_borders,
  cluster = "A",
  interface = "B",
  mode = "both"
)
#> Not enough core spots for pair(cluster=A, interface=B, mode=both). Requested:170, Available:0. Returning all available.
#> Error in `$<-.data.frame`(`*tmp*`, "interface", value = "B"): replacement has 1 row, data has 0
head(cores_both)
#> Error: object 'cores_both' not found
```
