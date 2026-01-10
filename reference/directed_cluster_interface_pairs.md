# Build all oriented cluster pairs (A -\> B) and their directed_pair labels

Given a vector of cluster labels, this function returns all
\*\*ordered\*\* (oriented) pairs of distinct clusters. For each pair
(\`cluster\`, \`interface\`), it also creates a \`directed_pair\` string
label such as \`"A-B"\`.

## Usage

``` r
directed_cluster_interface_pairs(
  cluster_labels,
  interface_separator = "-",
  sort_clusters = FALSE
)
```

## Arguments

- cluster_labels:

  A vector of cluster labels (character, factor, numeric, etc.). \`NA\`
  values are ignored.

- interface_separator:

  Character string used to join \`cluster\` and \`interface\` into the
  \`directed_pair\` label. Default is \`"-"\`.

- sort_clusters:

  Logical. If \`TRUE\`, unique cluster labels are sorted
  (lexicographically) before building pairs. If \`FALSE\`, preserves
  first appearance order. Default is \`FALSE\`.

## Value

A data.frame with columns:

- cluster:

  Source cluster label (character).

- interface:

  Target cluster label (character).

- directed_pair:

  Directed pair label, e.g. "A-B".

If fewer than 2 unique (non-NA) clusters are present, returns an empty
data.frame with the same columns.

## Examples

``` r
cl <- c("1","1","2","3", NA)
directed_cluster_interface_pairs(cl)
#>   cluster interface directed_pair
#> 1       2         1           2-1
#> 2       3         1           3-1
#> 3       1         2           1-2
#> 4       3         2           3-2
#> 5       1         3           1-3
#> 6       2         3           2-3
```
