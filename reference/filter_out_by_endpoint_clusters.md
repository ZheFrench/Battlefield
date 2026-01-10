# Filter lines by endpoint cluster membership

Filters trajectory data by checking the cluster labels at each
trajectory's endpoints. For every \`trajectory_id\`, the start endpoint
is defined as the spot with the smallest projection parameter
\`pos_on_seg\`, and the end endpoint as the spot with the largest
\`pos_on_seg\`. Only trajectories whose start cluster is in
\`allowed_start_clusters\` \*and\* whose end cluster is in
\`allowed_end_clusters\` are kept.

## Usage

``` r
filter_out_by_endpoint_clusters(
  out,
  allowed_start_clusters,
  allowed_end_clusters
)
```

## Arguments

- out:

  A data frame containing selected spots with columns \`trajectory_id\`,
  \`cluster\`, and \`pos_on_seg\`, typically the output of
  \`build_similar_trajectories()\`.

- allowed_start_clusters:

  Vector of allowed cluster labels for the start endpoint.

- allowed_end_clusters:

  Vector of allowed cluster labels for the end endpoint.

## Value

The same data frame as \`out\`, filtered to keep only the lines matching
the allowed endpoint cluster constraints.

## Details

The endpoint clusters are computed per \`trajectory_id\`:

- \`start_cluster = cluster\[which.min(pos_on_seg)\]\`

- \`end_cluster = cluster\[which.max(pos_on_seg)\]\`

If multiple spots share the same minimum/maximum \`pos_on_seg\`, the
first is taken (as per \`which.min()\` / \`which.max()\`).

This function uses \`dplyr\` (\`group_by\`, \`summarise\`, \`filter\`,
\`pull\`) and the base R pipe \`\|\>\`.

## Examples

``` r
# Minimal example
out <- data.frame(
  trajectory_id = c("L1","L1","L2","L2"),
  pos_on_seg = c(0.0, 1.0, 0.0, 1.0),
  cluster  = c("A", "B", "A", "C"),
  x = 1:4, y = 1:4
)

# Keep only trajectories starting in A and ending in B
filter_out_by_endpoint_clusters(out, allowed_start_clusters = "A",
allowed_end_clusters = "B")
#>   trajectory_id pos_on_seg cluster x y
#> 1            L1          0       A 1 1
#> 2            L1          1       B 2 2
```
