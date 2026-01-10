# Select spots near a segment and order them along the segment

Computes the distance from each spot in \`df\` to the segment \`\[p0,
p1\]\`, optionally filters by a maximum distance, keeps the \`top_n\`
closest spots, and finally orders the retained spots by their projection
position \`pos_on_seg\` along the segment (from \`p0\` to \`p1\`). This
is the foundation function for building trajectory lines across spatial
spots.

## Usage

``` r
build_one_trajectory(df, p0, p1, top_n = 100, max_dist = NULL)
```

## Arguments

- df:

  A data frame of spots containing at least columns \`x\` and \`y\`.
  Additional columns are preserved in the output.

- p0:

  A data frame with columns \`x\` and \`y\` defining the first endpoint
  of the segment. If it has multiple rows, only the first row is used.

- p1:

  A data frame with columns \`x\` and \`y\` defining the second endpoint
  of the segment. If it has multiple rows, only the first row is used.

- top_n:

  Integer; number of closest spots to keep (default \`100\`). If
  \`NULL\`, no top-N truncation is applied.

- max_dist:

  Numeric; if not \`NULL\`, only spots with distance to the segment
  \`\<= max_dist\` are kept.

## Value

A data frame containing the selected spots, with two extra columns:

- dist_to_seg:

  Euclidean distance from the spot to the segment.

- pos_on_seg:

  Clamped projection parameter in \`\[0, 1\]\` indicating position along
  the segment (0 at \`p0\`, 1 at \`p1\`).

The rows are ordered by \`pos_on_seg\` (i.e., from \`p0\` to \`p1\`).

## Details

Distances are computed using
[`point_segment_distance_vec()`](https://zhefrench.github.io/Battlefield/reference/point_segment_distance_vec.md).
Ordering uses the projection parameter \\pos_on_seg =
((P-A)\cdot(B-A))/\|\|B-A\|\|^2\\ clamped to \`\[0, 1\]\`.

This function uses `dplyr` verbs (\`mutate\`, \`filter\`, \`arrange\`,
\`slice_head\`) via the \`\|\>\` pipe.

## Examples

``` r
# Minimal example with base data.frame inputs
df <- data.frame(
  x = c(0, 1, 2, 3, 4),
  y = c(0, 1, 1, 2, 4),
  id = letters[1:5]
)
p0 <- data.frame(x = 0, y = 0)
p1 <- data.frame(x = 4, y = 0)

# Keep 3 closest spots (no distance threshold)
res <- build_one_trajectory(df, p0, p1, top_n = 3)
head(res)
#>   x y id dist_to_seg pos_on_seg trajectory_id
#> 1 0 0  a           0       0.00          main
#> 2 1 1  b           1       0.25          main
#> 3 2 1  c           1       0.50          main
```
