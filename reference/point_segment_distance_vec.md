# Point-to-segment distance (vectorized)

Computes the Euclidean distance from one or more points \`(px, py)\` to
the line segment defined by endpoints \`A(ax, ay)\` and \`B(bx, by)\`.
The projection is clamped to the segment (i.e., \`pos_on_seg\` is
restricted to \`\[0, 1\]\`).

## Usage

``` r
point_segment_distance_vec(px, py, ax, ay, bx, by)
```

## Arguments

- px, py:

  Numeric vectors (or scalars) giving the x/y coordinates of the query
  point(s).

- ax, ay:

  Numeric scalars giving the x/y coordinates of endpoint A.

- bx, by:

  Numeric scalars giving the x/y coordinates of endpoint B.

## Value

A numeric vector of distances, with length equal to \`max(length(px),
length(py))\` (after R's usual vector recycling rules).

## Details

Let \`v = B - A\` and \`pos_on_seg = ((P - A) · v) / \|\|v\|\|^2\`. The
closest point on the infinite line is \`A + pos_on_seg v\`; clamping
\`pos_on_seg\` to \`\[0, 1\]\` yields the closest point on the segment.

If \`A\` and \`B\` are identical (\`\|\|v\|\|^2 == 0\`), the function
returns the distance from \`P\` to \`A\`.

## Examples

``` r
# Single point to a horizontal segment
point_segment_distance_vec(px = 1, py = 2, ax = 0, ay = 0, bx = 3, by = 0)
#> [1] 2

# Multiple points (vectorized)
px <- c(0, 1, 2, 3)
py <- c(1, 1, 1, 1)
point_segment_distance_vec(px, py, ax = 0, ay = 0, bx = 3, by = 0)
#> [1] 1 1 1 1

# Degenerate segment (A == B)
point_segment_distance_vec(px = c(0, 1), py = c(0, 1), ax = 0, ay = 0, bx =
0, by = 0)
#> [1] 0.000000 1.414214
```
