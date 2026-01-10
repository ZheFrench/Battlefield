# Rasterize a line between two points using Bresenham's algorithm

Generates integer grid coordinates along the line segment connecting two
points using Bresenham's line algorithm.

## Usage

``` r
bresenham_line(p0, p1, snap = TRUE)
```

## Arguments

- p0:

  A \`data.frame\` with at least columns \`x\` and \`y\`. If it contains
  multiple rows, only the first row is used.

- p1:

  A \`data.frame\` with at least columns \`x\` and \`y\`. If it contains
  multiple rows, only the first row is used.

- snap:

  Logical; if \`TRUE\` (default), \`p0\` and \`p1\` coordinates are
  rounded to the nearest integer before running the algorithm. If
  \`FALSE\`, coordinates are truncated via \`as.integer()\`.

## Value

A \`data.frame\` with columns \`x\` and \`y\` giving the integer grid
points visited by the line, including both endpoints, in traversal order
from \`p0\` to \`p1\`.

## Details

The function:

- Validates inputs are data frames with \`x\`/\`y\` columns and at least
  one row.

- Converts endpoints to integer coordinates (optionally rounding first).

- Applies the classic Bresenham algorithm to enumerate grid points.

This is useful for tracing discrete paths on an integer lattice (e.g.,
pixels, tile grids, spatial transcriptomics spot grids).

## Examples

``` r
p0 <- data.frame(x = 1.2, y = 2.7)
p1 <- data.frame(x = 7.9, y = 5.1)

# Default (snap = TRUE): rounds endpoints first
bresenham_line(p0, p1)
#>   x y
#> 1 1 3
#> 2 2 3
#> 3 3 4
#> 4 4 4
#> 5 5 4
#> 6 6 4
#> 7 7 5
#> 8 8 5

# No rounding (snap = FALSE): integer coercion truncates toward zero
bresenham_line(p0, p1, snap = FALSE)
#>   x y
#> 1 1 2
#> 2 2 2
#> 3 3 3
#> 4 4 3
#> 5 5 4
#> 6 6 4
#> 7 7 5
```
