# Detect the grid type (square vs hexagonal) from spatial coordinates

This function inspects k-nearest-neighbor distances in a 2D coordinate
set to infer whether the points lie on a \*\*square\*\* or
\*\*hexagonal\*\* lattice.

## Usage

``` r
detect_grid_type(
  df,
  coords = c("x", "y"),
  k = 12,
  tolerance = 0.1,
  verbose = TRUE
)
```

## Arguments

- df:

  A data.frame containing spatial coordinates.

- coords:

  Character vector of length 2 giving the coordinate column names.
  Default is \`c("x","y")\`.

- k:

  Integer. Number of neighbors to consider (internally uses \`k + 1\` to
  include the point itself, then removes it). Default is 12.

- tolerance:

  Numeric. Tolerance for matching the characteristic ratio to
  \\\sqrt{2}\\ or \\\sqrt{3}\\. Default is 0.1.

- verbose:

  Logical. If \`TRUE\`, prints diagnostic messages. Default is \`TRUE\`.

## Value

A list with:

- grid_type:

  One of \`"square"\`, \`"hexagonal"\`, or \`"unknown"\`.

- step:

  Estimated grid step size (median nearest-neighbor distance).

- median_ratio:

  Median normalized distance ratio in the (1.2, 2.0) range.

## Details

The method:

1.  Computes the nearest-neighbor distance for each point and uses the
    median of those distances as an estimate of the grid step size.

2.  Normalizes all neighbor distances by this step size.

3.  Looks for a characteristic secondary distance peak:

    - square grid: ratio \\\approx \sqrt{2}\\

    - hexagonal grid: ratio \\\approx \sqrt{3}\\

The grid step size is estimated as the median of the per-point
nearest-neighbor distance. The "median distance ratio" is computed from
normalized neighbor distances restricted to the interval (1.2, 2.0),
which tends to capture the second-shell distances on regular grids.

## Examples

``` r
# --- Example 1: square grid ---
sq <- expand.grid(x = 0:9, y = 0:9)
res_sq <- detect_grid_type(sq, coords = c("x","y"), k = 12, tolerance = 0.1)
#> ---- Grid type detection ----
#> Estimated grid step: 1
#> Median distance ratio: 1.414
#> Detected grid type: square
#> --------------------------------
res_sq$grid_type
#> [1] "square"

# --- Example 2: hexagonal grid (pointy-top axial-like layout) ---
nx <- 10; ny <- 10
hex <- expand.grid(i = 0:(nx-1), j = 0:(ny-1))
hex$x <- hex$i + 0.5 * (hex$j %% 2)
hex$y <- (sqrt(3)/2) * hex$j
hex <- hex[, c("x","y")]
res_hex <- detect_grid_type(hex, coords = c("x","y"), k = 12, tolerance =
0.1)
#> ---- Grid type detection ----
#> Estimated grid step: 1
#> Median distance ratio: 1.732
#> Detected grid type: hexagonal
#> --------------------------------
res_hex$grid_type
#> [1] "hexagonal"
```
