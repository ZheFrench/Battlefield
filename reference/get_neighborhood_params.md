# Suggest neighborhood parameters based on detected grid type

This helper detects the spatial grid type (square vs hexagonal) and
returns recommended neighborhood parameters for building local adjacency
/ neighbor queries (e.g., for border detection or spatial graphs).

## Usage

``` r
get_neighborhood_params(
  df,
  coords = c("x", "y"),
  square_connectivity = c(4, 8),
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

- square_connectivity:

  Integer-like choice for square grids: \`4\` (Von Neumann) or \`8\`
  (Moore). Default is \`c(4, 8)\` which selects the first value (\`4\`).

- tolerance:

  Numeric. Passed to \[detect_grid_type()\] to match characteristic
  ratios (sqrt(2) / sqrt(3)). Default is 0.1.

- verbose:

  Logical. If \`TRUE\`, prints diagnostic messages. Default is \`TRUE\`.

## Value

A named list with neighborhood parameters: \`grid_type\`,
\`connectivity\`, \`radius\`, \`k\`, and \`comment\`.

## Details

It relies on \[detect_grid_type()\] to estimate the grid step size and
grid geometry, then chooses typical parameters:

- \*\*Hexagonal (Visium)\*\*: 6-neighborhood, radius ~ step

- \*\*Square (Visium HD)\*\*: 4- or 8-neighborhood, radius ~ step (4) or
  ~ sqrt(2)\*step (8)

Returned values:

- grid_type:

  Detected grid type: \`"hexagonal"\`, \`"square"\`, or \`"unknown"\`.

- connectivity:

  Nominal grid connectivity (6 for hex, 4 or 8 for square).

- r:

  Recommended distance threshold (radius) to define adjacency.

- k_max:

  A suggested upper bound for kNN queries (used as a safe cap).

- comment:

  Human-readable description of the choice.

Notes:

- \`radius\` is set to \`1.01 \* step\` (or \`1.01 \* sqrt(2) \* step\`
  for \#' 8-connectivity) to be slightly permissive under small
  numerical noise.

- \`k\` is deliberately larger than the nominal connectivity to make
  sure enough candidates are retrieved before applying distance-based
  filtering.

## Examples

``` r
# Square grid example
sq <- expand.grid(x = 0:9, y = 0:9)
get_neighborhood_params(sq, square_connectivity = 4, verbose = FALSE)

# Hexagonal grid example
nx <- 10; ny <- 10
hex <- expand.grid(i = 0:(nx-1), j = 0:(ny-1))
hex$x <- hex$i + 0.5 * (hex$j %% 2)
hex$y <- (sqrt(3)/2) * hex$j
hex <- hex[, c("x","y")]
get_neighborhood_params(hex, verbose = FALSE)
```
