# Shift a point along a given direction vector

Translates a point \`P\` by an amount \`offset\` in the direction \`(nx,
ny)\`. The returned point is: \$\$P' = (P_x + offset \cdot nx,\\ P_y +
offset \cdot ny).\$\$

## Usage

``` r
shift_point(P, nx, ny, offset)
```

## Arguments

- P:

  A data frame with columns \`x\` and \`y\` representing the point to
  shift. If it contains multiple rows, only the first row is used.

- nx:

  Numeric scalar; x-component of the direction vector.

- ny:

  Numeric scalar; y-component of the direction vector.

- offset:

  Numeric scalar; translation magnitude (in the same units as
  \`x\`/\`y\`). Positive values move in the \`(nx, ny)\` direction;
  negative values move in the opposite direction.

## Value

A one-row data frame with columns \`x\` and \`y\` giving the shifted
point.

## Examples

``` r
# Shift the point (1, 2) by 3 units in the direction (0, 1)
P <- data.frame(x = 1, y = 2)
shift_point(P, nx = 0, ny = 1, offset = 3)
#>   x y
#> 1 1 5
```
