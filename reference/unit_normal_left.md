# Compute the left unit normal of the directed segment A-\>B

Returns the unit-length normal vector pointing to the left of the
directed segment \\A \rightarrow B\\. For a direction vector \\v = (vx,
vy)\\, the left normal is \\(-vy, vx)\\ normalized by \\\|\|v\|\|\\.

## Usage

``` r
unit_normal_left(A, B)
```

## Arguments

- A:

  A data frame with columns \`x\` and \`y\` representing point A. If it
  contains multiple rows, only the first row is used.

- B:

  A data frame with columns \`x\` and \`y\` representing point B. If it
  contains multiple rows, only the first row is used.

## Value

A numeric vector of length 2 with named components:

- nx:

  x-component of the left unit normal.

- ny:

  y-component of the left unit normal.

## Examples

``` r
# Horizontal segment to the right: (0,0) -> (1,0)
#    Left normal points upward: (0, 1)
A <- data.frame(x = 0, y = 0)
B <- data.frame(x = 1, y = 0)
unit_normal_left(A, B)
#> nx ny 
#>  0  1 
```
