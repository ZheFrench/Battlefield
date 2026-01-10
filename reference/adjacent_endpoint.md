# Pick a point adjacent to a selected endpoint, on a given side of a segment

Given an endpoint (typically the start or end of a previously selected
path), this function computes a target point located one perpendicular
step away from the endpoint, on the "left" or "right" side relative to
the directed segment \\A \rightarrow B\\. It then returns the closest
spot to that target among the remaining candidates \`df_rest\`.

## Usage

``` r
adjacent_endpoint(df_rest, endpoint, A, B, spacing, side = c("left", "right"))
```

## Arguments

- df_rest:

  A data frame of candidate spots containing at least columns \`x\` and
  \`y\`.

- endpoint:

  A one-row data frame with columns \`x\` and \`y\` representing the
  endpoint from which to step sideways. If it has multiple rows, only
  the first row is used.

- A:

  A data frame with columns \`x\` and \`y\` representing point A
  defining the direction of the segment. Only the first row is used.

- B:

  A data frame with columns \`x\` and \`y\` representing point B
  defining the direction of the segment. Only the first row is used.

- spacing:

  Numeric scalar; step size (in the same coordinate units as
  \`x\`/\`y\`) used to move perpendicularly from \`endpoint\`.

- side:

  Character; which side to step to relative to the vector \\A
  \rightarrow B\\. Must be \`"left"\` or \`"right"\`.

## Value

A one-row data frame (same columns as \`df_rest\`) corresponding to the
spot closest to the computed perpendicular target point.

## Details

The unit direction vector is computed from \\A \rightarrow B\\ and a
unit left normal \\(-vy, vx)/\|\|v\|\|\\ is derived. The target point
is: \$\$T = endpoint + sign \cdot spacing \cdot n\$\$ where \`sign =
+1\` for \`"left"\` and \`-1\` for \`"right"\`.

The closest spot is selected with
[`closest_spot()`](https://zhefrench.github.io/Battlefield/reference/closest_spot.md)
using squared Euclidean distance.

## Examples

``` r
df <- data.frame(x = c(0, 1, 1, 2, 2), y = c(0, 0, 1, 0, 1), id = 1:5)
A <- data.frame(x = 0, y = 0)
B <- data.frame(x = 2, y = 0)
endpoint <- data.frame(x = 1, y = 0)

# Step "left" of A->B (here, toward positive y)
adjacent_endpoint(df, endpoint, A, B, spacing = 1, side = "left")
#>   x y id
#> 3 1 1  3

# Step "right" of A->B (here, toward negative y)
adjacent_endpoint(df, endpoint, A, B, spacing = 1, side = "right")
#>   x y id
#> 2 1 0  2
```
