# Build a single line of spots along a segment

Selects spots near the segment defined by endpoints \`A\` and \`B\`
using \`build_one_trajectory()\`, then orders the selected spots by
their projection parameter (from \`A\` to \`B\`). The result is returned
as a data frame.

## Usage

``` r
build_one_line(df_rest, A, B, top_n = 19, max_dist = NULL)
```

## Arguments

- df_rest:

  A data frame of candidate spots containing at least columns \`x\` and
  \`y\`.

- A:

  A data frame with columns \`x\` and \`y\` defining endpoint A of the
  segment. If it contains multiple rows, only the first row is used.

- B:

  A data frame with columns \`x\` and \`y\` defining endpoint B of the
  segment. If it contains multiple rows, only the first row is used.

- top_n:

  Integer; number of closest spots to keep (default \`19\`). If
  \`NULL\`, no top-N truncation is applied.

- max_dist:

  Numeric; if not \`NULL\`, only spots with distance to the segment
  \`\<= max_dist\` are kept.

## Value

A data frame of selected spots (subset of \`df_rest\`) ordered along the
segment (in increasing \`pos_on_seg\`). The output includes the extra
columns added by \`build_one_trajectory()\` (typically \`dist_to_seg\`
and \`pos_on_seg\`).

## Details

This function is a thin wrapper around \`build_one_trajectory()\` that
returns only the ordered \`selected\` data frame (and not the segment
metadata).

It uses the base R pipe \`\|\>\` and \`dplyr::arrange()\`.

## Examples

``` r
df <- data.frame(
  x = c(0, 1, 2, 3, 4, 2),
  y = c(0, 0, 0, 0, 0, 1),
  id = 1:6
)
A <- data.frame(x = 0, y = 0)
B <- data.frame(x = 4, y = 0)

build_one_line(df, A, B, top_n = 5)
#>   x y id dist_to_seg pos_on_seg trajectory_id
#> 1 0 0  1           0       0.00          main
#> 2 1 0  2           0       0.25          main
#> 3 2 0  3           0       0.50          main
#> 4 3 0  4           0       0.75          main
#> 5 4 0  5           0       1.00          main
```
