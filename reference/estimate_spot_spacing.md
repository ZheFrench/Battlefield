# Estimate spot spacing from nearest-neighbor distances

Estimates the typical spacing between spots by computing each spot's
nearest-neighbor distance (excluding itself) and returning the median of
those distances. For large datasets, the estimate is computed on a
random subsample for speed.

## Usage

``` r
estimate_spot_spacing(df, sample_n = 1000)
```

## Arguments

- df:

  A data frame containing at least columns \`x\` and \`y\` (numeric)
  representing spot coordinates.

- sample_n:

  Integer; maximum number of spots to sample (default \`1000\`). If
  \`nrow(df) \> sample_n\`, a random subset of size \`sample_n\` is
  used; otherwise all spots are used.

## Value

A numeric scalar: the median nearest-neighbor distance (in the same
units as \`x\` and \`y\`).

## Details

Nearest neighbors are computed with
[`RANN::nn2()`](https://jefferislab.github.io/RANN/reference/nn2.html)
using \`k = 2\`. The first neighbor is the point itself (distance 0), so
the function uses the second neighbor distance \`nn\$nn.dists\[, 2\]\`
as the true nearest-neighbor distance.

Missing values are handled via \`median(..., na.rm = TRUE)\`.

## Examples

``` r
set.seed(1)
df <- data.frame(
  x = rep(1:5, each = 5),
  y = rep(1:5, times = 5)
)
estimate_spot_spacing(df)
#> [1] 1
```
