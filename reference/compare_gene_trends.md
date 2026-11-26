# Compare gene trends along a fixed trajectory

Tests whether genes decrease or increase along trajectories already
stored by \[add_trajectories_to_spe()\], and reports how similar their
smoothed profile is to that of a reference gene. Spots are fixed; only
the genes change.

## Usage

``` r
compare_gene_trends(
  spe,
  genes,
  reference,
  assay_name = "logcounts",
  feature_names = NULL,
  trajectory_ids = NULL,
  min_detected = 0.1,
  span = 0.8,
  fdr_cutoff = 0.05
)
```

## Arguments

- spe:

  A SpatialExperiment with \`trajectory_id\` and \`pos_on_seg\` in its
  colData.

- genes:

  Character vector of genes to test.

- reference:

  Reference gene. It is reported but never tested.

- assay_name:

  Assay holding normalized expression.

- feature_names:

  Optional rowData column with gene names. Row names are used if
  \`NULL\`.

- trajectory_ids:

  Trajectories to pool. All are used if \`NULL\`.

- min_detected:

  Minimal fraction of trajectory spots expressing a gene.

- span:

  Fraction of spots used by the local smoother, in (0, 1\].

- fdr_cutoff:

  FDR threshold used for \`trend\`.

## Value

A data.frame with one row per gene: the reference first, then the tested
genes by increasing slope, then the genes that could not be tested
(statistics are \`NA\`), with columns \`gene\`, \`slope\`, \`r_curve\`,
\`fdr_down\`, \`fdr_up\`, \`trend\` and \`status\`. \`status\` is
"reference", "tested", "not found", "duplicated", "missing values" or
"too sparse or constant". Smoothed profiles, for plotting, are in
\`attr(result, "curves")\`.

## Details

\`slope\` is the estimated change in expression from the start to the
end of the trajectory (\`pos_on_seg\` from 0 to 1). It comes from a
generalized least squares model adjusted for the trajectory line.
Residuals of neighbouring spots are allowed to be correlated
(exponential covariance on the spatial coordinates, plus a nugget), with
parameters chosen per gene by REML. P-values are one-sided and adjusted
by Benjamini-Hochberg over tested genes. Smoothed profiles pool the
spots of all lines. \`r_curve\` is the Pearson correlation of each
profile with that of the reference; it is descriptive and has no
p-value. A spatial pattern about as long as the trajectory cannot be
told apart from a trend, so the FDR remain optimistic in that case: use
them to rank genes within a sample, and replicate across samples to
conclude. \`trend\` is "dec", "inc", "unknown" or "ref" (\`NA\` if not
tested). Do not use the gene that guided the choice of trajectory as
evidence.

## Examples

``` r
data("visiumHD_16um_simulated_spe", package = "Battlefield")
spe <- visiumHD_16um_simulated_spe
df <- data.frame(
  spot_id = colnames(spe),
  x = SpatialExperiment::spatialCoords(spe)[, 1],
  y = SpatialExperiment::spatialCoords(spe)[, 2],
  cluster = SummarizedExperiment::colData(spe)$cluster
)
centroids <- compute_centroids(df)
A <- centroids[centroids$cluster == 3, c("x", "y")]
B <- centroids[centroids$cluster == 4, c("x", "y")]
trajs <- build_similar_trajectories(df, A, B, top_n = 12, n_extra = 2,
  side = "both")

# The example data hold a single gene: simulate genes along the trajectory
set.seed(1)
pos <- rep(0.5, nrow(df))
pos[match(trajs$spot_id, df$spot_id)] <- trajs$pos_on_seg
sim <- function(slope) 3 + slope * pos + stats::rnorm(length(pos), sd = 0.3)
logcounts <- rbind(REF = sim(-2), DOWN = sim(-1.5), UP = sim(1.5),
  FLAT = sim(0))
colnames(logcounts) <- df$spot_id
toy <- SpatialExperiment::SpatialExperiment(
  assays = list(logcounts = logcounts),
  spatialCoords = as.matrix(df[, c("x", "y")])
)
toy <- add_trajectories_to_spe(toy, trajectory = trajs)

compare_gene_trends(toy, genes = c("DOWN", "UP", "FLAT"),
  reference = "REF")
#>   gene      slope    r_curve     fdr_down       fdr_up   trend    status
#> 1  REF -1.8237565  1.0000000           NA           NA     ref reference
#> 2 DOWN -1.6949907  0.9960631 3.499210e-13 1.000000e+00     dec    tested
#> 3 FLAT -0.1801978  0.8093483 1.212369e-01 1.000000e+00 unknown    tested
#> 4   UP  1.6096922 -0.9953335 1.000000e+00 1.261951e-15     inc    tested
```
