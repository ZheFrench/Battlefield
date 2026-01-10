# Add trajectory information to SpatialExperiment colData

This function takes trajectory results (e.g., from
\[build_similar_trajectories()\]) and adds trajectory metadata to the
colData of a SpatialExperiment object.

## Usage

``` r
add_trajectories_to_spe(spe, trajectory = NULL, erase = FALSE)
```

## Arguments

- spe:

  A SpatialExperiment object from which the trajectory dataframes were
  derived.

- trajectory:

  A data.frame of trajectory spots. Expected output from
  \[build_similar_trajectories()\]. Expected columns: spot_id,
  trajectory_id, offset, pos_on_seg, dist_to_seg.

- erase:

  Logical. If \`TRUE\`, erase pre-existing trajectory columns before
  adding new ones. If \`FALSE\` (default), issue a warning if columns
  already exist and skip adding them.

## Value

The SpatialExperiment object with updated colData.

## Details

The \`trajectories\` parameter must be provided and contain at least one
row.

The following columns are added: - \`trajectory_id\`: character,
identifier for each trajectory (e.g., "main", "left_1", "right_2") -
\`offset\`: numeric, the offset distance used to build each trajectory
line - \`pos_on_seg\`: numeric in \[0, 1\], the position along the
trajectory - \`dist_to_seg\`: numeric, the distance to the original
segment

## Examples

``` r
data("visiumHD_16um_simulated_spe", package = "Battlefield")
spe <- visiumHD_16um_simulated_spe
df <- data.frame(
  spot_id = colnames(spe),
  x = spatialCoords(spe)[, 1],
  y = spatialCoords(spe)[, 2],
  cluster = colData(spe)$cluster
)
#> Error in spatialCoords(spe): could not find function "spatialCoords"
centroids <- compute_centroids(df)
#> Error in model.frame.default(formula = cbind(x, y) ~ cluster, data = df): 'data' must be a data.frame, environment, or list
A <- centroids[centroids$cluster == 1, c("x", "y")]
#> Error: object 'centroids' not found
B <- centroids[centroids$cluster == 3, c("x", "y")]
#> Error: object 'centroids' not found
trajs <- build_similar_trajectories(df, A, B, top_n = 10, n_extra = 1, side =
"both")
#> Error in estimate_spot_spacing(df): all(c("x", "y") %in% names(df)) is not TRUE
spe <- add_trajectories_to_spe(spe, trajectory = trajs)
#> Error: object 'trajs' not found
```
