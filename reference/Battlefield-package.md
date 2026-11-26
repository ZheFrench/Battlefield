# Battlefield: Swiss-army toolkit for spatial transcriptomics region selection

Battlefield provides low-level, modular utilities to define and extract
spatial spots from specific regions in spatial transcriptomics data. It
supports four core workflows: identifying cluster interfaces,
characterizing intra-cluster layers, modeling inter-cluster
trajectories, and analyzing neighborhood composition.

## Key Features

- **Cluster Interfaces**:

  Identify and analyze boundary regions where distinct tissue or cell
  populations interact. See
  [`select_border_spots`](https://zhefrench.github.io/Battlefield/reference/select_border_spots.md),
  [`build_all_borders`](https://zhefrench.github.io/Battlefield/reference/build_all_borders.md),
  [`directed_cluster_interface_pairs`](https://zhefrench.github.io/Battlefield/reference/directed_cluster_interface_pairs.md).

- **Intra-cluster Layers**:

  Characterize spatial layers or gradients within a single cluster
  (border, intermediate, core). See
  [`create_cluster_layers`](https://zhefrench.github.io/Battlefield/reference/create_cluster_layers.md),
  [`create_all_layers`](https://zhefrench.github.io/Battlefield/reference/create_all_layers.md).

- **Inter-cluster Trajectories**:

  Model spatial transitions across multiple clusters via
  centroid-to-centroid Bresenham paths. See
  [`build_one_trajectory`](https://zhefrench.github.io/Battlefield/reference/build_one_trajectory.md),
  [`build_similar_trajectories`](https://zhefrench.github.io/Battlefield/reference/build_similar_trajectories.md).

- **Neighborhood Composition**:

  Analyze k-NN cluster composition with optional distance constraints.
  See
  [`count_neighborhood`](https://zhefrench.github.io/Battlefield/reference/count_neighborhood.md),
  [`get_neighborhood_spots`](https://zhefrench.github.io/Battlefield/reference/get_neighborhood_spots.md).

## Core Workflow

1\. Start with a
[`SpatialExperiment`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html)
object containing spot coordinates and cluster annotations. 2. Extract
spot metadata via `colData(spe)` to create a working `data.frame`. 3.
Apply selection functions (e.g.,
[`select_border_spots`](https://zhefrench.github.io/Battlefield/reference/select_border_spots.md),
[`create_cluster_layers`](https://zhefrench.github.io/Battlefield/reference/create_cluster_layers.md))
to identify regions of interest. 4. Integrate results back to `spe`
using helper functions:
[`add_borders_to_spe`](https://zhefrench.github.io/Battlefield/reference/add_borders_to_spe.md),
[`add_layers_to_spe`](https://zhefrench.github.io/Battlefield/reference/add_layers_to_spe.md),
[`add_trajectories_to_spe`](https://zhefrench.github.io/Battlefield/reference/add_trajectories_to_spe.md).
5. Use annotated `spe` for downstream analyses (e.g., ligand–receptor
interactions via BulkSignalR).

## Key Algorithms

- **Bresenham line rasterization**
  ([`bresenham_line`](https://zhefrench.github.io/Battlefield/reference/bresenham_line.md)):
  Traces discrete paths on a spot grid.

- **k-NN neighborhood detection** (via RANN): Identifies spatial
  neighbors and cluster composition.

- **Grid type inference**
  ([`detect_grid_type`](https://zhefrench.github.io/Battlefield/reference/detect_grid_type.md)):
  Distinguishes square vs. hexagonal layouts.

- **Layer classification**: Uses radial distance quantiles to stratify
  clusters into border/intermediate/core regions.

- **Centroid-based trajectories**: Connects cluster centroids and
  selects nearby spots for trajectory analysis.

## Important Data Shapes

All selection functions work with `data.frame`s containing:

- **Input columns (required)**:

  `x`, `y` (numeric coordinates), `cluster` (any grouping type),
  `spot_id` (must match `colnames(spe)` for SPE integration).

- **Output columns (added by selection functions)**:

  `interface`, `directed_pair` (e.g., "3-4"), `is_border`, `is_core`
  (logical flags), `layer` ("border"/"intermediate"/"core"). Keep these
  column names stable across pipelines.

## Example Datasets

Three simulated `SpatialExperiment` objects are included for testing:

- [`visium_simulated_spe`](https://zhefrench.github.io/Battlefield/reference/visium_simulated_spe.md):
  Square grid (Visium-like).

- [`visiumHD_8um_simulated_spe`](https://zhefrench.github.io/Battlefield/reference/visiumHD_8um_simulated_spe.md):
  Hexagonal grid (VisiumHD 8 µm binned).

- [`visiumHD_16um_simulated_spe`](https://zhefrench.github.io/Battlefield/reference/visiumHD_16um_simulated_spe.md):
  Hexagonal grid (VisiumHD 16 µm binned).

Load with: `data("visium_simulated_spe", package = "Battlefield")`

## Links & References

- Main vignette: `vignettes/Battlefield-Main.Rmd` — comprehensive
  examples for all four core workflows.

- GitHub: <https://github.com/ZheFrench/Battlefield>

## See also

Useful links:

- <https://github.com/ZheFrench/BattleField>

- <https://zhefrench.github.io/Battlefield/>

- Report bugs at <https://github.com/ZheFrench/BattleField/issues>

## Author

**Maintainer**: Jean-Philippe Villemin <jpvillemin@gmail.com>
([ORCID](https://orcid.org/0000-0002-1838-5880))

Other contributors:

- European Research Council (ERC-2022) \[funder\]
