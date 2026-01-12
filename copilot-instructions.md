# Guidance for AI coding agents — Battlefield

This file captures the minimal, actionable knowledge an AI coding agent
needs to be productive in the Battlefield R package repository.

## 1. Purpose & big picture

- **Battlefield**: low-level toolkit to select spatial spots (fronts,
  margins, trajectories) on spot-based spatial transcriptomics data
  (Visium / VisiumHD).
- **Key data flow**: raw `SpatialExperiment` / binned VisiumHD imports →
  cluster labels in `colData()` → spatial coordinate tables (`x`,`y`) →
  selection helpers in `R/` produce ordered spot lines and interfaces
  used by vignettes and downstream analyses.
- **Core use case**: Identify tissue regions (invasive margins, niche
  borders, cluster interfaces) for downstream analyses like
  ligand–receptor studies (e.g., ligand-receptor interactions via
  BulkSignalR).
- **Key algorithms**: Bresenham line rasterization, k-NN neighborhood
  detection, grid type inference (square vs. hexagonal), layer
  classification, trajectory building via centroid-to-centroid paths.
- **R version requirement**: R ≥ 4.5 (base pipe `|>` syntax required).

## 2. Important files & locations

- **Package metadata**:
  [DESCRIPTION](https://zhefrench.github.io/Battlefield/DESCRIPTION) (R
  \>= 4.5, Bioconductor package dependencies).
- **Core algorithm modules**:
  - [R/R-trajectory.R](https://zhefrench.github.io/Battlefield/R/R-trajectory.R)
    — geometry primitives (Bresenham), centroids, trajectory building.
  - [R/R-grid.R](https://zhefrench.github.io/Battlefield/R/R-grid.R) —
    grid type detection (square vs hexagonal), spacing estimation.
  - [R/R-border.R](https://zhefrench.github.io/Battlefield/R/R-border.R)
    — cluster interface/border/core spot selection logic.
  - [R/R-layer.R](https://zhefrench.github.io/Battlefield/R/R-layer.R) —
    layer classification (border/intermediate/core) within clusters.
  - [R/R-neighbourhood.R](https://zhefrench.github.io/Battlefield/R/R-neighbourhood.R)
    — k-NN neighborhood composition analysis.
  - [R/R-spe.R](https://zhefrench.github.io/Battlefield/R/R-spe.R) —
    SpatialExperiment integration, adding selections to colData.
- **Data**: [R/data.R](https://zhefrench.github.io/Battlefield/R/data.R)
  (documents included datasets: `visium_simulated_spe`,
  `visiumHD_8um_simulated_spe`, `visiumHD_16um_simulated_spe`).
- **Tests**:
  [tests/testthat](https://zhefrench.github.io/Battlefield/tests/testthat)
  (uses testthat edition 3; organized by module: trajectory, grid,
  frontline, layer, neighborhood, spe).
- **Developer scripts**:
  - [forTest/R-build.R](https://zhefrench.github.io/Battlefield/forTest/R-build.R)
    — full build, document, test, BiocCheck flow.
  - [forTest/R-visHD.R](https://zhefrench.github.io/Battlefield/forTest/R-visHD.R)
    — realistic end-to-end VisiumHD import & selection example.
  - [forTest/R-test-\*.R](https://zhefrench.github.io/Battlefield/forTest/)
    — interactive exploration scripts per module (R-test-traj.R,
    R-test-frontline.R, etc.).
- **Main vignette**:
  [vignettes/Battlefield-Main.Rmd](https://zhefrench.github.io/Battlefield/vignettes/Battlefield-Main.Rmd)
  — primary user documentation with examples for all four core
  workflows.

## 3. Project-specific conventions & patterns

- **Pipe & verbs**: Code uses the base pipe `|>` (R 4.5+) and `dplyr`
  verbs (`mutate`, `filter`, `arrange`, `bind_rows`, `group_by`,
  `slice_*`, `pull`, `summarise`). Use modern pipe syntax when editing;
  avoid magrittr `%>%`.
- **Data shapes**: All selection functions work with `data.frame`s
  containing:
  - Spatial coords: `x`, `y` (numeric; required for geometry helpers)
  - Cluster identity: `cluster` column (factor/character/numeric;
    grouped by cluster)
  - Spot identity: `spot_id` (optional; used when integrating with
    SpatialExperiment)
  - Selection metadata: `interface`, `mode` (added by frontline/layer
    functions)
  - Layer metadata: `layer` (added by layer functions; values: “border”,
    “intermediate”, “core”)
  - Directed pair: `directed_pair` (e.g., “3-4” added by
    [`directed_cluster_interface_pairs()`](https://zhefrench.github.io/Battlefield/reference/directed_cluster_interface_pairs.md))
  - Keep these column names stable across refactoring; many helpers
    assume exact names.
- **Keying for deduplication**: Coordinate deduplication uses
  `.xy_key <- paste0(round(x), "_", round(y))`. This key identifies
  spatial duplicates; be very careful when changing rounding logic or
  keying scheme, as intersection detection and nearest-neighbor logic
  depend on consistent keying.
- **Numeric robustness**: Functions defensively check edge cases
  (zero-length segments, empty clusters, insufficient neighbors) via
  [`stopifnot()`](https://rdrr.io/r/base/stopifnot.html). Examples:
  [`unit_normal_left()`](https://zhefrench.github.io/Battlefield/reference/unit_normal_left.md),
  [`build_one_trajectory()`](https://zhefrench.github.io/Battlefield/reference/build_one_trajectory.md),
  [`bresenham_line()`](https://zhefrench.github.io/Battlefield/reference/bresenham_line.md),
  [`detect_grid_type()`](https://zhefrench.github.io/Battlefield/reference/detect_grid_type.md).
  Pattern: Always validate `is.data.frame(df)`, presence of required
  columns, non-empty rows, numeric coords. Preserve these guards when
  refactoring.
- **Roxygen2 workflow**: All exports are declared via roxygen tags.
  After any code changes, **always** run
  [`devtools::document()`](https://devtools.r-lib.org/reference/document.html)
  to regenerate
  [NAMESPACE](https://zhefrench.github.io/Battlefield/NAMESPACE) and man
  pages. Changes to function signatures or `@export` tags require this
  step. Never manually edit NAMESPACE.
- **Grid-aware code**: Functions like
  [`detect_grid_type()`](https://zhefrench.github.io/Battlefield/reference/detect_grid_type.md)
  and
  [`estimate_spot_spacing()`](https://zhefrench.github.io/Battlefield/reference/estimate_spot_spacing.md)
  probe grid geometry. Square grids have ~4 neighbors; hexagonal grids
  have ~6. Some functions behave differently per grid type; always test
  with both square and hexagonal layouts (using included test datasets).
- **Layer depth control**: Layer classification functions use
  `intermediate_quantile` (default 0.5) to control the depth of
  intermediate layers relative to distance to border. Lower values
  (0.33) create narrow intermediate layers; higher values (0.75) create
  wider ones. This is a tunable parameter users may adjust.

## 4. Development & test workflows

### Interactive development (recommended)

In an R session inside the repository:

``` r
devtools::load_all()
devtools::document()
devtools::test()
```

### Full package build & checks

``` r
devtools::build(vignettes = TRUE)
R CMD check <built_tar.gz>
BiocCheck::BiocCheck(<built_tar.gz>)  # used in forTest/R-build.R
```

### Vignettes

``` r
devtools::build_vignettes()  # generates vignettes/ HTML output
browseVignettes("Battlefield")  # view rendered vignettes
```

The
[forTest/R-build.R](https://zhefrench.github.io/Battlefield/forTest/R-build.R)
script automates the full workflow (load, document, test, build, check,
BiocCheck).

## 5. Tests & expectations

- Unit tests live in
  [tests/testthat](https://zhefrench.github.io/Battlefield/tests/testthat).
  They validate geometry (`bresenham_line`), nearest-neighbor spacing
  (`estimate_spot_spacing`), and selection logic
  (`build_similar_trajectories`, `build_one_trajectory`).

- Run tests with
  [`devtools::test()`](https://devtools.r-lib.org/reference/test.html)
  for rapid iteration. For debugging individual tests:

  ``` r
  devtools::load_all()
  testthat::test_file("tests/testthat/test-trajectory.R")
  ```

- CI should run `R CMD check`; tests use testthat edition 3.

- After adding/modifying input validation
  ([`stopifnot()`](https://rdrr.io/r/base/stopifnot.html) calls), write
  corresponding test cases covering empty data frames, missing columns,
  and boundary conditions.

## 6. External integrations & runtime deps

- Key runtime packages referenced in code and forTest scripts:
  `SpatialExperiment`, `VisiumIO`, `RANN`, `dplyr`, `ggplot2`,
  `testthat`. Check
  [DESCRIPTION](https://zhefrench.github.io/Battlefield/DESCRIPTION) for
  the authoritative list.
- Large-data: for VisiumHD imports the repository uses VisiumIO and
  OSF-hosted artifacts; some forTest examples download/untar sample
  data.

## 7. When editing code

- **Data flow preservation**: Functions accept and return `data.frame`s
  with stable column names (`x`, `y`, `cluster`, `spot_id`, `interface`,
  `mode`, `layer`). Keep column order and naming consistent across the
  pipeline to avoid silent failures in downstream consumers.

- **Input validation pattern**: Always include defensive checks
  (`stopifnot`) at function entry points:

  ``` r
  stopifnot(is.data.frame(df))
  stopifnot(all(c("x","y","cluster") %in% names(df)))
  stopifnot(nrow(df) >= 1)
  ```

  This catches user errors early and makes debugging easier.

- **Roxygen docs + regeneration**: After any signature or export
  changes, run
  [`devtools::document()`](https://devtools.r-lib.org/reference/document.html)
  to update
  [NAMESPACE](https://zhefrench.github.io/Battlefield/NAMESPACE) and man
  pages. Failing to do so breaks package checks.

- **Grid-aware updates**: When modifying geometry helpers, test both
  square and hexagonal grids using built-in datasets
  (`visium_simulated_spe`, `visiumHD_8um_simulated_spe`,
  `visiumHD_16um_simulated_spe`).

- **SpatialExperiment integration**: When adding layer classifications,
  use
  [`add_layers_to_spe()`](https://zhefrench.github.io/Battlefield/reference/add_layers_to_spe.md)
  which auto-manages colData columns (`is_layer`). Always keep `spot_id`
  columns aligned with `colnames(spe)`.

- **Run tests**: After edits, run
  [`devtools::test()`](https://devtools.r-lib.org/reference/test.html)
  to catch breakage early. Focus on geometry and boundary-case tests.
  Use
  [`devtools::load_all()`](https://devtools.r-lib.org/reference/load_all.html)
  first to load changes without reinstalling.

## 8. Common gotchas & debugging

- **Silent failures from column mismatches**: Functions like
  [`build_all_borders()`](https://zhefrench.github.io/Battlefield/reference/build_all_borders.md)
  silently drop rows if `spot_id` is missing or has mismatched names.
  Always verify [`colnames()`](https://rdrr.io/r/base/colnames.html) and
  [`nrow()`](https://rdrr.io/r/base/nrow.html) before passing data.frame
  results downstream.
- **RANN k-NN with small clusters**: When a cluster has fewer than `k+1`
  spots,
  [`RANN::nn2()`](https://jefferislab.github.io/RANN/reference/nn2.html)
  behaves unexpectedly. Functions defensively check cluster sizes; if
  tests fail with tiny clusters, verify `k < min_cluster_size`.
- **Roxygen sync issues**: After editing function signatures
  (parameters, return types) or `@export` tags, running
  [`devtools::test()`](https://devtools.r-lib.org/reference/test.html)
  without first running
  [`devtools::document()`](https://devtools.r-lib.org/reference/document.html)
  will fail because the NAMESPACE and `.Rd` files are stale. Always
  document first.
- **Grid type detection edge cases**:
  [`detect_grid_type()`](https://zhefrench.github.io/Battlefield/reference/detect_grid_type.md)
  may return `"unknown"` for sparse or irregular grids. Check diagnostic
  output (`verbose = TRUE`) and test result before downstream code
  branches on grid type.

## 9. Quick pointers for common tasks

- Find geometry helpers: search for `bresenham_line`,
  `build_one_trajectory`, `build_similar_trajectories` inside `R/`.
- Investigate examples: open
  [forTest/R-visHD.R](https://zhefrench.github.io/Battlefield/forTest/R-visHD.R)
  — it demonstrates realistic end-to-end usage loading VisiumHD data and
  running selection routines.
- For interactive debugging: use
  [forTest/R-test-traj.R](https://zhefrench.github.io/Battlefield/forTest/R-test-traj.R),
  [forTest/R-test-frontline.R](https://zhefrench.github.io/Battlefield/forTest/R-test-frontline.R)
  as templates; they load example data and step through workflows.

If any part is unclear (missing external data access, CI steps, or
runtime environment), tell me which area to expand and I will iterate.
