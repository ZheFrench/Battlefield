# Simulated VisiumHD 8 µm binned SpatialExperiment dataset

A simulated `SpatialExperiment` object mimicking a 10x Genomics Visium
spatial transcriptomics layout. The dataset contains spot-level spatial
coordinates and cluster annotations, and is intended for testing.

## Usage

``` r
data(visiumHD_simulated_spe)
```

## Format

A
[`SpatialExperiment`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html)
object with:

- assays:

  One dummy assay (required container structure)

- colData:

  Spot metadata including `barcode` and `cluster`

- spatialCoords:

  Numeric matrix of x/y spot coordinates

## Source

Simulated internally for package development.

## Details

This dataset contain only gene expression values for one gene named
FAKE_GENE. The dummy assay is included solely to satisfy the
`SummarizedExperiment` container requirements.

Spot coordinates follow a Visium-like squared grid. Cluster labels are
simulated and have no biological meaning.

The layout is intended to mimic a Visium HD \*binned\* layer (e.g., 8 µm
bin size), where each spot represents one bin arranged on an adjacent
square grid.

## See also

- [`SpatialExperiment`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html)

- [`compute_centroids`](https://zhefrench.github.io/Battlefield/reference/compute_centroids.md)

"visiumHD_simulated_spe"

## Examples

``` r
data(visiumHD_simulated_spe)
#> Warning: data set ‘visiumHD_simulated_spe’ not found

visiumHD_simulated_spe
#> Error: object 'visiumHD_simulated_spe' not found
```
