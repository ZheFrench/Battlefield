# Simulated VisiumHD 16 µm binned SpatialExperiment dataset

A simulated `SpatialExperiment` object mimicking a 10x Genomics Visium
spatial transcriptomics layout. The dataset contains spot-level spatial
coordinates and cluster annotations, and is intended for testing.

## Usage

``` r
data(visiumHD_16um_simulated_spe)
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

The layout is intended to mimic a Visium HD \*binned\* layer (e.g., 16
µm bin size), where each spot represents one bin arranged on an adjacent
square grid.

## See also

- [`SpatialExperiment`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html)

- [`compute_centroids`](https://zhefrench.github.io/Battlefield/reference/compute_centroids.md)

"visiumHD_16um_simulated_spe"

## Examples

``` r
data(visiumHD_16um_simulated_spe)

visiumHD_16um_simulated_spe
#> class: SpatialExperiment 
#> dim: 1 1000 
#> metadata(0):
#> assays(1): counts
#> rownames(1): FAKE_GENE
#> rowData names(0):
#> colnames(1000): bin16_1 bin16_41 ... bin16_960 bin16_1000
#> colData names(3): barcode_id cluster sample_id
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
#> imgData names(0):
```
