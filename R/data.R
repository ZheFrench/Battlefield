#' Simulated Visium SpatialExperiment dataset
#'
#' A simulated \code{SpatialExperiment} object mimicking a 10x Genomics
#' Visium spatial transcriptomics layout. The dataset contains spot-level
#' spatial coordinates and cluster annotations, and is intended for
#' testing.
#'
#' @format
#' A \code{\link[SpatialExperiment]{SpatialExperiment}} object with:
#' \describe{
#'   \item{assays}{One dummy assay (required container structure)}
#'   \item{colData}{Spot metadata including \code{barcode} and
#'   \code{cluster}}
#'   \item{spatialCoords}{Numeric matrix ofx/y spot coordinates}
#' }
#'
#' @details
#' This dataset  contain only gene expression values for one gene
#' named FAKE_GENE.
#' The dummy assay is included solely to satisfy the
#' \code{SummarizedExperiment} container requirements.
#'
#' Spot coordinates follow a Visium-like hexagonal grid.
#' Cluster labels are simulated and have no biological meaning.
#'
#' @source
#' Simulated internally for package development.
#'
#' @usage data(visium_simulated_spe)
#'
#' @examples
#' data(visium_simulated_spe)
#'
#' visium_simulated_spe
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[SpatialExperiment]{SpatialExperiment}}
#'   \item \code{\link{compute_centroids}}
#' }
#'
#' visium_simulated_spe
"visium_simulated_spe"

#' Simulated VisiumHD 8 µm binned SpatialExperiment dataset
#'
#' A simulated \code{SpatialExperiment} object mimicking a 10x Genomics
#' Visium spatial transcriptomics layout. The dataset contains spot-level
#' spatial coordinates and cluster annotations, and is intended for
#' testing.
#'
#' @format
#' A \code{\link[SpatialExperiment]{SpatialExperiment}} object with:
#' \describe{
#'   \item{assays}{One dummy assay (required container structure)}
#'   \item{colData}{Spot metadata including \code{barcode} and
#'   \code{cluster}}
#'   \item{spatialCoords}{Numeric matrix of x/y spot coordinates}
#' }
#'
#' @details
#' This dataset  contain only gene expression values for one gene
#' named FAKE_GENE.
#' The dummy assay is included solely to satisfy the
#' \code{SummarizedExperiment} container requirements.
#'
#' Spot coordinates follow a Visium-like squared grid.
#' Cluster labels are simulated and have no biological meaning.
#' 
#' The layout is intended to mimic a Visium HD *binned* layer
#' (e.g., 8 µm bin size), where each spot represents one bin
#' arranged on an adjacent square grid.
#' @source
#' Simulated internally for package development.
#'
#' @usage data(visiumHD_8um_simulated_spe)
#'
#' @examples
#' data(visiumHD_8um_simulated_spe)
#'
#' visiumHD_8um_simulated_spe
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[SpatialExperiment]{SpatialExperiment}}
#'   \item \code{\link{compute_centroids}}
#' }
#'
#' visiumHD_8um_simulated_spe
"visiumHD_8um_simulated_spe"

#' Simulated VisiumHD 16 µm binned SpatialExperiment dataset
#'
#' A simulated \code{SpatialExperiment} object mimicking a 10x Genomics
#' Visium spatial transcriptomics layout. The dataset contains spot-level
#' spatial coordinates and cluster annotations, and is intended for
#' testing.
#'
#' @format
#' A \code{\link[SpatialExperiment]{SpatialExperiment}} object with:
#' \describe{
#'   \item{assays}{One dummy assay (required container structure)}
#'   \item{colData}{Spot metadata including \code{barcode} and
#'   \code{cluster}}
#'   \item{spatialCoords}{Numeric matrix of x/y spot coordinates}
#' }
#'
#' @details
#' This dataset  contain only gene expression values for one gene
#' named FAKE_GENE.
#' The dummy assay is included solely to satisfy the
#' \code{SummarizedExperiment} container requirements.
#'
#' Spot coordinates follow a Visium-like squared grid.
#' Cluster labels are simulated and have no biological meaning.
#' 
#' The layout is intended to mimic a Visium HD *binned* layer
#' (e.g., 16 µm bin size), where each spot represents one bin
#' arranged on an adjacent square grid.
#' @source
#' Simulated internally for package development.
#'
#' @usage data(visiumHD_16um_simulated_spe)
#'
#' @examples
#' data(visiumHD_16um_simulated_spe)
#'
#' visiumHD_16um_simulated_spe
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[SpatialExperiment]{SpatialExperiment}}
#'   \item \code{\link{compute_centroids}}
#' }
#'
#' visiumHD_16um_simulated_spe
"visiumHD_16um_simulated_spe"