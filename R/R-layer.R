#' Create layer classification for spots in a cluster
#'
#' This function classifies spots within a target cluster into three layers:
#' core, intermediate, and border. The border layer consists of spots at the
#' cluster interface (touching other clusters). The core layer consists of spots
#' far from any border. The intermediate layer bridges the two, with depth
#' defined by proximity to the border.
#'
#' @param df A data.frame containing at least the coordinate columns and a
#' cluster label column.
#' @param target_cluster Cluster label for which to create layers.
#' @param k Integer. Number of nearest neighbors to consider (excluding self)
#' when
#'   identifying border spots. Default is 6.
#' @param max_dist Optional numeric. Maximum Euclidean distance for neighbors.
#'   If `NULL`, no distance filtering is applied.
#' @param intermediate_quantile Numeric between 0 and 1. The quantile of
#' distances to border
#' used to define the intermediate layer threshold. Default is 0.5 (median).
#' #' Spots with
#' distance to border <= this quantile threshold are classified as intermediate.
#' Use lower values (e.g., 0.33) for narrower intermediate layers, higher values
#' (e.g., 0.75)
#'   for wider intermediate layers.
#' @param coord_cols Character vector of length 2 giving the coordinate column
#'names.
#'   Default is `c("x","y")`.
#' @param cluster_col Character. Name of the column containing cluster labels.
#'   Default is `"cluster"`.
#'
#' @details
#' Layer definitions:
#' - **Border**: Spots in the target cluster that touch at least one spot from a
#' different cluster (kNN-based).
#' - **Intermediate**: Non-border spots whose distance to the nearest border
#' spot is within
#'   the `intermediate_quantile` threshold of all non-border spots.
#' - **Core**: All remaining non-border spots that are far from the border.
#'
#' The `intermediate_quantile` parameter controls the depth of the intermediate
#' layer:
#' - 0.33: Narrow intermediate layer (only very close to border)
#' - 0.50: Moderate intermediate layer (median distance threshold)
#' - 0.75: Wide intermediate layer (most non-border spots included)
#'
#' @return A data.frame based on `df` filtered to contain only spots from
#' `target_cluster`, with an additional column:
#' - `layer`: character, one of "border", "intermediate", or "core"
#'
#' @examples
#' data("visiumHD_16um_simulated_spe", package = "Battlefield")
#' spe <- visiumHD_16um_simulated_spe
#' df <- data.frame(
#'   spot_id = colnames(spe),
#'   x = spatialCoords(spe)[, 1],
#'   y = spatialCoords(spe)[, 2],
#'   cluster = colData(spe)$cluster
#' )
#' # Default: moderate depth
#' layers <- create_cluster_layers(df, target_cluster = 1, k = 6)
#' table(layers$layer)
#'
#' # Narrow intermediate layer
#' layers_narrow <- create_cluster_layers(df, target_cluster = 1, k = 6,
#' intermediate_quantile = 0.33)
#' table(layers_narrow$layer)
#'
#' @importFrom RANN nn2
#' @export
create_cluster_layers <- function(df,
                                target_cluster,
                                k = 6,
                                max_dist = NULL,
                                intermediate_quantile = 0.5,
                                coord_cols = c("x", "y"),
                                cluster_col = "cluster") {

stopifnot(all(coord_cols %in% colnames(df)))
stopifnot(cluster_col %in% colnames(df))
stopifnot(intermediate_quantile >= 0 && intermediate_quantile <= 1)

# Filter to target cluster
idx_target <- df[[cluster_col]] == target_cluster
target_df <- df[idx_target, , drop = FALSE]

if (nrow(target_df) == 0) {
    stop("No spots found in target_cluster.")
}

coords <- as.matrix(df[, coord_cols])
cl <- df[[cluster_col]]

# Identify border spots: spots in target cluster touching other clusters
idx_from <- which(idx_target)
nn <- RANN::nn2(data = coords, query = coords[idx_from, , drop = FALSE],
k = k + 1)
nn_idx  <- nn$nn.idx[, -1, drop = FALSE]     # remove self
nn_dist <- nn$nn.dists[, -1, drop = FALSE]

neigh_cl <- matrix(cl[nn_idx], nrow = nrow(nn_idx))

# Optional distance filter
if (!is.null(max_dist)) {
    neigh_cl[nn_dist > max_dist] <- NA
}

# A spot is a border spot if it touches at least one spot from a different
# cluster
touches_other <- apply(!(neigh_cl == target_cluster) & !is.na(neigh_cl),
                        1, any)

border_spot_ids <- target_df$spot_id[touches_other]

# Identify intermediate spots: spots neighboring border spots but not themselves
# borders
# Get all border spots (including those from other clusters)
all_border_idx <- which(touches_other)
border_spot_indices_global <- idx_from[all_border_idx]

# For non-border spots in target cluster, check if they neighbor border spots
non_border_idx <- idx_from[!touches_other]

if (length(non_border_idx) == 0) {
    # All spots are borders
    target_df$layer <- "border"
    message("Cluster ", target_cluster,
    ": Only 'border' layer detected. All spots are at cluster interfaces.")} else if (length(border_spot_indices_global) == 0) {
    # No border spots found (entire cluster is homogeneous)
    target_df$layer <- "core"
    message("Cluster ", target_cluster,
    ": Only 'core' layer detected. No cluster interfaces found.")} else {
    # Query non-border spots for neighbors in border region
    nn_check <- RANN::nn2(data = coords[border_spot_indices_global, , drop = FALSE],                        query = coords[non_border_idx, , drop = FALSE],
                        k = 1)

    # Calculate depth: distance to nearest border spot
    distances_to_border <- nn_check$nn.dists[, 1]

    # Use quantile to define intermediate layer threshold
    depth_threshold <- stats::quantile(distances_to_border,
    intermediate_quantile)
    is_intermediate <- distances_to_border <= depth_threshold

    # Initialize layer column
    target_df$layer <- "core"
    target_df$layer[touches_other] <- "border"
    non_border_rows <- which(!touches_other)
    target_df$layer[non_border_rows[is_intermediate]] <- "intermediate"

    # Check if we only have 2 layers (border and intermediate, no core)
    if (sum(target_df$layer == "core") == 0) {
    # Verify we have both border and intermediate (not just border)
    if (sum(target_df$layer == "border") > 0 && sum(target_df$layer == "intermediate") > 0) {        # Transform intermediate spots to core
        target_df$layer[target_df$layer == "intermediate"] <- "core"
        message("Cluster ", target_cluster, ": Only 'border' and 'intermediate' layers detected. ",                "Setting 'intermediate' spots to 'core'.")
    }
    }
}

target_df
}

#' Create and visualize layers for multiple clusters
#'
#' Batch create layer classifications for multiple clusters and optionally
#' generate individual plots for each cluster.
#'
#' @param df A data.frame containing at least coordinate and cluster columns.
#' @param clusters Optional vector of cluster labels to process. If `NULL`,
#'   all unique clusters are used.
#' @param k Integer. Number of nearest neighbors to consider. Default is 6.
#' @param max_dist Optional numeric. Maximum Euclidean distance for neighbors.
#' @param intermediate_quantile Numeric between 0 and 1. The quantile of
#' distances to border
#'   used to define the intermediate layer threshold. Default is 0.5 (median).
#' @param coord_cols Character vector of length 2 giving the coordinate column
#' names.
#'   Default is `c("x","y")`.
#' @param cluster_col Character. Name of the column containing cluster labels.
#'   Default is `"cluster"`.
#'
#' @return A data.frame combining all processed clusters with the `layer` column
#' added.
#'
#' @examples
#' data("visiumHD_16um_simulated_spe", package = "Battlefield")
#' spe <- visiumHD_16um_simulated_spe
#' df <- data.frame(
#'   spot_id = colnames(spe),
#'   x = spatialCoords(spe)[, 1],
#'   y = spatialCoords(spe)[, 2],
#'   cluster = colData(spe)$cluster
#' )
#' all_layers <- create_all_layers(df, k = 6, intermediate_quantile = 0.5)
#' head(all_layers)
#'
#' @importFrom dplyr bind_rows
#' @export
create_all_layers <- function(df,
                            clusters = NULL,
                            k = 6,
                            max_dist = NULL,
                            intermediate_quantile = 0.5,
                            coord_cols = c("x", "y"),
                            cluster_col = "cluster") {

stopifnot(all(coord_cols %in% colnames(df)))
stopifnot(cluster_col %in% colnames(df))

if (is.null(clusters)) {
    clusters <- unique(df[[cluster_col]])
    clusters <- clusters[!is.na(clusters)]
}

res <- lapply(clusters, function(c) {
    create_cluster_layers(
    df,
    target_cluster = c,
    k = k,
    max_dist = max_dist,
    intermediate_quantile = intermediate_quantile,
    coord_cols = coord_cols,
    cluster_col = cluster_col
    )
})

dplyr::bind_rows(res)
}
