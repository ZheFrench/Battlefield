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
#' @param center Optional spot used as the center of the core: its identifier,
#'   or a one-row data.frame with its exact coordinates (e.g., `A` copied from
#'   [view_spots()]). It must be an existing non-border spot of
#'   `target_cluster` that does not touch the border layer either. If `NULL`
#'   (default), layers are defined by the distance to the border.
#' @param core_radius Optional numeric, only used with `center`. Maximum
#'   distance to `center` of core spots, in coordinate units. If `NULL`
#'   (default), it is set from `intermediate_quantile`.
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
#' Spots touching a border spot (kNN-based, in either direction) are always
#' intermediate, so that core and border are never in contact. Narrow clusters
#' may thus have no core.
#'
#' If `center` is provided, border spots are unchanged, but core spots are the
#' non-border spots closest to `center`: those within `core_radius` or, by
#' default, the closest `1 - intermediate_quantile` proportion.
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
#'   x = SpatialExperiment::spatialCoords(spe)[, 1],
#'   y = SpatialExperiment::spatialCoords(spe)[, 2],
#'   cluster = SummarizedExperiment::colData(spe)$cluster
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
#' # Core centered on a chosen spot: here the core spot closest to the
#' # cluster centroid
#' inner <- layers[layers$layer == "core", ]
#' center <- closest_spot(inner, mean(inner$x), mean(inner$y))$spot_id
#' layers_center <- create_cluster_layers(df, target_cluster = 1, k = 6,
#' center = center)
#' table(layers_center$layer)
#'
#' @importFrom RANN nn2
#' @export
create_cluster_layers <- function(df,
                                target_cluster,
                                k = 6,
                                max_dist = NULL,
                                intermediate_quantile = 0.5,
                                coord_cols = c("x", "y"),
                                cluster_col = "cluster",
                                center = NULL,
                                core_radius = NULL) {

stopifnot(all(coord_cols %in% colnames(df)))
stopifnot(cluster_col %in% colnames(df))
stopifnot(intermediate_quantile >= 0 && intermediate_quantile <= 1)

# Filter to target cluster
idx_target <- df[[cluster_col]] == target_cluster
center_row <- NULL
if (!is.null(core_radius)) {
    if (is.null(center)) stop("'core_radius' requires 'center'.")
    stopifnot(is.numeric(core_radius), length(core_radius) == 1,
    !is.na(core_radius), core_radius >= 0)
}
if (!is.null(center)) {
    center_row <- .spot_row(df, center, coord_cols, "center")
    if (!isTRUE(idx_target[center_row])) {
    stop("'center' is not in cluster ", target_cluster, ".")
    }
}
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

if (!is.null(center_row) && touches_other[match(center_row, idx_from)]) {
    stop("'center' is a border spot of cluster ", target_cluster,
    ": choose a spot that does not touch another cluster.")
}

# Identify intermediate spots: spots neighboring border spots but not themselves
# borders
# Get all border spots (including those from other clusters)
all_border_idx <- which(touches_other)
border_spot_indices_global <- idx_from[all_border_idx]

# For non-border spots in target cluster, check if they neighbor border spots
non_border_idx <- idx_from[!touches_other]

# Spots touching a border spot, in either kNN direction: never core
near <- nn_idx
if (!is.null(max_dist)) near[nn_dist > max_dist] <- NA
touches_border <- !touches_other &
    (rowSums(matrix(near %in% idx_from[touches_other], nrow(near))) > 0 |
    idx_from %in% near[touches_other, ])

if (!is.null(center_row)) {
    center_pos <- match(center_row, idx_from)
    if (touches_border[center_pos]) {
    stop("'center' touches the border layer of cluster ", target_cluster,
        ": choose a spot further inside.")
    }
    # Core = other spots closest to the center
    delta <- sweep(coords[idx_from, , drop = FALSE], 2, coords[center_row, ])
    dist_center <- sqrt(rowSums(delta^2))
    if (is.null(core_radius)) {
    core_radius <- stats::quantile(dist_center[!touches_other],
        1 - intermediate_quantile)
    }
    target_df$layer <- ifelse(touches_other, "border",
    ifelse(dist_center <= core_radius & !touches_border, "core",
        "intermediate"))
    return(target_df)
}

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
    # An intermediate spot always separates core from border
    target_df$layer[touches_border] <- "intermediate"

    # No core: use the intermediate spots that do not touch the border, if any
    if (!any(target_df$layer == "core")) {
    inner <- target_df$layer == "intermediate" & !touches_border
    target_df$layer[inner] <- "core"
    if (any(inner)) {
        message("Cluster ", target_cluster, ": Only 'border' and 'intermediate' layers detected. ",
        "Setting 'intermediate' spots not touching the border to 'core'.")
    } else {
        message("Cluster ", target_cluster, ": No 'core' layer, all inner spots touch the border.")
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
#'   x = SpatialExperiment::spatialCoords(spe)[, 1],
#'   y = SpatialExperiment::spatialCoords(spe)[, 2],
#'   cluster = SummarizedExperiment::colData(spe)$cluster
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
