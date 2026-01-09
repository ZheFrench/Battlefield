#' Add border or core spot selections to SpatialExperiment colData
#'
#' This function takes spot selection results (either border or core spots)
#' and adds them as annotations in the colData of a SpatialExperiment object.
#'
#' @param spe A SpatialExperiment object from which the selection dataframes
#'were derived.
#' @param border Optional data.frame of border spots. Can be output from either
#'   [select_border_spots()] or [build_all_borders()].
#'   Expected columns: spot_id, interface, mode.
#' @param core Optional data.frame of core spots. Can be output from either
#'   [select_core_spots()] or [build_all_cores()].
#'   Expected columns: spot_id, interface, mode.
#' @param erase Logical. If `TRUE`, erase pre-existing battlefield columns
#'   before adding new ones. If `FALSE` (default), issue a warning if columns
#'   already exist and skip adding them.
#'
#' @details
#' At least `border` must be provided. `core` is optional. If `core` is provided
#' without `border`, an error is raised.
#'
#' The following unified columns are added:
#' - `is_border`: logical, TRUE if spot is a border spot, FALSE or NA otherwise
#' - `is_core`: logical, TRUE if spot is a core spot, FALSE or NA otherwise
#' - `interface`: the target interface cluster
#' - `border_mode`: character, "inner", "outer", or NA (from the `mode` column
#'in input data)
#'
#' @return The SpatialExperiment object with updated colData.
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
#' # Using build_all_borders and build_all_cores
#' all_borders <- build_all_borders(df, k = 6)
#' all_cores <- build_all_cores(df, all_borders, region = "inner")
# spe <- add_borders_to_spe(spe, border = all_borders, core = all_cores)
#'
#' # Or using individual select functions
#' border_3_to_4 <- select_border_spots(df, cluster = 3, interface = 4, k = 6)
#' core_3_to_4 <- select_core_spots(df, all_borders, cluster = 3, interface = 4)
# spe <- add_borders_to_spe(spe, border = border_3_to_4, core = core_3_to_4)
#'
#' @importFrom SummarizedExperiment colData "colData<-"
#' @importFrom methods is
#' @export
add_borders_to_spe <- function(spe,
                                border = NULL,
                                core = NULL,
                                erase = FALSE) {

# Validate input
stopifnot(methods::is(spe, "SpatialExperiment"))

if (is.null(border)) {
    stop("'border' must be provided. Cannot add only core spots without border spots.")}

# Get spot IDs
if ("spot_id" %in% colnames(colData(spe))) {
    spot_ids <- colData(spe)$spot_id
} else {
    spot_ids <- colnames(spe)
}

cd <- colData(spe)

# Check for pre-existing battlefield columns
battlefield_cols <- c("is_border", "is_core", "interface", "border_mode")
existing_battlefield_cols <- battlefield_cols[battlefield_cols %in% colnames(cd)]
if (length(existing_battlefield_cols) > 0) {
    if (!erase) {
    message("Battlefield columns already exist: ",
    paste(existing_battlefield_cols, collapse = ", "),            ". Set erase=TRUE to overwrite them.")
    border <- NULL
    core <- NULL
    } else {
    cd <- cd[, !(colnames(cd) %in% existing_battlefield_cols), drop = FALSE]
    message("Erasing pre-existing battlefield columns.")
    }
}

# Initialize unified battlefield columns for all spots
is_border_col <- rep(FALSE, length(spot_ids))
is_core_col <- rep(FALSE, length(spot_ids))
interface_col <- rep(NA_character_, length(spot_ids))
border_mode_col <- rep(NA_character_, length(spot_ids))

# Add border data
if (!is.null(border) && nrow(border) > 0) {
    idx <- match(border$spot_id, spot_ids)
    valid_idx <- !is.na(idx)

    is_border_col[idx[valid_idx]] <- TRUE

    if ("interface" %in% colnames(border)) {
    interface_col[idx[valid_idx]] <- as.character(border$interface[valid_idx])
    }

    if ("mode" %in% colnames(border)) {
    border_mode_col[idx[valid_idx]] <- as.character(border$mode[valid_idx])
    }
}

# Add core data
if (!is.null(core) && nrow(core) > 0) {
    idx <- match(core$spot_id, spot_ids)
    valid_idx <- !is.na(idx)

    # Only set to core if not already marked as border
    is_core_only <- is.na(is_border_col[idx[valid_idx]])
    is_core_col[idx[valid_idx][is_core_only]] <- TRUE

    # Only set interface/mode for spots marked as core (not already border)
    if ("interface" %in% colnames(core)) {
    interface_col[idx[valid_idx][is_core_only]] <- as.character(core$interface[valid_idx][is_core_only])    }

    if ("mode" %in% colnames(core)) {
    border_mode_col[idx[valid_idx][is_core_only]] <- as.character(core$mode[valid_idx][is_core_only])    }
}

# Add columns to colData
cd$is_border <- is_border_col
cd$is_core <- is_core_col
cd$interface <- interface_col
cd$border_mode <- border_mode_col

# Update SPE with new colData
colData(spe) <- cd

spe
}

#' Add layer classifications to SpatialExperiment colData
#'
#' This function takes layer classification results (border, intermediate, core)
#' and adds them as annotations in the colData of a SpatialExperiment object.
#'
#' @param spe A SpatialExperiment object from which the layer dataframes were
#'derived.
#' @param layer A data.frame of layer classifications. Can be output from either
#'   [create_cluster_layers()] or [create_all_layers()].
#'   Expected columns: spot_id, layer.
#' @param erase Logical. If `TRUE`, erase pre-existing layer columns
#'   before adding new ones. If `FALSE` (default), issue a warning if columns
#'   already exist and skip adding them.
#'
#' @details
#' The `layer` parameter must be provided and contain at least one row.
#'
#' The following column is added:
#' - `layer`: character, one of "border", "intermediate", "core", or NA
#'
#' @return The SpatialExperiment object with updated colData.
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
#' # Or for a single cluster
#' cluster_1_layers <- create_cluster_layers(df, target_cluster = 1, k = 6)
#' spe <- add_layers_to_spe(spe, layer = cluster_1_layers)
#'
#' @importFrom SummarizedExperiment colData "colData<-"
#' @importFrom methods is
#' @export
add_layers_to_spe <- function(spe,
                            layer = NULL,
                            erase = FALSE) {

# Validate input
stopifnot(methods::is(spe, "SpatialExperiment"))

if (is.null(layer)) {
    stop("'layer' must be provided.")
}

stopifnot(is.data.frame(layer))
stopifnot("spot_id" %in% colnames(layer))
stopifnot("layer" %in% colnames(layer))

# Get spot IDs
if ("spot_id" %in% colnames(colData(spe))) {
    spot_ids <- colData(spe)$spot_id
} else {
    spot_ids <- colnames(spe)
}

cd <- colData(spe)

# Check for pre-existing layer column
if ("layer" %in% colnames(cd)) {
    if (!erase) {
    message("Layer column already exists. Set erase=TRUE to overwrite it.")
    return(spe)
    } else {
    cd <- cd[, colnames(cd) != "layer", drop = FALSE]
    message("Erasing pre-existing layer column.")
    }
}

# Initialize layer column for all spots
layer_col <- rep(NA_character_, length(spot_ids))

# Add layer data
if (nrow(layer) > 0) {
    idx <- match(layer$spot_id, spot_ids)
    valid_idx <- !is.na(idx)

    if (any(valid_idx)) {
    layer_col[idx[valid_idx]] <- as.character(layer$layer[valid_idx])
    }
}

# Add layer column to colData
cd$layer <- layer_col

# Update SPE with new colData
colData(spe) <- cd

spe
}

#' Add trajectory information to SpatialExperiment colData
#'
#' This function takes trajectory results (e.g., from
#'[build_similar_trajectories()])
#' and adds trajectory metadata to the colData of a SpatialExperiment object.
#'
#' @param spe A SpatialExperiment object from which the trajectory dataframes
#'were derived.
#' @param trajectory A data.frame of trajectory spots. Expected output from
#' [build_similar_trajectories()]. Expected columns: spot_id, trajectory_id,
#'offset, pos_on_seg, dist_to_seg.
#' @param erase Logical. If `TRUE`, erase pre-existing trajectory columns
#'   before adding new ones. If `FALSE` (default), issue a warning if columns
#'   already exist and skip adding them.
#'
#' @details
#' The `trajectories` parameter must be provided and contain at least one row.
#'
#' The following columns are added:
#' - `trajectory_id`: character, identifier for each trajectory (e.g., "main",
#'"left_1", "right_2")
#' - `offset`: numeric, the offset distance used to build each trajectory line
#' - `pos_on_seg`: numeric in [0, 1], the position along the trajectory
#' - `dist_to_seg`: numeric, the distance to the original segment
#'
#' @return The SpatialExperiment object with updated colData.
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
#' centroids <- compute_centroids(df)
#' A <- centroids[centroids$cluster == 1, c("x", "y")]
#' B <- centroids[centroids$cluster == 3, c("x", "y")]
#' trajs <- build_similar_trajectories(df, A, B, top_n = 10, n_extra = 1, side =
#'"both")
#' spe <- add_trajectories_to_spe(spe, trajectory = trajs)
#'
#' @importFrom SummarizedExperiment colData "colData<-"
#' @importFrom methods is
#' @export
add_trajectories_to_spe <- function(spe,
                                    trajectory = NULL,
                                    erase = FALSE) {

# Validate input
stopifnot(methods::is(spe, "SpatialExperiment"))

if (is.null(trajectory)) {
    stop("'trajectory' must be provided.")
}

stopifnot(is.data.frame(trajectory))
stopifnot("spot_id" %in% colnames(trajectory))
stopifnot("trajectory_id" %in% colnames(trajectory))

# Get spot IDs
if ("spot_id" %in% colnames(colData(spe))) {
    spot_ids <- colData(spe)$spot_id
} else {
    spot_ids <- colnames(spe)
}

cd <- colData(spe)

# Check for pre-existing trajectory columns
trajectory_cols <- c("trajectory_id", "offset", "pos_on_seg", "dist_to_seg")
existing_trajectory_cols <- trajectory_cols[trajectory_cols %in% colnames(cd)]

if (length(existing_trajectory_cols) > 0) {
    if (!erase) {
    message("Trajectory columns already exist: ",
    paste(existing_trajectory_cols, collapse = ", "),            ". Set erase=TRUE to overwrite them.")
    return(spe)
    } else {
    cd <- cd[, !(colnames(cd) %in% existing_trajectory_cols), drop = FALSE]
    message("Erasing pre-existing trajectory columns.")
    }
}

# Initialize trajectory columns for all spots
trajectory_id_col <- rep(NA_character_, length(spot_ids))
offset_col <- rep(NA_real_, length(spot_ids))
pos_on_seg_col <- rep(NA_real_, length(spot_ids))
dist_to_seg_col <- rep(NA_real_, length(spot_ids))

# Add trajectory data
if (nrow(trajectory) > 0) {
    idx <- match(trajectory$spot_id, spot_ids)
    valid_idx <- !is.na(idx)

    if (any(valid_idx)) {
    trajectory_id_col[idx[valid_idx]] <- as.character(trajectory$trajectory_id[valid_idx])
    if ("offset" %in% colnames(trajectory)) {
        offset_col[idx[valid_idx]] <- as.numeric(trajectory$offset[valid_idx])
    }

    if ("pos_on_seg" %in% colnames(trajectory)) {
        pos_on_seg_col[idx[valid_idx]] <- as.numeric(trajectory$pos_on_seg[valid_idx])    }

    if ("dist_to_seg" %in% colnames(trajectory)) {
        dist_to_seg_col[idx[valid_idx]] <- as.numeric(trajectory$dist_to_seg[valid_idx])    }
    }
}

# Add columns to colData
cd$trajectory_id <- trajectory_id_col
cd$offset <- offset_col
cd$pos_on_seg <- pos_on_seg_col
cd$dist_to_seg <- dist_to_seg_col

# Update SPE with new colData
colData(spe) <- cd

spe
}
