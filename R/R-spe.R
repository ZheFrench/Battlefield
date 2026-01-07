#' Add border or inner spot selections to SpatialExperiment colData
#'
#' This function takes spot selection results (either border or inner spots)
#' and adds them as annotations in the colData of a SpatialExperiment object.
#'
#' @param spe A SpatialExperiment object from which the selection dataframes were derived.
#' @param border Optional data.frame of border spots. Can be output from either
#'   [select_border_spots()] or [build_all_borders()].
#'   Expected columns: spot_id, interface, is_border_multiple, other_adjacent_borders,
#'   directed_pair, undirected_pair.
#' @param inner Optional data.frame of inner spots. Can be output from either
#'   [select_inner_spots()] or [build_all_inners()].
#'   Expected columns: spot_id, interface, is_inner, mode.
#' @param erase Logical. If `TRUE`, erase pre-existing border or inner columns
#'   before adding new ones. If `FALSE` (default), issue a warning if columns
#'   already exist and skip adding them.
#'
#' @details
#' Either `border` or `inner` (or both) should be provided. If pre-existing
#' selection columns are detected and `erase=FALSE`, a message is printed
#' and those columns are not overwritten. Set `erase=TRUE` to overwrite them.
#'
#' For border spots, the following columns are added:
#' - `is_border`: logical, TRUE if spot is in border
#' - `border_interface`: the target interface cluster
#' - `border_directed_pair`: the directed pair (e.g., "3-4")
#' - `border_undirected_pair`: undirected pair (e.g., "3-4 / 4-3")
#' - `border_multiple`: logical, TRUE if spot touches multiple interfaces
#' - `border_other_adjacent`: other clusters touched (if multiple interface)
#'
#' For inner spots, the following columns are added:
#' - `is_inner`: logical, TRUE if spot is inner control
#' - `inner_interface`: the target interface cluster
#' - `inner_mode`: mode type ("inner", "outer", or "both")
#' - `inner_directed_pair`: the directed pair (if available)
#'
#' @return The SpatialExperiment object with updated colData.
#'
#' @examples
#' library(SpatialExperiment)
#' spe <- visiumHD_16um_simulated_spe
#' df <- data.frame(
#'   spot_id = colnames(spe),
#'   x = spatialCoords(spe)[, 1],
#'   y = spatialCoords(spe)[, 2],
#'   cluster = colData(spe)$cluster
#' )
#' # Using build_all_borders and build_all_inners
#' all_borders <- build_all_borders(df, k = 6)
#' all_inners <- build_all_inners(df, all_borders, region = "inner")
#' spe <- add_selections_to_spe(spe, border = all_borders, inner = all_inners)
#'
#' # Or using individual select functions
#' border_3_to_4 <- select_border_spots(df, cluster = 3, interface = 4, k = 6)
#' inner_3_to_4 <- select_inner_spots(df, all_borders, cluster = 3, interface = 4, region = "inner")
#' spe <- add_selections_to_spe(spe, border = border_3_to_4, inner = inner_3_to_4)
#'
#' @importFrom SpatialExperiment colData
#' @importFrom methods is
#' @export
add_selections_to_spe <- function(spe,
                                   border = NULL,
                                   inner = NULL,
                                   erase = FALSE) {

  # Validate input
  stopifnot(methods::is(spe, "SpatialExperiment"))

  if (is.null(border) && is.null(inner)) {
    stop("At least one of 'border' or 'inner' must be provided.")
  }

  # Get spot IDs
  if ("spot_id" %in% colnames(colData(spe))) {
    spot_ids <- colData(spe)$spot_id
  } else {
    spot_ids <- colnames(spe)
  }

  cd <- colData(spe)

  # Check for pre-existing border columns
  border_cols <- c("is_border", "border_interface", "border_directed_pair",
                   "border_undirected_pair", "border_multiple", "border_other_adjacent")
  existing_border_cols <- border_cols[border_cols %in% colnames(cd)]

  if (length(existing_border_cols) > 0 && !is.null(border)) {
    if (!erase) {
      message("Border columns already exist: ", paste(existing_border_cols, collapse = ", "),
              ". Set erase=TRUE to overwrite them.")
      border <- NULL
    } else {
      cd <- cd[, !(colnames(cd) %in% existing_border_cols), drop = FALSE]
      message("Erasing pre-existing border columns.")
    }
  }

  # Check for pre-existing inner columns
  inner_cols <- c("is_inner", "inner_interface", "inner_mode", "inner_directed_pair")
  existing_inner_cols <- inner_cols[inner_cols %in% colnames(cd)]

  if (length(existing_inner_cols) > 0 && !is.null(inner)) {
    if (!erase) {
      message("Inner columns already exist: ", paste(existing_inner_cols, collapse = ", "),
              ". Set erase=TRUE to overwrite them.")
      inner <- NULL
    } else {
      cd <- cd[, !(colnames(cd) %in% inner_cols), drop = FALSE]
      message("Erasing pre-existing inner columns.")
    }
  }

  # Add border columns
  if (!is.null(border) && nrow(border) > 0) {
    cd$is_border <- spot_ids %in% border$spot_id

    # Add border metadata
    border_interface <- rep(NA_character_, length(spot_ids))
    border_directed_pair <- rep(NA_character_, length(spot_ids))
    border_undirected_pair <- rep(NA_character_, length(spot_ids))
    border_multiple <- rep(FALSE, length(spot_ids))
    border_other_adjacent <- rep(NA_character_, length(spot_ids))

    idx <- match(border$spot_id, spot_ids)
    valid_idx <- !is.na(idx)

    if ("interface" %in% colnames(border)) {
      border_interface[idx[valid_idx]] <- as.character(border$interface[valid_idx])
    }

    if ("directed_pair" %in% colnames(border)) {
      border_directed_pair[idx[valid_idx]] <- as.character(border$directed_pair[valid_idx])
    }

    if ("undirected_pair" %in% colnames(border)) {
      border_undirected_pair[idx[valid_idx]] <- as.character(border$undirected_pair[valid_idx])
    }

    if ("is_border_multiple" %in% colnames(border)) {
      border_multiple[idx[valid_idx]] <- border$is_border_multiple[valid_idx]
    }

    if ("other_adjacent_borders" %in% colnames(border)) {
      border_other_adjacent[idx[valid_idx]] <- as.character(border$other_adjacent_borders[valid_idx])
    }

    cd$border_interface <- border_interface
    cd$border_directed_pair <- border_directed_pair
    cd$border_undirected_pair <- border_undirected_pair
    cd$border_multiple <- border_multiple
    cd$border_other_adjacent <- border_other_adjacent
  }

  # Add inner columns
  if (!is.null(inner) && nrow(inner) > 0) {
    cd$is_inner <- spot_ids %in% inner$spot_id

    # Add inner metadata
    inner_interface <- rep(NA_character_, length(spot_ids))
    inner_mode <- rep(NA_character_, length(spot_ids))
    inner_directed_pair <- rep(NA_character_, length(spot_ids))

    idx <- match(inner$spot_id, spot_ids)
    valid_idx <- !is.na(idx)

    if ("interface" %in% colnames(inner)) {
      inner_interface[idx[valid_idx]] <- as.character(inner$interface[valid_idx])
    }

    if ("region" %in% colnames(inner)) {
      inner_mode[idx[valid_idx]] <- as.character(inner$region[valid_idx])
    }

    if ("directed_pair" %in% colnames(inner)) {
      inner_directed_pair[idx[valid_idx]] <- as.character(inner$directed_pair[valid_idx])
    }

    cd$is_inner <- cd$is_inner
    cd$inner_interface <- inner_interface
    cd$inner_mode <- inner_mode
    cd$inner_directed_pair <- inner_directed_pair
  }

  # Update SPE with new colData
  colData(spe) <- cd

  spe
}
