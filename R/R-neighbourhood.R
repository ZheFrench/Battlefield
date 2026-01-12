#' Get neighborhood spots around a target cluster or point
#'
#' This function identifies all spots in the neighborhood of a target cluster
#' or a specific point using k-nearest neighbors search. It returns all neighbor
#' spots
#' that were identified, marking which ones are neighbors.
#'
#' @param df A data.frame containing at least the columns `x`, `y`,
#'   cluster identifier column, and `spot_id`.
#'   - `x`, `y`: numeric coordinates of spots
#' - cluster column (name specified by `cluster_col`): cluster assignment for
#' each spot
#'   - `spot_id`: unique identifier for each spot
#' @param cluster The cluster label for which to analyze the neighborhood. Either
#' `cluster`
#'   or `spot_id` must be provided, but not both.
#' @param spot_id Character. The spot identifier to find neighbors for. Either
#' `cluster`
#'   or `spot_id` must be provided, but not both.
#' @param k Integer. Number of nearest neighbors to consider. Default is 100.
#'   Larger values capture more distant neighbors. If k exceeds the number of 
#'   spots, it will be capped at nrow(df) - 1.
#' @param max_dist Numeric. Maximum distance from the target to consider.
#' If `NULL` (default), all neighbors up to k are included without distance
#' filtering.
#' @param inlaid_col Character. Name of the column containing inlaid/annotation
#' values to return in the `neighborhood` column. If `NULL` (default), uses the
#'   `cluster_col` value. This allows counting neighborhood composition by any
#'   column (e.g., cell type, tissue type).
#' @param coords Character vector of length 2 giving the coordinate column
#' names.
#'   Default is `c("x", "y")`.
#' @param cluster_col Character. Name of the column containing cluster
#' assignments.
#'   Default is `"cluster"`.
#'
#' @details
#' The function can work in two modes:
#' - **Cluster mode**: If `cluster` is provided, computes k-nearest neighbors for
#' all spots
#' in that cluster. Returns all neighbors found (excluding source cluster
#' spots).
#' - **Point mode**: If `spot_id` is provided, finds the k-nearest neighbors to
#that
#'   specific spot.
#'
#' @return A data.frame with columns:
#' \describe{
#'   \item{spot_id}{Unique spot identifier.}
#'   \item{x}{X coordinate.}
#'   \item{y}{Y coordinate.}
#'   \item{neighborhood}{Value from `inlaid_col` (or cluster if inlaid_col is NULL) of the neighbor spot.}
#'   \item{cluster}{The source cluster or spot_id being analyzed.}
#'   \item{is_neighborhood}{Logical, always TRUE for returned rows.}
#' }
#' Rows are sorted by neighborhood and spot_id.
#'
#' @examples
#' data("visium_simulated_spe", package = "Battlefield")
#' spe <- visium_simulated_spe
#' df <- data.frame(
#'   spot_id = colnames(spe),
#'   x = spatialCoords(spe)[, 1],
#'   y = spatialCoords(spe)[, 2],
#'   cluster = colData(spe)$cluster
#' )
#' # Get all neighbor spots for cluster 1
#' neighbors_cluster <- get_neighborhood_spots(df, cluster = 1, k = 50)
#' head(neighbors_cluster)
#'
#' # Get neighbors for a specific point
#' first_spot <- df$spot_id[1]
#' neighbors_point <- get_neighborhood_spots(df, spot_id = first_spot, k = 10)
#' head(neighbors_point)
#'
#' @importFrom RANN nn2
#' @importFrom dplyr arrange
#' @export
get_neighborhood_spots <- function(df,
                                cluster = NULL,
                                spot_id = NULL,
                                k = 100,
                                max_dist = NULL,
                                inlaid_col = NULL,
                                coords = c("x", "y"),
                                cluster_col = "cluster") {

# Validation
stopifnot(is.data.frame(df))
stopifnot(all(coords %in% colnames(df)))
stopifnot(cluster_col %in% colnames(df))
stopifnot("spot_id" %in% colnames(df))
stopifnot(nrow(df) >= 1)

# If inlaid_col is not provided, use cluster_col
if (is.null(inlaid_col)) {
    inlaid_col <- cluster_col
} else {
    stopifnot(inlaid_col %in% colnames(df))
}

# Check that exactly one of cluster or spot_id is provided
if ((is.null(cluster) && is.null(spot_id)) || (!is.null(cluster) && !is.null(spot_id))) {    stop("Exactly one of 'cluster' or 'spot_id' must be provided,
    not both or neither.")}

# Get coordinates
X <- as.matrix(df[, coords])

# Determine query mode: cluster or point
if (!is.null(cluster)) {
    # Cluster mode: find neighbors for all spots in the cluster
    target_idx <- which(df[[cluster_col]] == cluster)
    if (length(target_idx) == 0) {
    stop("No spots found for cluster = ", cluster)
    }
    source_label <- cluster
    exclude_cluster <- cluster
} else {
    # Point mode: find neighbors for a specific spot
    target_idx <- which(df$spot_id == spot_id)
    if (length(target_idx) == 0) {
    stop("No spot found with spot_id = ", spot_id)
    }
    source_label <- spot_id
    exclude_cluster <- NULL  # Don't exclude any cluster in point mode
}

# Find neighbors for target spot(s)
target_coords <- X[target_idx, , drop = FALSE]

# Compute k-nearest neighbors
# Guardrail: ensure k_request doesn't exceed available points (RANN requires k <
#  nrow(X))
# We request more neighbors to account for points that will be filtered
max_available <- nrow(df) - 1
k_request <- min(k + length(target_idx), max_available)
if (k_request < 1) {
    stop("Insufficient points in dataset to compute neighbors (need at least 2 points).")}
nn <- RANN::nn2(X, target_coords, k = k_request)

# Collect neighbor indices and distances
# Note: RANN::nn2 returns nn.idx (not nn.index) and nn.dists
neighbor_indices <- as.vector(nn$nn.idx)
neighbor_dists <- as.vector(nn$nn.dists)

# Filter by max_dist if specified
if (!is.null(max_dist)) {
    valid <- neighbor_dists <= max_dist
    neighbor_indices <- neighbor_indices[valid]
}

# Remove the query spots from neighbors (self)
neighbor_indices <- neighbor_indices[!neighbor_indices %in% target_idx]

# In cluster mode, also exclude spots from the source cluster
if (!is.null(exclude_cluster)) {
    neighbor_indices <- neighbor_indices[df[[cluster_col]][neighbor_indices] != exclude_cluster]}

if (length(neighbor_indices) == 0) {
    message("No neighbors found for cluster = ", source_label,
    " within specified parameters.")    
    return(data.frame(
    spot_id = character(),
    x = numeric(),
    y = numeric(),
    cluster = character(),
    neighborhood = character(),
    is_neighborhood = logical()
    ))
}

# Keep only the top k results
neighbor_indices <- head(neighbor_indices, k * length(target_idx))

# Remove duplicates: keep only unique neighbor spots
neighbor_indices <- unique(neighbor_indices)

# Limit to k unique neighbors
neighbor_indices <- head(neighbor_indices, k)

# Build result dataframe
result <- df[neighbor_indices, c("spot_id", coords, inlaid_col), drop = FALSE]
colnames(result)[colnames(result) == inlaid_col] <- "neighborhood"
result$cluster <- source_label
result$is_neighborhood <- TRUE
result <- result[, c("spot_id", coords, "cluster", "neighborhood", 
"is_neighborhood")]
# Sort by neighborhood and spot_id
result <- dplyr::arrange(result, neighborhood, spot_id)
result <- as.data.frame(result)

# Set rownames to spot_id for uniqueness
rownames(result) <- result$spot_id

result
}



#' Count annotated spot types in the neighborhood of a cluster or point
#'
#' This function computes summary statistics of cluster types in the
#' neighborhood
#' of a target cluster or a specific point. It uses [get_neighborhood_spots()]
#' internally
#' to identify neighbors, then counts how many belong to each cluster type.
#'
#' @param df A data.frame containing at least the columns `x`, `y`, cluster
#' identifier column, and `spot_id`.
#'   - `x`, `y`: numeric coordinates of spots
#' - cluster column (name specified by `cluster_col`): cluster assignment for
#' each spot
#' - `spot_id`: unique identifier for each spot (optional, required for detailed
#' analysis)
#' @param cluster The cluster label for which to analyze the neighborhood. Either
#' `cluster`
#'   or `spot_id` must be provided, but not both.
#' @param spot_id Character. The spot identifier to find neighbors for. Either
#' `cluster`
#'   or `spot_id` must be provided, but not both.
#' @param k Integer. Number of nearest neighbors to consider. Default is 100.
#'   Larger values capture more distant neighbors.
#' @param max_dist Numeric. Maximum distance from the target to consider.
#' If `NULL` (default), all neighbors up to k are included without distance
#' #' filtering.
#' @param inlaid_col Character. Name of the column containing inlaid/annotation
#' values to count in neighborhoods. If `NULL` (default), uses the
#'   `cluster_col` value. This allows counting neighborhood composition by any
#'   column (e.g., cell type, tissue type).
#' @param coords Character vector of length 2 giving the coordinate column
#' names.
#'   Default is `c("x", "y")`.
#' @param cluster_col Character. Name of the column containing cluster
#' assignments.
#'   Default is `"cluster"`.
#'
#' @details
#' This function wraps [get_neighborhood_spots()] and summarizes the results by
#' counting spots from each neighboring cluster or inlaid type. In cluster mode, the source
#' cluster
#' itself is excluded from the counts. In point mode, all neighboring clusters
#' are counted.
#'
#' The result shows the composition of the neighborhood: how many spots of each
#' cluster
#' type surround the target.
#'
#' @return A data.frame with columns:
#' \describe{
#'   \item{neighborhood}{Cluster label or inlaid value (from neighboring spots).}
#'   \item{count}{Number of times this cluster/type appears in the neighborhood.}
#'   \item{proportion}{Proportion of this cluster/type relative to all neighbors.}
#' }
#' Rows are sorted by count in descending order.
#'
#' @examples
#' data("visium_simulated_spe", package = "Battlefield")
#' spe <- visium_simulated_spe
#' df <- data.frame(
#'   spot_id = colnames(spe),
#'   x = spatialCoords(spe)[, 1],
#'   y = spatialCoords(spe)[, 2],
#'   cluster = colData(spe)$cluster
#' )
#' # Count cluster types in neighborhood of cluster 1
#' neighbor_counts <- count_neighborhood(df, cluster = 1, k = 50)
#' neighbor_counts
#'
#' # Count cluster types around a specific point
#' first_spot <- df$spot_id[1]
#' point_neighbors <- count_neighborhood(df, spot_id = first_spot, k = 10)
#' point_neighbors
#'
#' @importFrom RANN nn2
#' @importFrom dplyr group_by summarise arrange desc
#' @export
count_neighborhood <- function(df,
                            cluster = NULL,
                            spot_id = NULL,
                            k = 100,
                            max_dist = NULL,
                            inlaid_col = NULL,
                            coords = c("x", "y"),
                            cluster_col = "cluster") {

# Use get_neighborhood_spots to identify neighbors
neighbors <- get_neighborhood_spots(df, cluster, spot_id, k, max_dist, inlaid_col, coords,
cluster_col)
if (nrow(neighbors) == 0) {
    return(data.frame(cluster = character(), neighborhood = character(), count = integer(),
    proportion = numeric()))}

# Create summary
result <- data.frame(neighborhood = neighbors$neighborhood) |>
    dplyr::group_by(neighborhood) |>
    dplyr::summarise(count = dplyr::n(), .groups = "drop") |>
    dplyr::mutate(proportion = count / sum(count)) |>
    dplyr::arrange(dplyr::desc(count)) |>
    as.data.frame()

# Add cluster column to track source
source_label <- if (!is.null(cluster)) cluster else spot_id
result$cluster <- source_label
result <- result[, c("cluster", "neighborhood", "count", "proportion")]

result
}


#' Count neighborhood composition for all clusters
#'
#' This function computes neighborhood statistics for all clusters in a dataset
#' in a single call. It returns a dataframe with counts of neighboring cluster
#' types
#' for each source cluster.
#'
#' @param df A data.frame containing at least the columns `x`, `y`, cluster
#' identifier column, and `spot_id`.
#'   - `x`, `y`: numeric coordinates of spots
#' - cluster column (name specified by `cluster_col`): cluster assignment for
#' each spot
#' - `spot_id`: unique identifier for each spot (optional, required for detailed
#' analysis)
#' @param clusters Optional vector of cluster labels to process. If `NULL`,
#'   all unique clusters in `df` are used.
#' @param k Integer. Number of nearest neighbors to consider. Default is 100.
#'   Larger values capture more distant neighbors.
#' @param max_dist Numeric. Maximum distance from the target to consider.
#' If `NULL` (default), all neighbors up to k are included without distance
#' filtering.
#' @param inlaid_col Character. Name of the column containing inlaid/annotation
#' values to count in neighborhoods. If `NULL` (default), uses the
#'   `cluster_col` value. This allows counting neighborhood composition by any
#'   column (e.g., cell type, tissue type).
#' @param coords Character vector of length 2 giving the coordinate column
#' names.
#'   Default is `c("x", "y")`.
#' @param cluster_col Character. Name of the column containing cluster
#' assignments.
#'   Default is `"cluster"`.
#'
#' @details
#' This is a batch processing function that applies [count_neighborhood()] to
#' all
#' clusters and combines the results into a single dataframe. Each row
#' represents
#' a neighbor cluster found around a source cluster, with counts and
#' proportions.
#'
#' @return A data.frame with columns:
#' \describe{
#'   \item{cluster}{The source cluster being analyzed.}
#'   \item{neighborhood}{The neighbor cluster type or inlaid value.}
#'   \item{count}{Number of neighbor spots of this cluster type.}
#' \item{proportion}{Proportion of this cluster among all neighbors of the
#' cluster.}
#' }
#' Rows are grouped by source cluster and sorted by count within each group
#' (descending).
#'
#' @examples
#' data("visium_simulated_spe", package = "Battlefield")
#' spe <- visium_simulated_spe
#' df <- data.frame(
#'   spot_id = colnames(spe),
#'   x = spatialCoords(spe)[, 1],
#'   y = spatialCoords(spe)[, 2],
#'   cluster = colData(spe)$cluster
#' )
#' # Get neighborhood statistics for all clusters
#' all_neighbors <- count_all_neighborhoods(df, k = 50)
#' head(all_neighbors)
#'
#' @importFrom dplyr bind_rows
#' @export
count_all_neighborhoods <- function(df,
                                    clusters = NULL,
                                    k = 100,
                                    max_dist = NULL,
                                    inlaid_col = NULL,
                                    coords = c("x", "y"),
                                    cluster_col = "cluster") {

# Validation
stopifnot(is.data.frame(df))
stopifnot(all(coords %in% colnames(df)))
stopifnot(cluster_col %in% colnames(df))
stopifnot(nrow(df) >= 1)

# Determine clusters to process
if (is.null(clusters)) {
    clusters <- unique(df[[cluster_col]])
    clusters <- clusters[!is.na(clusters)]
}

# Compute neighborhood stats for each cluster
results <- lapply(clusters, function(clust) {
    neighbor_counts <- count_neighborhood(df, cluster = clust, k = k,
    max_dist = max_dist, inlaid_col = inlaid_col, coords = coords, cluster_col = cluster_col)
    # Add cluster column
    neighbor_counts$cluster <- clust
    neighbor_counts <- neighbor_counts[, c("cluster", "neighborhood", "count",
    "proportion")]
    neighbor_counts
})

# Combine results
combined <- dplyr::bind_rows(results)

# Sort by cluster and count (descending)
combined <- combined |>
    dplyr::arrange(cluster, dplyr::desc(count)) |>
    as.data.frame()

combined
}


#' Get inlaid spots within a source cluster
#'
#' This function retrieves all spots from a source cluster along with their
#' inlaid/annotation values. Unlike [get_neighborhood_spots()], this returns
#' spots
#' *inside* the source cluster itself, not around it.
#'
#' @param df A data.frame containing at least the columns `x`, `y`, cluster
#' identifier column, inlaid column, and `spot_id`.
#'   - `x`, `y`: numeric coordinates of spots
#' - cluster column (name specified by `cluster_col`): cluster assignment for
#' each spot
#' - inlaid column (name specified by `inlaid_col`): inlaid/annotation for each
#' spot
#'   - `spot_id`: unique identifier for each spot
#' @param cluster The cluster label to retrieve inlaid spots for.
#' @param inlaid_col Character. Name of the column containing inlaid/annotation
#' values.
#'   Default is `"cluster"`.
#' @param cluster_col Character. Name of the column containing cluster
#' assignments.
#'   Default is `"cluster"`.
#' @param coords Character vector of length 2 giving the coordinate column
#' names.
#'   Default is `c("x", "y")`.
#'
#' @return A data.frame with columns:
#' \describe{
#'   \item{spot_id}{Unique spot identifier.}
#'   \item{x}{X coordinate.}
#'   \item{y}{Y coordinate.}
#'   \item{inlaid}{Inlaid/annotation value.}
#'   \item{cluster}{The source cluster being analyzed.}
#'   \item{is_inlaid}{Logical, always TRUE for returned rows.}
#' }
#' Rows are sorted by inlaid and spot_id.
#'
#' @examples
#' data("visium_simulated_spe", package = "Battlefield")
#' spe <- visium_simulated_spe
#' df <- data.frame(
#'   spot_id = colnames(spe),
#'   x = spatialCoords(spe)[, 1],
#'   y = spatialCoords(spe)[, 2],
#'   cluster = colData(spe)$cluster,
#' inlaid = sample(paste0("type_", 1:3), length(colnames(spe)), replace = TRUE)
#' )
#' # Get all inlaid spots within cluster 1
#' inlaid_spots <- get_inlaid_spots(df, cluster = 1, inlaid_col = "inlaid")
#' head(inlaid_spots)
#'
#' @importFrom dplyr arrange
#' @export
get_inlaid_spots <- function(df,
                            cluster,
                            inlaid_col = "cluster",
                            cluster_col = "cluster",
                            coords = c("x", "y")) {

# Validation
stopifnot(is.data.frame(df))
stopifnot(all(coords %in% colnames(df)))
stopifnot(cluster_col %in% colnames(df))
stopifnot(inlaid_col %in% colnames(df))
stopifnot("spot_id" %in% colnames(df))
stopifnot(nrow(df) >= 1)

# Get spots from cluster
source_idx <- which(df[[cluster_col]] == cluster)

if (length(source_idx) == 0) {
    stop("No spots found for cluster = ", cluster)
}

# Build result dataframe
result <- df[source_idx, c("spot_id", coords, inlaid_col), drop = FALSE]
colnames(result)[colnames(result) == inlaid_col] <- "inlaid"
result$cluster <- cluster
result$is_inlaid <- TRUE
result <- result[, c("spot_id", coords, "cluster", "inlaid", "is_inlaid")]

# Sort by inlaid and spot_id
result <- dplyr::arrange(result, inlaid, spot_id)
result <- as.data.frame(result)

# Set rownames to spot_id for uniqueness
rownames(result) <- result$spot_id

result
}


#' Count inlaid composition within a source cluster
#'
#' This function counts the composition of an inlaid/annotation column within
#' a target source cluster. Unlike [count_neighborhood()], this counts types
#' *inside* the source cluster itself, not around it.
#'
#' @param df A data.frame containing at least the columns `x`, `y`, cluster
#' identifier column, and `spot_id`.
#'   - `x`, `y`: numeric coordinates of spots
#' - cluster column (name specified by `cluster_col`): cluster assignment for
#' each spot
#' - inlaid column (name specified by `inlaid_col`): inlaid/annotation for each
#' spot
#'   - `spot_id`: unique identifier for each spot
#' @param cluster The cluster label to analyze inlaid composition for.
#' @param inlaid_col Character. Name of the column containing inlaid/annotation
#' values.
#' Default is `"cluster"` (to analyze cluster composition within another
#' grouping).
#' @param cluster_col Character. Name of the column containing cluster
#' assignments.
#'   Default is `"cluster"`.
#'
#' @return A data.frame with columns:
#' \describe{
#'   \item{inlaid}{Inlaid/annotation value.}
#'   \item{count}{Number of spots with this inlaid type in the cluster.}
#' \item{proportion}{Proportion of this inlaid type among all spots in cluster.}
#' }
#' Rows are sorted by count in descending order.
#'
#' @examples
#' data("visium_simulated_spe", package = "Battlefield")
#' spe <- visium_simulated_spe
#' df <- data.frame(
#'   spot_id = colnames(spe),
#'   x = spatialCoords(spe)[, 1],
#'   y = spatialCoords(spe)[, 2],
#'   cluster = colData(spe)$cluster,
#' inlaid = sample(paste0("type_", 1:3), length(colnames(spe)), replace = TRUE)
#' )
#' # Count inlaid types within cluster 1
#' inlaid_counts <- count_inlaid(df, cluster = 1, inlaid_col = "inlaid")
#' inlaid_counts
#'
#' @importFrom dplyr group_by summarise arrange desc
#' @export
count_inlaid <- function(df,
                        cluster,
                        inlaid_col = "cluster",
                        cluster_col = "cluster") {

# Validation
stopifnot(is.data.frame(df))
stopifnot(cluster_col %in% colnames(df))
stopifnot(inlaid_col %in% colnames(df))
stopifnot(nrow(df) >= 1)

# Use get_inlaid_spots to retrieve inlaid spots
inlaid_spots <- get_inlaid_spots(df, cluster, inlaid_col, cluster_col)

if (nrow(inlaid_spots) == 0) {
    return(data.frame(cluster = character(), inlaid = character(), count = integer(),
    proportion = numeric()))}

# Count inlaid composition
result <- data.frame(inlaid = inlaid_spots$inlaid) |>
    dplyr::group_by(inlaid) |>
    dplyr::summarise(count = dplyr::n(), .groups = "drop") |>
    dplyr::mutate(proportion = count / sum(count)) |>
    dplyr::arrange(dplyr::desc(count)) |>
    as.data.frame()

# Add cluster column to track source
result$cluster <- cluster
result <- result[, c("cluster", "inlaid", "count", "proportion")]

result
}


#' Count inlaid composition for all clusters
#'
#' This function counts inlaid/annotation composition within each cluster in a
#' single call. It returns a dataframe with counts of inlaid types for each
#' source cluster.
#'
#' @param df A data.frame containing at least the columns `x`, `y`, cluster
#' identifier column, and `spot_id`.
#'   - `x`, `y`: numeric coordinates of spots
#' - cluster column (name specified by `cluster_col`): cluster assignment for
#' #' each spot
#' - inlaid column (name specified by `inlaid_col`): inlaid/annotation for each
#' spot
#'   - `spot_id`: unique identifier for each spot
#' @param clusters Optional vector of cluster labels to process. If `NULL`,
#'   all unique clusters in `df` are used.
#' @param inlaid_col Character. Name of the column containing inlaid/annotation
#' values.
#'   Default is `"cluster"`.
#' @param cluster_col Character. Name of the column containing cluster
#' assignments.
#'   Default is `"cluster"`.
#'
#' @return A data.frame with columns:
#' \describe{
#'   \item{cluster}{The source cluster being analyzed.}
#'   \item{inlaid}{The inlaid/annotation type.}
#'   \item{count}{Number of spots with this inlaid type in the cluster.}
#' \item{proportion}{Proportion of this inlaid type among all spots in cluster.}
#' }
#' Rows are grouped by source cluster and sorted by count within each group
#' (descending).
#'
#' @examples
#' data("visium_simulated_spe", package = "Battlefield")
#' spe <- visium_simulated_spe
#' df <- data.frame(
#'   spot_id = colnames(spe),
#'   x = spatialCoords(spe)[, 1],
#'   y = spatialCoords(spe)[, 2],
#'   cluster = colData(spe)$cluster,
#' inlaid = sample(paste0("type_", 1:3), length(colnames(spe)), replace = TRUE)
#' )
#' # Get inlaid statistics for all clusters
#' all_inlaids <- count_all_inlaids(df, inlaid_col = "inlaid")
#' head(all_inlaids)
#'
#' @importFrom dplyr bind_rows
#' @export
count_all_inlaids <- function(df,
                            clusters = NULL,
                            inlaid_col = "cluster",
                            cluster_col = "cluster") {

# Validation
stopifnot(is.data.frame(df))
stopifnot(cluster_col %in% colnames(df))
stopifnot(inlaid_col %in% colnames(df))
stopifnot(nrow(df) >= 1)

# Determine clusters to process
if (is.null(clusters)) {
    clusters <- unique(df[[cluster_col]])
    clusters <- clusters[!is.na(clusters)]
}

# Compute inlaid stats for each cluster
results <- lapply(clusters, function(clust) {
    inlaid_counts <- count_inlaid(df, cluster = clust, inlaid_col = inlaid_col,
    cluster_col = cluster_col)
    # Add cluster column
    inlaid_counts$cluster <- clust
    inlaid_counts <- inlaid_counts[, c("cluster", "inlaid", "count",
    "proportion")]
    inlaid_counts
})

# Combine results
combined <- dplyr::bind_rows(results)

# Sort by cluster and count (descending)
combined <- combined |>
    dplyr::arrange(cluster, dplyr::desc(count)) |>
    as.data.frame()

combined
}
