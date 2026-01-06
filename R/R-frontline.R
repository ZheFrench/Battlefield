#' Build all oriented cluster pairs (A -> B) and their directed_pair labels
#'
#' Given a vector of cluster labels, this function returns all **ordered**
#' (oriented) pairs of distinct clusters. For each pair (`cluster`,
#' `interface`), it also creates a `directed_pair` string label such as `"A-B"`.
#'
#' @param cluster_labels A vector of cluster labels (character, factor, numeric, etc.).
#'   `NA` values are ignored.
#' @param interface_separator Character string used to join `cluster` and
#'   `interface` into the `directed_pair` label. Default is `"-"`.
#' @param sort_clusters Logical. If `TRUE`, unique cluster labels are sorted
#'   (lexicographically) before building pairs. If `FALSE`, preserves first
#'   appearance order. Default is `FALSE`.
#'
#' @return A data.frame with columns:
#' \describe{
#'   \item{cluster}{Source cluster label (character).}
#'   \item{interface}{Target cluster label (character).}
#'   \item{directed_pair}{Directed pair label, e.g. `"A-B"`.
#' }
#' If fewer than 2 unique (non-NA) clusters are present, returns an empty data.frame
#' with the same columns.
#'
#' @examples
#' cl <- c("1","1","2","3", NA)
#' directed_cluster_interface_pairs(cl)
#'
#' @export
directed_cluster_interface_pairs <- function(cluster_labels,
                                             interface_separator = "-",
                                             sort_clusters = FALSE) {

  cl <- unique(as.character(cluster_labels))
  cl <- cl[!is.na(cl)]

  if (sort_clusters) {
    cl <- sort(cl)
  }

  if (length(cl) < 2) {
    return(data.frame(
      cluster = character(0),
      interface   = character(0),
      directed_pair    = character(0),
      stringsAsFactors = FALSE
    ))
  }

  pairs <- expand.grid(
    cluster = cl,
    interface   = cl,
    stringsAsFactors = FALSE
  )

  pairs <- subset(pairs, cluster != interface)

  pairs$directed_pair <- paste0(pairs$cluster, interface_separator, pairs$interface)
  
  rownames(pairs) <- NULL

  pairs
}

#' Select border spots from cluster A that touch cluster B (and flag junctions)
#'
#' This function identifies **directed interface spots**: spots belonging to
#' `cluster` (A) that have at least one neighbor in `interface` (B) within
#' a k-nearest-neighbors neighborhood (optionally constrained by a distance cutoff).
#'
#' In addition, it flags spots that also touch **other clusters** (i.e., junction
#' / multi-interface spots) and reports which other clusters are touched.
#'
#' @param df A data.frame containing at least the coordinate columns and a cluster label column.
#' @param cluster Cluster label for the source cluster (A). Spots in this cluster
#'   are tested for adjacency to `interface`.
#' @param interface Cluster label for the target cluster (B).
#' @param k Integer. Number of nearest neighbors to consider (excluding self).
#'   Internally uses `k + 1` to include self then removes it. Default is 7.
#' @param max_dist Optional numeric. Maximum Euclidean distance for a neighbor to be
#'   considered (distance cutoff). If `NULL`, no distance filtering is applied.
#' @param coord_cols Character vector of length 2 giving the coordinate column names.
#'   Default is `c("x","y")`.
#' @param cluster_col Character. Name of the column containing cluster labels.
#'   Default is `"cluster"`.
#'
#' @details
#' Steps:
#' \enumerate{
#'   \item Builds a kNN neighborhood for each spot in `cluster`.
#'   \item Optionally removes neighbors farther than `max_dist`.
#'   \item Marks a spot as a border spot if it has at least one neighbor in `interface`.
#'   \item Flags a spot as multi-interface if it also has a neighbor belonging to a
#'   cluster different from `cluster` and `interface`.
#' }
#'
#' The returned data.frame contains only spots from `cluster` that touch `interface`,
#' with columns:
#' \describe{
#'   \item{spot_id}{Spot identifier.}
#'   \item{x}{X coordinate.}
#'   \item{y}{Y coordinate.}
#'   \item{cluster}{Cluster label (same as `cluster`).}
#'   \item{interface}{The value of `interface`.}
#'   \item{directed_pair}{Interface label combining cluster and `interface`, e.g. `"A-B"`.
#'   \item{is_border}{Always `TRUE` for returned rows (kept for clarity).}
#'   \item{is_border_multiple}{`TRUE` if the spot also touches other clusters.}
#'   \item{other_adjacent_borders}{Comma-separated list of other clusters touched, or `NA`.}
#' }
#'
#' @return A data.frame subset of `df` containing border spots from `cluster` to
#' `interface` with columns: spot_id, x, y, cluster, interface, directed_pair,
#' is_border, is_border_multiple, other_adjacent_borders.
#'
#' @examples
#' # Minimal example (requires RANN and dplyr)
#' set.seed(1)
#' df_ex <- data.frame(
#'   x = rnorm(200),
#'   y = rnorm(200),
#'   cluster = sample(c("A","B","C"), 200, replace = TRUE)
#' )
#' res <- select_border_spots(df_ex, cluster = "A", interface = "B", k = 7)
#' head(res)
#'
#' @importFrom RANN nn2
#' @export
select_border_spots <- function(df,
                                cluster,
                                interface,
                                k = 7,
                                max_dist = NULL,
                                coord_cols = c("x", "y"),
                                cluster_col = "cluster") {
  stopifnot(all(coord_cols %in% colnames(df)))
  stopifnot(cluster_col %in% colnames(df))

  coords <- as.matrix(df[, coord_cols])
  cl <- df[[cluster_col]]

  idx_from <- which(cl == cluster)
  if (length(idx_from) == 0) {
    return(df[0, , drop = FALSE])
  }

  nn <- RANN::nn2(data = coords, query = coords[idx_from, , drop = FALSE], k = k + 1)

  nn_idx  <- nn$nn.idx[, -1, drop = FALSE]     # remove self
  nn_dist <- nn$nn.dists[, -1, drop = FALSE]

  neigh_cl <- matrix(cl[nn_idx], nrow = nrow(nn_idx))

  # Optional distance filter
  if (!is.null(max_dist)) {
    neigh_cl[nn_dist > max_dist] <- NA
  }

  # Does a spot in A touch B?
  touches_to <- apply(neigh_cl == interface, 1, any, na.rm = TRUE)

  # Does it touch any other cluster (not A, not B)?
  touches_other <- apply(!(neigh_cl %in% c(cluster, interface)) & !is.na(neigh_cl),
                         1, any)

  # Which other clusters are touched?
  other_clusters <- apply(neigh_cl, 1, function(v) {
    v <- v[!is.na(v)]
    v <- unique(v[!(v %in% c(cluster, interface))])
    if (length(v) == 0) NA_character_ else paste(sort(v), collapse = ",")
  })

  # Keep all A->B border points
  out <- df[idx_from[touches_to], , drop = FALSE]

  # Add annotations and reorder columns
  out <- out |>
    dplyr::mutate(
      interface = interface,
      directed_pair = paste0(cluster, "-", interface),
      is_border = TRUE,
      is_border_multiple = touches_other[touches_to],
      other_adjacent_borders = other_clusters[touches_to]
    ) |>
    dplyr::select(
      spot_id, x, y, directed_pair, cluster, interface,
      is_border, is_border_multiple, other_adjacent_borders
    )

  out
}

#' Build border spots for all oriented cluster pairs
#'
#' This function iterates over all **ordered** cluster pairs (A -> B) and returns
#' a single data.frame containing the border spots for each interface, as computed
#' by a border-detection function (e.g. [select_border_spots()]).
#'
#' If `pairs_df` is not provided, it is generated from the cluster labels using
#' [directed_cluster_interface_pairs()].
#'
#' @param df A data.frame containing at least coordinate columns and a cluster label column.
#' @param k Integer. Number of nearest neighbors to consider (excluding self) when
#'   computing borders. Default is 6.
#' @param max_dist Optional numeric. Maximum Euclidean distance for neighbors to be
#'   considered. If `NULL`, no distance filtering is applied.
#' @param pairs Optional data.frame of oriented cluster pairs, typically produced by
#'   [directed_cluster_interface_pairs()]. Must contain `cluster` and `interface`
#'   columns. If `NULL`, it is computed from `df[[cluster_col]]`.
#' @param coord_cols Character vector of length 2 giving the coordinate column names.
#'   Default is `c("x","y")`.
#' @param cluster_col Character. Name of the column containing cluster labels.
#'   Default is `"cluster"`.
#'
#' @return A data.frame produced by row-binding the result of border selection for
#' each oriented pair. Typically contains the original columns of `df` plus interface
#' annotation columns from the border selector.
#'
#' @examples
#' # Example with synthetic data
#' set.seed(1)
#' df_ex <- data.frame(
#'   x = rnorm(200),
#'   y = rnorm(200),
#'   cluster = sample(c("A","B","C"), 200, replace = TRUE)
#' )
#' all_borders <- build_all_borders(df_ex, k = 6)
#' head(all_borders)
#' @importFrom dplyr bind_rows
#' @export
build_all_borders <- function(df,
                              k = 6,
                              max_dist = NULL,
                              pairs = NULL,
                              coord_cols = c("x", "y"),
                              cluster_col = "cluster") {

  if (is.null(pairs)) {
    pairs <- directed_cluster_interface_pairs(df[[cluster_col]])
  }

  if (nrow(pairs) == 0) {
    return(df[0, , drop = FALSE])
  }

  res <- lapply(seq_len(nrow(pairs)), function(i) {
    a <- pairs$cluster[i]
    b <- pairs$interface[i]

    select_border_spots(
      df,
      cluster = a,
      interface = b,
      k = k,
      max_dist = max_dist,
      coord_cols = coord_cols,
      cluster_col = cluster_col
    )

  })

  out <- dplyr::bind_rows(res)
  
  # Create undirected pair label showing both directions (e.g., "3-4 / 4-3")
  out <- out |>
    dplyr::rowwise() |>
    dplyr::mutate(
      undirected_pair = paste0(
        min(as.character(cluster), as.character(interface)),
        "-",
        max(as.character(cluster), as.character(interface)),
        " / ",
        max(as.character(cluster), as.character(interface)),
        "-",
        min(as.character(cluster), as.character(interface))
      )
    ) |>
    dplyr::ungroup() |> as.data.frame()
  
  out
}

#' Build inner spots for all oriented cluster pairs
#'
#' This function iterates over all **ordered** cluster pairs (A -> B) and returns
#' a single data.frame containing the inner (control) spots for each pair, as computed
#' by [select_inner_spots()].
#'
#' If `pairs` is not provided, it is generated from the cluster labels using
#' [directed_cluster_interface_pairs()].
#'
#' @param df A data.frame containing at least coordinate columns and a cluster label column.
#' @param border_df A data.frame of border spots from [build_all_borders()].
#' @param region Character. One of "inner", "outer", or "both". See [select_inner_spots()] for details.
#'   Default is "inner".
#' @param pairs Optional data.frame of oriented cluster pairs, typically produced by
#'   [directed_cluster_interface_pairs()]. Must contain `cluster` and `interface`
#'   columns. If `NULL`, it is computed from `df[[cluster_col]]`.
#' @param coord_cols Character vector of length 2 giving the coordinate column names.
#'   Default is `c("x","y")`.
#' @param cluster_col Character. Name of the column containing cluster labels.
#'   Default is `"cluster"`.
#'
#' @return A data.frame produced by row-binding the result of inner spot selection for
#' each oriented pair. Typically contains the original columns of `df` plus annotation
#' columns (is_inner, interface, region) from the inner spot selector.
#'
#' @examples
#' # Example with synthetic data
#' set.seed(1)
#' df_ex <- data.frame(
#'   x = rnorm(200),
#'   y = rnorm(200),
#'   cluster = sample(c("A","B","C"), 200, replace = TRUE)
#' )
#' all_borders <- build_all_borders(df_ex, k = 6)
#' all_inners <- build_all_inners(df_ex, all_borders, region = "inner")
#' head(all_inners)
#'
#' @importFrom dplyr bind_rows
#' @export
build_all_inners <- function(df,
                             border_df,
                             region = "inner",
                             pairs = NULL,
                             coord_cols = c("x", "y"),
                             cluster_col = "cluster") {

  # If pairs not provided, extract unique cluster/interface pairs from border_df
  if (is.null(pairs)) {
    pairs <- unique(border_df[, c("cluster", "interface")])
  } else {
    # Validate that provided pairs exist in border_df
    border_pairs <- unique(border_df[, c("cluster", "interface")])
    
    # Check which provided pairs are in border_pairs
    valid_pairs <- rep(FALSE, nrow(pairs))
    for (i in seq_len(nrow(pairs))) {
      valid_pairs[i] <- any(
        border_pairs$cluster == pairs$cluster[i] & border_pairs$interface == pairs$interface[i]
      )
    }
    
    if (!all(valid_pairs)) {
      invalid_count <- sum(!valid_pairs)
      warn <- paste0("Warning: ", invalid_count, " provided pair(s) not found in border_df.
       Filtering them out.")
      warning(warn)
      pairs <- pairs[valid_pairs, ]
    }
  }

  if (nrow(pairs) == 0) {
    return(df[0, , drop = FALSE])
  }

  res <- lapply(seq_len(nrow(pairs)), function(i) {
    a <- pairs$cluster[i]
    b <- pairs$interface[i]

    select_inner_spots(
      df,
      border_df = border_df,
      cluster = a,
      interface = b,
      region = region,
      coord_cols = coord_cols,
      cluster_col = cluster_col
    )
    
  })

  dplyr::bind_rows(res)
}

#' Select inner (non-interface) spots for a directed pair
#'
#' This function selects inner (non-interface) spots from a cluster that match
#' the count of border spots for a specific directed pair. The `region` parameter
#' controls which direction(s) of the interface to consider when counting borders.
#'
#' For example with `cluster="1"` and `interface="2"`:
#' - If `region="inner"`: counts only border spots from 1→2
#' - If `region="outer"`: counts only border spots from 2→1
#' - If `region="both"`: counts border spots from both 1→2 AND 2→1
#'
#' This returns N random inner spots from the cluster, where N equals the border count.
#'
#' @param df A data.frame containing at least the coordinate columns and a cluster label column.
#' @param border_df A data.frame of border spots from [build_all_borders()], containing
#'   columns `spot_id`, `cluster`, `interface`, `directed_pair`, etc.
#' @param cluster Character. Label of the cluster from which to select inner spots.
#' @param interface Character. Label of the target cluster for the directed pair.
#' @param region Character. One of "inner", "outer", or "both".
#'   - "inner": counts border spots only from cluster → interface
#'   - "outer": counts border spots only from interface → cluster
#'   - "both": counts border spots from both directions
#'   Default is "inner".
#' @param coord_cols Character vector of length 2 giving the coordinate column names.
#'   Default is `c("x","y")`.
#' @param cluster_col Character. Name of the column containing cluster labels.
#'   Default is `"cluster"`.
#'
#' @details
#' Steps:
#' \enumerate{
#'   \item Counts directed border spots based on `region` parameter.
#'   \item Gets all spots in `cluster` that are not at any border.
#'   \item Randomly samples the same number of inner spots as the border count.
#'   \item Returns these sampled inner spots with `inner_cluster` annotation.
#' }
#'
#' If not enough inner spots exist to match the border count, all available inner
#' spots are returned with a warning.
#'
#' @return A data.frame of sampled inner spots from `cluster`, with all columns from `df`
#'   plus `interface` annotation.
#'
#' @examples
#' # Assuming df and big_border_df from build_all_borders(df, k=4)
#' # Inner direction (1→2):
#' # inners_inner <- select_inner_spots(df, big_border_df, "1", "2", region="inner")
#' # Outer direction (2→1):
#' # inners_outer <- select_inner_spots(df, big_border_df, "1", "2", region="outer")
#' # Both directions (1→2 AND 2→1):
#' # inners_both <- select_inner_spots(df, big_border_df, "1", "2", region="both")
#'
#' @importFrom dplyr filter mutate slice_sample pull
#' @export
select_inner_spots <- function(df,
                               border_df,
                               cluster,
                               interface,
                               region = "inner",
                               coord_cols = c("x", "y"),
                               cluster_col = "cluster") {

  stopifnot(all(coord_cols %in% colnames(df)))
  stopifnot(cluster_col %in% colnames(df))
  stopifnot("spot_id" %in% colnames(df))
  stopifnot("spot_id" %in% colnames(border_df))
  stopifnot("cluster" %in% colnames(border_df))
  stopifnot("interface" %in% colnames(border_df))
  stopifnot(region %in% c("inner", "outer", "both"))

  # Check if the cluster and interface pair exists in border_df
  pair_exists <- any(
    (border_df$cluster == cluster & border_df$interface == interface) |
    (border_df$cluster == interface & border_df$interface == cluster)
  )
  
  if (!pair_exists) {
    warn <- paste0("No border spots found for pair (", cluster, ", ", interface, "). Returning NULL.")
    warning(warn)
    return(NULL)
  }

  # Count border spots based on region parameter
  if (region == "inner") {
    # Inner: cluster -> interface only
    idx_count <- (border_df$cluster == cluster & border_df$interface == interface)
    border_count <- sum(idx_count)
  } else if (region == "outer") {
    # Outer: interface -> cluster only (reverse direction)
    idx_count <- (border_df$cluster == interface & border_df$interface == cluster)
    border_count <- sum(idx_count)
  } else {
    # Both: cluster -> interface AND interface -> cluster
    idx_inner <- (border_df$cluster == cluster & border_df$interface == interface)
    idx_outer <- (border_df$cluster == interface & border_df$interface == cluster)
    border_count <- sum(idx_inner) + sum(idx_outer)
  }

  # Get all spots from cluster
  idx_cluster <- df[[cluster_col]] == cluster
  all_from <- df[idx_cluster, ]

  # Get all border spot IDs
  border_spot_ids <- unique(border_df$spot_id)

  # Keep only inner spots (not in border)
  idx_inner <- !(all_from$spot_id %in% border_spot_ids)
  inner_candidates <- all_from[idx_inner, ]

  n_inner <- nrow(inner_candidates)

  # If not enough inner spots, return all with warning
  if (n_inner < border_count) {
    warn <- paste0("Not enough inner spots. Requested: ", border_count,
                   ", Available: ", n_inner, ". Returning all available.")
    warning(warn)
    sample_n <- n_inner
  } else {
    sample_n <- border_count
  }

  # Sample and return
  if (sample_n > 0) {
    idx_sample <- sample(seq_len(nrow(inner_candidates)), size = sample_n, replace = FALSE)
    out <- inner_candidates[idx_sample, ]
  } else {
    out <- inner_candidates[0, ]
  }

  out$interface <- interface
  out$is_inner <- TRUE

  out$region <- region
  out
}

