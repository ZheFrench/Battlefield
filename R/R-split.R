#' Split one scattered cluster into contiguous pieces
#'
#' Spots of the cluster that touch each other form a piece. Pieces with at least
#' `min_size` spots are renamed `<cluster>_1`, `<cluster>_2`, ... from largest
#' to smallest. Smaller pieces (orphan spots) keep the original cluster name.
#' A cluster made of a single piece, and all other clusters, are unchanged.
#'
#' @param df A data.frame with coordinates and a cluster column.
#' @param clusters The cluster to split (a single label).
#' @param min_size Minimal number of spots for a piece to be renamed.
#' @param max_dist Maximal distance between touching spots. `NULL` uses
#'   1.5 times the spot spacing (direct grid neighbors).
#' @param coord_cols,cluster_col Names of the coordinate and cluster columns.
#'
#' @return `df` with an additional character column `cluster_split`.
#'
#' @details
#' Each call starts from `cluster_col`. To split several clusters, use
#' `cluster_col = "cluster_split"` from the second call. Use
#' [add_split_to_spe()] to store the result in a SpatialExperiment object.
#'
#' @examples
#' df <- expand.grid(x = 1:20, y = 1:10)
#' df$spot_id <- paste0("s", seq_len(nrow(df)))
#' df$cluster <- "B"
#' df$cluster[df$x <= 4 | df$x >= 18] <- "A"   # two pieces
#' df$cluster[df$x == 10 & df$y == 5] <- "A"   # one orphan spot
#' split <- split_clusters(df, "A", min_size = 5)
#' table(split$cluster_split)
#'
#' @seealso [clean_orphans()], [add_split_to_spe()]
#' @importFrom RANN nn2
#' @importFrom igraph graph_from_edgelist components
#' @export
split_clusters <- function(df, clusters, min_size = 20,
    max_dist = NULL, coord_cols = c("x", "y"), cluster_col = "cluster") {
    stopifnot(all(c(coord_cols, cluster_col) %in% colnames(df)),
        length(clusters) == 1L, min_size >= 1)
    xy <- as.matrix(df[, coord_cols])
    cl <- as.character(df[[cluster_col]])
    k <- as.character(clusters)
    i <- which(cl == k)
    if (length(i) == 0L) stop("No spots found in cluster ", k, ".")
    df$cluster_split <- cl
    if (length(i) == 1L) return(df)
    if (is.null(max_dist)) {
        # All spots: subsampling would overestimate the spacing
        max_dist <- 1.5 * stats::median(RANN::nn2(xy, xy, k = 2)$nn.dists[, 2])
    }
    nn <- RANN::nn2(xy[i, , drop = FALSE], k = min(9L, length(i)))
    touch <- which(nn$nn.dists <= max_dist, arr.ind = TRUE)
    # self edges keep isolated spots in the graph
    edges <- rbind(cbind(seq_along(i), seq_along(i)),
        cbind(touch[, 1], nn$nn.idx[touch]))
    comp <- igraph::components(igraph::graph_from_edgelist(edges, directed = FALSE))
    if (comp$no == 1L) return(df)
    big <- which(comp$csize >= min_size)
    big <- big[order(-comp$csize[big])]
    rank <- match(comp$membership, big)
    ok <- !is.na(rank)
    df$cluster_split[i[ok]] <- paste0(k, "_", rank[ok])
    df
}
