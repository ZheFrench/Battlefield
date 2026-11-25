#' Remove the orphan spots of one cluster
#'
#' Removes the spots of the cluster that belong to a contiguous piece of fewer
#' than `min_size` spots. Spots of other clusters are always kept.
#'
#' To clean border spots, use it on the output of [create_cluster_layers()]:
#' layers are then computed with all spots, and labels are not changed.
#'
#' @inheritParams split_clusters
#' @param min_size Pieces with fewer spots are removed.
#'
#' @return `df` without the orphan spots. If every piece is smaller than
#'   `min_size`, `df` is returned unchanged with a warning.
#'
#' @examples
#' df <- expand.grid(x = 1:20, y = 1:10)
#' df$spot_id <- paste0("s", seq_len(nrow(df)))
#' df$cluster <- "B"
#' df$cluster[df$x <= 4 | df$x >= 18] <- "A"   # two pieces
#' df$cluster[df$x == 10 & df$y == 5] <- "A"   # one orphan spot
#' cleaned <- clean_orphans(df, "A", min_size = 5)
#' nrow(df) - nrow(cleaned)
#'
#' # Clean layers: computed with all spots, then orphans are removed
#' layers <- create_cluster_layers(df, "A", k = 4)
#' layers <- clean_orphans(layers, "A", min_size = 5)
#' table(layers$layer)
#'
#' @seealso [split_clusters()]
#' @export
clean_orphans <- function(df, clusters, min_size = 20, max_dist = NULL,
    coord_cols = c("x", "y"), cluster_col = "cluster") {
    sp <- split_clusters(df, clusters, min_size, max_dist, coord_cols,
        cluster_col)$cluster_split
    in_cl <- which(as.character(df[[cluster_col]]) == as.character(clusters))
    orphan <- in_cl[sp[in_cl] == as.character(clusters)]
    # a cluster in a single piece is not renamed and has no orphan
    if (length(orphan) == length(in_cl)) {
        one <- split_clusters(df, clusters, 1, max_dist, coord_cols,
            cluster_col)$cluster_split[in_cl]
        if (length(unique(one)) > 1L) {
            warning("All pieces are smaller than min_size; nothing removed.")
        }
        return(df)
    }
    if (length(orphan) > 0L) df <- df[-orphan, , drop = FALSE]
    df
}
