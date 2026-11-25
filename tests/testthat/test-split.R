test_that("split_clusters renames big pieces and keeps orphans", {
g <- expand.grid(x = 1:20, y = 1:10)
g$cluster <- "B"
g$cluster[g$x <= 4] <- "A"                    # 40 spots
g$cluster[g$x >= 18] <- "A"                   # 30 spots
g$cluster[g$x == 10 & g$y == 5] <- "A"        # orphan
res <- split_clusters(g, clusters = "A", min_size = 5)
expect_equal(sum(res$cluster_split == "A_1"), 40)
expect_equal(sum(res$cluster_split == "A_2"), 30)
expect_equal(sum(res$cluster_split == "A"), 1)
expect_true(all(res$cluster_split[g$cluster == "B"] == "B"))
expect_true(all(split_clusters(g, "B")$cluster_split == g$cluster))
expect_error(split_clusters(g, "Z"))
expect_error(split_clusters(g, c("A", "B")))
})

test_that("split_clusters can be chained and stored in a SpatialExperiment", {
g <- expand.grid(x = 1:20, y = 1:10)
g$spot_id <- paste0("s", seq_len(nrow(g)))
g$cluster <- ifelse(g$x <= 4 | g$x >= 18, "A", "B")
g$cluster[g$x == 10] <- "C"
sp <- split_clusters(g, "A", min_size = 5)
sp <- split_clusters(sp, "B", min_size = 5, cluster_col = "cluster_split")
expect_setequal(sp$cluster_split, c("A_1", "A_2", "B_1", "B_2", "C"))
spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = matrix(1, 1, 200,
        dimnames = list("g", g$spot_id))),
    spatialCoords = as.matrix(g[, c("x", "y")]))
res <- add_split_to_spe(spe, sp[200:1, ])
expect_equal(res$cluster_split, sp$cluster_split)
})
