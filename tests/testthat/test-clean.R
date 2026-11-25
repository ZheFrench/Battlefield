test_that("clean_orphans removes small pieces of the cluster only", {
g <- expand.grid(x = 1:20, y = 1:10)
g$cluster <- "B"
g$cluster[g$x <= 4] <- "A"
g$cluster[g$x == 10 & g$y == 5] <- "A"
res <- clean_orphans(g, clusters = "A", min_size = 5)
expect_equal(sum(res$cluster == "A"), 40)
expect_equal(sum(res$cluster == "B"), sum(g$cluster == "B"))
expect_equal(names(res), names(g))
expect_equal(clean_orphans(g, "B", min_size = 5), g)
expect_warning(same <- clean_orphans(g, "A", min_size = 1000), "nothing removed")
expect_equal(same, g)
})
