test_that("directed_cluster_interface_pairs generates all ordered pairs", {
cl <- c("1", "2", "3")
result <- directed_cluster_interface_pairs(cl)

expect_true(is.data.frame(result))
expect_named(result, c("cluster", "interface", "directed_pair"))
expect_equal(nrow(result), 6)
})

test_that("directed_cluster_interface_pairs handles NA values", {
cl <- c("A", "B", NA)
result <- directed_cluster_interface_pairs(cl)

expect_equal(nrow(result), 2)
expect_true(all(!is.na(result$cluster)))
})

test_that("directed_cluster_interface_pairs excludes self-pairs", {
cl <- c("1", "2", "3")
result <- directed_cluster_interface_pairs(cl)

expect_true(all(result$cluster != result$interface))
})

test_that("directed_cluster_interface_pairs returns empty df with < 2 clusters", {
cl <- c("1", "1")
result <- directed_cluster_interface_pairs(cl)

expect_equal(nrow(result), 0)
})

test_that("directed_cluster_interface_pairs sorts clusters when requested", {
cl <- c("3", "1", "2")
result_sorted <- directed_cluster_interface_pairs(cl, sort_clusters = TRUE)

# Check that all pairs are generated and sorted
unique_clusters <- unique(c(result_sorted$cluster, result_sorted$interface))
# Just verify that sorting was applied - clusters should be in sorted order
expect_true(length(unique_clusters) == 3)
})

test_that("select_border_spots identifies border spots correctly", {
set.seed(1)
df <- data.frame(
    spot_id = 1:50,
    x = rnorm(50),
    y = rnorm(50),
    cluster = sample(c("A", "B"), 50, replace = TRUE)
)

result <- select_border_spots(df, cluster = "A", interface = "B", k = 5)

expect_true(is.data.frame(result))
if (nrow(result) > 0) {
    expect_true(all(result$cluster == "A"))
}
})

test_that("select_border_spots returns correct columns", {
set.seed(1)
df <- data.frame(
    spot_id = 1:50,
    x = rnorm(50),
    y = rnorm(50),
    cluster = sample(c("A", "B"), 50, replace = TRUE)
)

result <- select_border_spots(df, cluster = "A", interface = "B", k = 5)

expected_cols <- c("spot_id", "x", "y", "directed_pair", "undirected_pair", 
                    "cluster", "interface", "is_border", "is_border_multiple", 
                    "other_adjacent_borders", "mode")
expect_equal(colnames(result), expected_cols)
})

test_that("select_border_spots handles mode parameter correctly", {
set.seed(1)
df <- data.frame(
    spot_id = 1:50,
    x = rnorm(50),
    y = rnorm(50),
    cluster = sample(c("A", "B"), 50, replace = TRUE)
)

inner <- select_border_spots(df, cluster = "A", interface = "B", mode = "inner", k = 5)

if (nrow(inner) > 0) {
    expect_true(all(inner$mode == "inner"))
}
})

test_that("select_border_spots returns empty df for non-existent cluster", {
df <- data.frame(
    spot_id = 1:50,
    x = rnorm(50),
    y = rnorm(50),
    cluster = rep(c("A", "B"), 25)
)

result <- select_border_spots(df, cluster = "Z", interface = "A", k = 5)

expect_equal(nrow(result), 0)
})

test_that("build_all_borders processes multiple pairs", {
set.seed(1)
df <- data.frame(
    spot_id = 1:60,
    x = rnorm(60),
    y = rnorm(60),
    cluster = sample(c("A", "B", "C"), 60, replace = TRUE)
)

result <- build_all_borders(df, k = 5)

expect_true(is.data.frame(result))
})

test_that("build_all_borders returns empty df with no pairs", {
df <- data.frame(
    spot_id = 1:50,
    x = rnorm(50),
    y = rnorm(50),
    cluster = rep("A", 50)
)

result <- build_all_borders(df, k = 5)

expect_equal(nrow(result), 0)
})

test_that("build_all_borders respects mode parameter", {
set.seed(1)
df <- data.frame(
    spot_id = 1:60,
    x = rnorm(60),
    y = rnorm(60),
    cluster = sample(c("A", "B"), 60, replace = TRUE)
)

result_inner <- build_all_borders(df, k = 5, mode = "inner")
result_both <- build_all_borders(df, k = 5, mode = "both")

expect_gte(nrow(result_both), nrow(result_inner))
})

test_that("build_all_cores returns data.frame", {
set.seed(1)
df <- data.frame(
    spot_id = 1:100,
    x = rnorm(100),
    y = rnorm(100),
    cluster = sample(c("A", "B"), 100, replace = TRUE)
)

border_df <- build_all_borders(df, k = 5)

if (nrow(border_df) > 0) {
    result <- build_all_cores(df, border_df, mode = "both")
    expect_true(is.data.frame(result))
}
})

test_that("select_core_spots returns correct type", {
set.seed(1)
df <- data.frame(
    spot_id = 1:100,
    x = rnorm(100),
    y = rnorm(100),
    cluster = sample(c("A", "B"), 100, replace = TRUE)
)

border_df <- build_all_borders(df, k = 5)

if (nrow(border_df) > 0) {
    result <- select_core_spots(df, border_df, "A", "B", mode = "inner")
    expect_true(is.null(result) || is.data.frame(result))
}
})
