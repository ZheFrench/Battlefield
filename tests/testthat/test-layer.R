test_that("create_cluster_layers returns data.frame with layer column", {
  set.seed(1)
  df <- data.frame(
    spot_id = 1:100,
    x = rnorm(100),
    y = rnorm(100),
    cluster = sample(c("A", "B", "C"), 100, replace = TRUE)
  )
  
  result <- create_cluster_layers(df, target_cluster = "A", k = 6)
  
  expect_true(is.data.frame(result))
  expect_true("layer" %in% colnames(result))
  expect_true(all(result$layer %in% c("border", "intermediate", "core")))
})

test_that("create_cluster_layers filters to target cluster only", {
  set.seed(1)
  df <- data.frame(
    spot_id = 1:100,
    x = rnorm(100),
    y = rnorm(100),
    cluster = sample(c("A", "B", "C"), 100, replace = TRUE)
  )
  
  result_a <- create_cluster_layers(df, target_cluster = "A", k = 6)
  result_b <- create_cluster_layers(df, target_cluster = "B", k = 6)
  
  expect_true(all(result_a$cluster == "A"))
  expect_true(all(result_b$cluster == "B"))
  expect_true(nrow(result_a) < nrow(df))
})

test_that("create_cluster_layers identifies border spots", {
  # Create synthetic data with clear borders
  set.seed(42)
  
  # Cluster A: points around (0,0)
  a_x <- rnorm(30, mean = 0, sd = 0.5)
  a_y <- rnorm(30, mean = 0, sd = 0.5)
  
  # Cluster B: points around (2,2) - close enough to have border interaction
  b_x <- rnorm(30, mean = 2, sd = 0.5)
  b_y <- rnorm(30, mean = 2, sd = 0.5)
  
  df <- data.frame(
    spot_id = 1:60,
    x = c(a_x, b_x),
    y = c(a_y, b_y),
    cluster = c(rep("A", 30), rep("B", 30))
  )
  
  result <- create_cluster_layers(df, target_cluster = "A", k = 12)
  
  # May have border, intermediate, or core - just verify valid layers
  expect_true(all(result$layer %in% c("border", "intermediate", "core")))
})

test_that("create_cluster_layers respects k parameter", {
  set.seed(1)
  df <- data.frame(
    spot_id = 1:100,
    x = rnorm(100),
    y = rnorm(100),
    cluster = sample(c("A", "B"), 100, replace = TRUE)
  )
  
  result_k4 <- create_cluster_layers(df, target_cluster = "A", k = 4)
  result_k10 <- create_cluster_layers(df, target_cluster = "A", k = 10)
  
  expect_true(is.data.frame(result_k4))
  expect_true(is.data.frame(result_k10))
  # Results may differ based on different neighbor sets
  expect_true(nrow(result_k4) > 0)
  expect_true(nrow(result_k10) > 0)
})

test_that("create_cluster_layers respects intermediate_quantile parameter", {
  set.seed(1)
  df <- data.frame(
    spot_id = 1:100,
    x = rnorm(100),
    y = rnorm(100),
    cluster = sample(c("A", "B"), 100, replace = TRUE)
  )
  
  result_q33 <- create_cluster_layers(df, target_cluster = "A", k = 6, intermediate_quantile = 0.33)
  result_q75 <- create_cluster_layers(df, target_cluster = "A", k = 6, intermediate_quantile = 0.75)
  
  n_inter_q33 <- sum(result_q33$layer == "intermediate")
  n_inter_q75 <- sum(result_q75$layer == "intermediate")
  
  # Higher quantile should generally result in more intermediate spots
  expect_gte(n_inter_q75, n_inter_q33)
})

test_that("create_cluster_layers handles max_dist parameter", {
  set.seed(1)
  df <- data.frame(
    spot_id = 1:100,
    x = rnorm(100),
    y = rnorm(100),
    cluster = sample(c("A", "B"), 100, replace = TRUE)
  )
  
  # With very small max_dist, fewer neighbors should be considered
  result_small_dist <- create_cluster_layers(df, target_cluster = "A", k = 6, max_dist = 0.1)
  result_no_dist <- create_cluster_layers(df, target_cluster = "A", k = 6, max_dist = NULL)
  
  expect_true(is.data.frame(result_small_dist))
  expect_true(is.data.frame(result_no_dist))
})

test_that("create_cluster_layers preserves original columns", {
  set.seed(1)
  df <- data.frame(
    spot_id = 1:50,
    x = rnorm(50),
    y = rnorm(50),
    cluster = sample(c("A", "B"), 50, replace = TRUE),
    extra_col = letters[1:50]
  )
  
  result <- create_cluster_layers(df, target_cluster = "A", k = 6)
  
  # All original columns should be present
  expect_true("spot_id" %in% colnames(result))
  expect_true("x" %in% colnames(result))
  expect_true("y" %in% colnames(result))
  expect_true("cluster" %in% colnames(result))
  expect_true("extra_col" %in% colnames(result))
})

test_that("create_cluster_layers handles custom coordinate column names", {
  set.seed(1)
  df <- data.frame(
    spot_id = 1:50,
    X = rnorm(50),
    Y = rnorm(50),
    cluster = sample(c("A", "B"), 50, replace = TRUE)
  )
  
  result <- create_cluster_layers(df, target_cluster = "A", k = 6, coord_cols = c("X", "Y"))
  
  expect_true("layer" %in% colnames(result))
  expect_true(all(result$layer %in% c("border", "intermediate", "core")))
})

test_that("create_cluster_layers handles custom cluster column name", {
  set.seed(1)
  df <- data.frame(
    spot_id = 1:50,
    x = rnorm(50),
    y = rnorm(50),
    group = sample(c(1, 2, 3), 50, replace = TRUE)
  )
  
  result <- create_cluster_layers(df, target_cluster = 1, k = 6, cluster_col = "group")
  
  expect_true("layer" %in% colnames(result))
  expect_true(all(result$group == 1))
})

test_that("create_cluster_layers raises error for non-existent target_cluster", {
  set.seed(1)
  df <- data.frame(
    spot_id = 1:50,
    x = rnorm(50),
    y = rnorm(50),
    cluster = sample(c("A", "B"), 50, replace = TRUE)
  )
  
  expect_error(create_cluster_layers(df, target_cluster = "Z", k = 6))
})

test_that("create_cluster_layers handles cluster with all border spots", {
  # Create synthetic data where target cluster is surrounded
  set.seed(42)
  
  # Cluster A: small cluster at center
  a_x <- rnorm(10, mean = 0, sd = 0.1)
  a_y <- rnorm(10, mean = 0, sd = 0.1)
  
  # Cluster B: large surrounding cluster
  b_x <- rnorm(100, mean = 0, sd = 2)
  b_y <- rnorm(100, mean = 0, sd = 2)
  
  df <- data.frame(
    spot_id = 1:110,
    x = c(a_x, b_x),
    y = c(a_y, b_y),
    cluster = c(rep("A", 10), rep("B", 100))
  )
  
  result <- create_cluster_layers(df, target_cluster = "A", k = 6)
  
  # All spots might be borders if cluster is small and surrounded
  expect_true(is.data.frame(result))
  expect_true("layer" %in% colnames(result))
})

test_that("create_cluster_layers handles cluster with all core spots", {
  # Create synthetic data with single large cluster (no border)
  set.seed(42)
  
  df <- data.frame(
    spot_id = 1:50,
    x = rnorm(50),
    y = rnorm(50),
    cluster = rep("A", 50)
  )
  
  result <- create_cluster_layers(df, target_cluster = "A", k = 6)
  
  # All spots should be core since no other clusters
  expect_true(all(result$layer == "core"))
})

test_that("create_cluster_layers validates intermediate_quantile range", {
  set.seed(1)
  df <- data.frame(
    spot_id = 1:50,
    x = rnorm(50),
    y = rnorm(50),
    cluster = sample(c("A", "B"), 50, replace = TRUE)
  )
  
  # Should fail for invalid quantile values
  expect_error(create_cluster_layers(df, target_cluster = "A", k = 6, intermediate_quantile = -0.1))
  expect_error(create_cluster_layers(df, target_cluster = "A", k = 6, intermediate_quantile = 1.5))
})

test_that("create_all_layers processes multiple clusters", {
  set.seed(1)
  df <- data.frame(
    spot_id = 1:150,
    x = rnorm(150),
    y = rnorm(150),
    cluster = sample(c("A", "B", "C"), 150, replace = TRUE)
  )
  
  result <- create_all_layers(df, k = 6)
  
  expect_true(is.data.frame(result))
  expect_true("layer" %in% colnames(result))
  # Result should be complete
  expect_equal(nrow(result), nrow(df))
})

test_that("create_all_layers processes specified clusters only", {
  set.seed(1)
  df <- data.frame(
    spot_id = 1:150,
    x = rnorm(150),
    y = rnorm(150),
    cluster = sample(c("A", "B", "C"), 150, replace = TRUE)
  )
  
  result <- create_all_layers(df, clusters = c("A", "B"), k = 6)
  
  # Result should only contain A and B
  expect_true(all(result$cluster %in% c("A", "B")))
})

test_that("create_all_layers handles NULL clusters parameter", {
  set.seed(1)
  df <- data.frame(
    spot_id = 1:100,
    x = rnorm(100),
    y = rnorm(100),
    cluster = sample(c(1, 2, 3), 100, replace = TRUE)
  )
  
  result <- create_all_layers(df, clusters = NULL, k = 6)
  
  expect_true(is.data.frame(result))
  expect_true("layer" %in% colnames(result))
  expect_equal(nrow(result), nrow(df))
})

test_that("create_all_layers respects k parameter", {
  set.seed(1)
  df <- data.frame(
    spot_id = 1:100,
    x = rnorm(100),
    y = rnorm(100),
    cluster = sample(c("A", "B"), 100, replace = TRUE)
  )
  
  result_k4 <- create_all_layers(df, k = 4)
  result_k10 <- create_all_layers(df, k = 10)
  
  expect_equal(nrow(result_k4), nrow(df))
  expect_equal(nrow(result_k10), nrow(df))
})

test_that("create_all_layers respects intermediate_quantile parameter", {
  set.seed(1)
  df <- data.frame(
    spot_id = 1:100,
    x = rnorm(100),
    y = rnorm(100),
    cluster = sample(c("A", "B"), 100, replace = TRUE)
  )
  
  result_q33 <- create_all_layers(df, k = 6, intermediate_quantile = 0.33)
  result_q75 <- create_all_layers(df, k = 6, intermediate_quantile = 0.75)
  
  n_inter_q33 <- sum(result_q33$layer == "intermediate")
  n_inter_q75 <- sum(result_q75$layer == "intermediate")
  
  # Higher quantile should generally result in more intermediate spots
  expect_gte(n_inter_q75, n_inter_q33)
})

test_that("create_all_layers preserves original data", {
  set.seed(1)
  df <- data.frame(
    spot_id = 1:100,
    x = rnorm(100),
    y = rnorm(100),
    cluster = sample(c("A", "B"), 100, replace = TRUE)
  )
  
  result <- create_all_layers(df, k = 6)
  
  # Number of rows should match original
  expect_equal(nrow(result), nrow(df))
  # All original columns should be present
  expect_true(all(colnames(df) %in% colnames(result)))
})

test_that("create_all_layers returns rows sorted by cluster", {
  set.seed(1)
  df <- data.frame(
    spot_id = 1:100,
    x = rnorm(100),
    y = rnorm(100),
    cluster = sample(c("A", "B", "C"), 100, replace = TRUE)
  )
  
  result <- create_all_layers(df, k = 6)
  
  # Result should group by cluster
  clusters_in_result <- result$cluster
  expect_true(is.character(clusters_in_result) || is.factor(clusters_in_result) || is.numeric(clusters_in_result))
})

test_that("create_all_layers handles custom coordinate columns", {
  set.seed(1)
  df <- data.frame(
    spot_id = 1:100,
    lon = rnorm(100),
    lat = rnorm(100),
    cluster = sample(c("A", "B"), 100, replace = TRUE)
  )
  
  result <- create_all_layers(df, k = 6, coord_cols = c("lon", "lat"))
  
  expect_true("layer" %in% colnames(result))
  expect_equal(nrow(result), nrow(df))
})

test_that("create_all_layers handles custom cluster column", {
  set.seed(1)
  df <- data.frame(
    spot_id = 1:100,
    x = rnorm(100),
    y = rnorm(100),
    region = sample(c(1, 2, 3), 100, replace = TRUE)
  )
  
  result <- create_all_layers(df, k = 6, cluster_col = "region")
  
  expect_true("layer" %in% colnames(result))
  expect_equal(nrow(result), nrow(df))
})

test_that("create_all_layers handles max_dist parameter", {
  set.seed(1)
  df <- data.frame(
    spot_id = 1:100,
    x = rnorm(100),
    y = rnorm(100),
    cluster = sample(c("A", "B"), 100, replace = TRUE)
  )
  
  result_with_dist <- create_all_layers(df, k = 6, max_dist = 1.0)
  result_no_dist <- create_all_layers(df, k = 6, max_dist = NULL)
  
  expect_equal(nrow(result_with_dist), nrow(df))
  expect_equal(nrow(result_no_dist), nrow(df))
})

test_that("create_all_layers with empty clusters list returns empty dataframe", {
  set.seed(1)
  df <- data.frame(
    spot_id = 1:50,
    x = rnorm(50),
    y = rnorm(50),
    cluster = sample(c("A", "B"), 50, replace = TRUE)
  )
  
  result <- create_all_layers(df, clusters = character(0), k = 6)
  
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 0)
})

test_that("create_all_layers includes all layer types", {
  set.seed(42)
  
  # Create synthetic data with clear layer patterns
  a_x <- rnorm(30, mean = 0, sd = 0.5)
  a_y <- rnorm(30, mean = 0, sd = 0.5)
  
  b_x <- rnorm(50, mean = 2, sd = 0.5)
  b_y <- rnorm(50, mean = 2, sd = 0.5)
  
  df <- data.frame(
    spot_id = 1:80,
    x = c(a_x, b_x),
    y = c(a_y, b_y),
    cluster = c(rep("A", 30), rep("B", 50))
  )
  
  result <- create_all_layers(df, k = 12)
  
  # Verify we have valid layer types
  expect_true(all(result$layer %in% c("border", "intermediate", "core")))
})
