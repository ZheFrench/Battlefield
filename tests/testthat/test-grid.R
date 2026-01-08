test_that("detect_grid_type identifies square grids correctly", {
  # Create a simple square grid
  sq <- expand.grid(x = 0:9, y = 0:9)
  
  result <- detect_grid_type(sq, coords = c("x", "y"), k = 12, verbose = FALSE)
  
  expect_true(is.list(result))
  expect_equal(result$grid_type, "square")
  expect_true(is.numeric(result$step))
  expect_true(result$step > 0)
  expect_true(is.numeric(result$median_ratio))
})

test_that("detect_grid_type identifies hexagonal grids correctly", {
  # Create a hexagonal grid
  nx <- 10
  ny <- 10
  hex <- expand.grid(i = 0:(nx - 1), j = 0:(ny - 1))
  hex$x <- hex$i + 0.5 * (hex$j %% 2)
  hex$y <- (sqrt(3) / 2) * hex$j
  hex <- hex[, c("x", "y")]
  
  result <- detect_grid_type(hex, coords = c("x", "y"), k = 12, verbose = FALSE)
  
  expect_true(is.list(result))
  expect_equal(result$grid_type, "hexagonal")
  expect_true(is.numeric(result$step))
  expect_true(result$step > 0)
  expect_true(is.numeric(result$median_ratio))
})

test_that("detect_grid_type returns correct list structure", {
  sq <- expand.grid(x = 0:5, y = 0:5)
  
  result <- detect_grid_type(sq, coords = c("x", "y"), k = 6, verbose = FALSE)
  
  expect_named(result, c("grid_type", "step", "median_ratio"))
  expect_type(result$grid_type, "character")
  expect_type(result$step, "double")
  expect_type(result$median_ratio, "double")
})

test_that("detect_grid_type validates coordinate columns", {
  df <- data.frame(a = 1:5, b = 1:5)
  
  expect_error(detect_grid_type(df, coords = c("x", "y"), verbose = FALSE))
})

test_that("detect_grid_type works with different k values", {
  sq <- expand.grid(x = 0:9, y = 0:9)
  
  result_k6 <- detect_grid_type(sq, coords = c("x", "y"), k = 6, verbose = FALSE)
  result_k12 <- detect_grid_type(sq, coords = c("x", "y"), k = 12, verbose = FALSE)
  
  expect_equal(result_k6$grid_type, "square")
  expect_equal(result_k12$grid_type, "square")
})

test_that("detect_grid_type handles tolerance parameter", {
  sq <- expand.grid(x = 0:9, y = 0:9)
  
  result_strict <- detect_grid_type(sq, coords = c("x", "y"), tolerance = 0.01, verbose = FALSE)
  result_lenient <- detect_grid_type(sq, coords = c("x", "y"), tolerance = 0.2, verbose = FALSE)
  
  expect_true(is.character(result_strict$grid_type))
  expect_true(is.character(result_lenient$grid_type))
})

test_that("get_neighborhood_params returns correct structure for square grid", {
  sq <- expand.grid(x = 0:9, y = 0:9)
  
  params <- get_neighborhood_params(sq, square_connectivity = 4, verbose = FALSE)
  
  expect_true(is.list(params))
  expect_named(params, c("grid_type", "connectivity", "radius", "k", "comment"))
  expect_equal(params$grid_type, "square")
  expect_equal(params$connectivity, 4)
  expect_true(is.numeric(params$radius))
  expect_true(params$radius > 0)
  expect_true(is.numeric(params$k))
})

test_that("get_neighborhood_params returns correct structure for hexagonal grid", {
  nx <- 10
  ny <- 10
  hex <- expand.grid(i = 0:(nx - 1), j = 0:(ny - 1))
  hex$x <- hex$i + 0.5 * (hex$j %% 2)
  hex$y <- (sqrt(3) / 2) * hex$j
  hex <- hex[, c("x", "y")]
  
  params <- get_neighborhood_params(hex, verbose = FALSE)
  
  expect_equal(params$grid_type, "hexagonal")
  expect_equal(params$connectivity, 6)
  expect_true(is.numeric(params$radius))
  expect_true(params$radius > 0)
})

test_that("get_neighborhood_params respects square_connectivity parameter", {
  sq <- expand.grid(x = 0:9, y = 0:9)
  
  params_4 <- get_neighborhood_params(sq, square_connectivity = 4, verbose = FALSE)
  params_8 <- get_neighborhood_params(sq, square_connectivity = 8, verbose = FALSE)
  
  expect_equal(params_4$connectivity, 4)
  expect_equal(params_8$connectivity, 8)
  # 8-connectivity should have larger radius
  expect_gt(params_8$radius, params_4$radius)
})

test_that("get_neighborhood_params returns different k values for connectivity", {
  sq <- expand.grid(x = 0:9, y = 0:9)
  
  params_4 <- get_neighborhood_params(sq, square_connectivity = 4, verbose = FALSE)
  params_8 <- get_neighborhood_params(sq, square_connectivity = 8, verbose = FALSE)
  
  expect_lt(params_4$k, params_8$k)
})

test_that("get_neighborhood_params includes descriptive comment", {
  sq <- expand.grid(x = 0:9, y = 0:9)
  
  params <- get_neighborhood_params(sq, verbose = FALSE)
  
  expect_true(is.character(params$comment))
  expect_true(nchar(params$comment) > 0)
})

test_that("get_neighborhood_params validates coordinate columns", {
  df <- data.frame(a = 1:5, b = 1:5)
  
  expect_error(get_neighborhood_params(df, verbose = FALSE))
})

test_that("detect_grid_type verbose output doesn't affect results", {
  sq <- expand.grid(x = 0:5, y = 0:5)
  
  result_quiet <- detect_grid_type(sq, verbose = FALSE)
  result_verbose <- detect_grid_type(sq, verbose = TRUE)
  
  expect_equal(result_quiet$grid_type, result_verbose$grid_type)
  expect_equal(result_quiet$step, result_verbose$step)
  expect_equal(result_quiet$median_ratio, result_verbose$median_ratio)
})

test_that("get_neighborhood_params verbose output doesn't affect results", {
  sq <- expand.grid(x = 0:5, y = 0:5)
  
  params_quiet <- get_neighborhood_params(sq, verbose = FALSE)
  params_verbose <- get_neighborhood_params(sq, verbose = TRUE)
  
  expect_equal(params_quiet$grid_type, params_verbose$grid_type)
  expect_equal(params_quiet$connectivity, params_verbose$connectivity)
  expect_equal(params_quiet$radius, params_verbose$radius)
})

test_that("detect_grid_type step size is reasonable", {
  # Square grid with step size 1
  sq <- expand.grid(x = 0:9, y = 0:9)
  result <- detect_grid_type(sq, verbose = FALSE)
  
  # Step should be close to 1 for a unit-spaced square grid
  expect_true(result$step > 0.9)
  expect_true(result$step < 1.1)
})

test_that("detect_grid_type ratio is reasonable for square grid", {
  sq <- expand.grid(x = 0:9, y = 0:9)
  result <- detect_grid_type(sq, tolerance = 0.15, verbose = FALSE)
  
  # For square grid, ratio should be around sqrt(2) ~= 1.414
  sqrt2 <- sqrt(2)
  expect_true(result$median_ratio > sqrt2 - 0.2)
  expect_true(result$median_ratio < sqrt2 + 0.2)
})

test_that("detect_grid_type ratio is reasonable for hexagonal grid", {
  nx <- 10
  ny <- 10
  hex <- expand.grid(i = 0:(nx - 1), j = 0:(ny - 1))
  hex$x <- hex$i + 0.5 * (hex$j %% 2)
  hex$y <- (sqrt(3) / 2) * hex$j
  hex <- hex[, c("x", "y")]
  
  result <- detect_grid_type(hex, tolerance = 0.15, verbose = FALSE)
  
  # For hexagonal grid, ratio should be around sqrt(3) ~= 1.732
  sqrt3 <- sqrt(3)
  expect_true(result$median_ratio > sqrt3 - 0.2)
  expect_true(result$median_ratio < sqrt3 + 0.2)
})
