testthat::test_that("compute_centroids returns mean x/y per cluster", {
  df <- data.frame(
    cluster = c(1, 1, 2, 2),
    x = c(0, 2, 10, 12),
    y = c(1, 3, 5, 7)
  )

  res <- compute_centroids(df)

  # order may follow cluster; keep robust anyway
  res <- res[order(res$cluster), , drop = FALSE]

  testthat::expect_equal(res$cluster, c(1, 2))
  testthat::expect_equal(res$x, c(1, 11))
  testthat::expect_equal(res$y, c(2, 6))
})

testthat::test_that("bresenham_line produces correct grid points on simple lines", {
  # Horizontal line from (0,0) to (3,0)
  p0 <- data.frame(x = 0, y = 0)
  p1 <- data.frame(x = 3, y = 0)
  res <- bresenham_line(p0, p1)

  testthat::expect_equal(res, data.frame(x = 0:3, y = c(0, 0, 0, 0)))

  # Vertical line from (2,1) to (2,4)
  p0 <- data.frame(x = 2, y = 1)
  p1 <- data.frame(x = 2, y = 4)
  res <- bresenham_line(p0, p1)

  testthat::expect_equal(res, data.frame(x = c(2, 2, 2, 2), y = 1:4))
})

testthat::test_that("bresenham_line snaps/does not snap endpoints as requested", {
  p0 <- data.frame(x = 0.49, y = 0.49)
  p1 <- data.frame(x = 2.51, y = 0.49)

  # snap = TRUE: rounds (0,0) -> (3,0)
  res_snap <- bresenham_line(p0, p1, snap = TRUE)
  testthat::expect_equal(res_snap, data.frame(x = 0:3, y = c(0, 0, 0, 0)))

  # snap = FALSE: integer truncation -> (0,0) -> (2,0)
  res_nosnap <- bresenham_line(p0, p1, snap = FALSE)
  testthat::expect_equal(res_nosnap, data.frame(x = 0:2, y = c(0, 0, 0)))
})

testthat::test_that("point_segment_distance_vec matches known distances", {
  # Segment A(0,0) to B(3,0)
  ax <- 0; ay <- 0; bx <- 3; by <- 0

  # Point above middle
  d <- point_segment_distance_vec(px = 1, py = 2, ax, ay, bx, by)
  testthat::expect_equal(d, 2)

  # Point beyond segment endpoint: closest to B(3,0)
  d2 <- point_segment_distance_vec(px = 10, py = 0, ax, ay, bx, by)
  testthat::expect_equal(d2, 7)

  # Vectorized input
  px <- c(0, 1, 2, 3)
  py <- c(1, 1, 1, 1)
  dv <- point_segment_distance_vec(px, py, ax, ay, bx, by)
  testthat::expect_equal(dv, rep(1, 4))
})

testthat::test_that("point_segment_distance_vec handles degenerate segment (A==B)", {
  d <- point_segment_distance_vec(px = c(0, 3), py = c(4, 0), ax = 1, ay = 2, bx = 1, by = 2)
  # distances to A(1,2): sqrt((x-1)^2+(y-2)^2)
  testthat::expect_equal(d, c(sqrt((0 - 1)^2 + (4 - 2)^2), sqrt((3 - 1)^2 + (0 - 2)^2)))
})

testthat::test_that("spots_near_segment returns selected + segment_df and orders by t", {
  df <- data.frame(
    x = 0:4,
    y = c(0, 1, 1, 2, 4),
    id = letters[1:5]
  )
  p0 <- data.frame(x = 0, y = 0)
  p1 <- data.frame(x = 4, y = 0)

  res <- spots_near_segment(df, p0, p1, top_n = 3)

  testthat::expect_true(is.list(res))
  testthat::expect_true(all(c("selected", "segment_df") %in% names(res)))

  sel <- res$selected
  testthat::expect_true(all(c("dist_to_seg", "t") %in% names(sel)))

  # ordered by t
  testthat::expect_true(all(diff(sel$t) >= 0))

  # top_n = 3
  testthat::expect_lte(nrow(sel), 3)

  # segment_df endpoints correct
  testthat::expect_equal(res$segment_df, data.frame(x = c(0, 4), y = c(0, 0)))
})

testthat::test_that("spots_near_segment errors on zero-length segment", {
  df <- data.frame(x = 1:3, y = 1:3)
  p0 <- data.frame(x = 0, y = 0)
  p1 <- data.frame(x = 0, y = 0)

  testthat::expect_error(
    spots_near_segment(df, p0, p1),
    "zero-length segment|coincide|Invalid segment"
  )
})

testthat::test_that("remove_used_points removes rows matching rounded (x,y) keys", {
  df <- data.frame(
    x = c(0.1, 1.2, 1.6, 2.0),
    y = c(0.2, 3.4, 3.4, 9.5),
    id = 1:4
  )
  used_df <- data.frame(
    x = c(1.49, 2.1),
    y = c(3.49, 9.49)
  )

  res <- remove_used_points(df, used_df)

  # rows rounding to (1,3) and (2,9) removed
  # df rows: (0,0), (1,3), (2,3), (2,10)
  # used:    (1,3), (2,9) so only (1,3) removed here; (2,9) matches none
  testthat::expect_false(any(round(res$x) == 1 & round(res$y) == 3))
  testthat::expect_equal(res$id, c(1, 3, 4))
})

testthat::test_that("closest_spot returns the nearest row", {
  df <- data.frame(x = c(0, 2, 5), y = c(0, 2, 1), id = c("a", "b", "c"))

  res <- closest_spot(df, tx = 1.9, ty = 2.1)
  testthat::expect_equal(res$id, "b")
  testthat::expect_equal(nrow(res), 1)
})

testthat::test_that("adjacent_endpoint steps left/right as expected on horizontal segment", {
  df <- data.frame(
    x = c(1, 1, 1),
    y = c(-1, 0, 1),
    id = c("down", "mid", "up")
  )
  A <- data.frame(x = 0, y = 0)
  B <- data.frame(x = 2, y = 0)
  endpoint <- data.frame(x = 1, y = 0)

  # left of A->B should be +y
  res_left <- adjacent_endpoint(df, endpoint, A, B, spacing = 1, side = "left")
  testthat::expect_equal(res_left$id, "up")

  # right of A->B should be -y
  res_right <- adjacent_endpoint(df, endpoint, A, B, spacing = 1, side = "right")
  testthat::expect_equal(res_right$id, "down")
})

testthat::test_that("estimate_spot_spacing returns ~1 on unit grid", {
  set.seed(123)
  df <- data.frame(
    x = rep(1:5, each = 5),
    y = rep(1:5, times = 5)
  )

  s <- estimate_spot_spacing(df, sample_n = 1000)
  testthat::expect_equal(s, 1)
})

testthat::test_that("unit_normal_left returns correct unit normal", {
  A <- data.frame(x = 0, y = 0)
  B <- data.frame(x = 2, y = 0)

  n <- unit_normal_left(A, B)
  testthat::expect_equal(unname(n["nx"]), 0)
  testthat::expect_equal(unname(n["ny"]), 1)

  # length is 1
  testthat::expect_equal(sqrt(sum(unname(n)^2)), 1)
})

testthat::test_that("shift_point shifts coordinates by offset * (nx, ny)", {
  P <- data.frame(x = 1, y = 2)

  P2 <- shift_point(P, nx = 0, ny = 1, offset = 3)
  testthat::expect_equal(P2, data.frame(x = 1, y = 5))

  P3 <- shift_point(P, nx = 1, ny = 0, offset = -2)
  testthat::expect_equal(P3, data.frame(x = -1, y = 2))
})

testthat::test_that("build_one_line returns spots ordered by t", {
  df <- data.frame(
    x = c(0, 1, 2, 3, 4, 2),
    y = c(0, 0, 0, 0, 0, 1),
    id = 1:6
  )
  A <- data.frame(x = 0, y = 0)
  B <- data.frame(x = 4, y = 0)

  sel <- build_one_line(df, A, B, top_n = 5)
  testthat::expect_true("t" %in% names(sel))
  testthat::expect_true(all(diff(sel$t) >= 0))
})

testthat::test_that("build_parallel_lines returns expected structure and line_ids", {
  set.seed(42)

  df <- data.frame(
    x = rep(1:10, each = 3),
    y = rep(1:3, times = 10),
    id = seq_len(30),
    cluster = rep(c("A", "B", "C"), times = 10)
  )
  A <- data.frame(x = 1, y = 2)
  B <- data.frame(x = 10, y = 2)

  out <- build_parallel_lines(df, A, B, top_n = 5, n_extra = 1, side = "both")

  testthat::expect_true(is.list(out))
  testthat::expect_true(all(c("lines", "spacing", "lane_width") %in% names(out)))
  testthat::expect_true(is.data.frame(out$lines))
  testthat::expect_true(all(c("line_id", "offset") %in% names(out$lines)))

  ids <- unique(out$lines$line_id)
  # should contain center and up to one left/right line
  testthat::expect_true("center" %in% ids)
  testthat::expect_true(any(grepl("^left_1$", ids)) || any(grepl("^right_1$", ids)) || length(ids) >= 1)
})

testthat::test_that("filter_out_by_endpoint_clusters keeps only matching line_ids", {
  out <- list(
    lines = data.frame(
      line_id = c("L1","L1","L2","L2"),
      t       = c(0.0, 1.0, 0.0, 1.0),
      cluster = c("A", "B", "A", "C"),
      x = 1:4, y = 1:4
    )
  )

  res <- filter_out_by_endpoint_clusters(out, allowed_start_clusters = "A", allowed_end_clusters = "B")

  testthat::expect_true(all(res$lines$line_id %in% "L1"))
  testthat::expect_false(any(res$lines$line_id %in% "L2"))
})
