#' Compute cluster centroids (mean x/y per cluster)
#'
#' Computes the centroid of each cluster by taking the mean of the `x` and `y`
#' coordinates for each `cluster` level.
#'
#' @param df A `data.frame` containing at least the columns `cluster`, `x`, and
#' `y`.
#' `cluster` can be any type accepted by `aggregate()` grouping (e.g., integer,
#'   factor, character). `x` and `y` must be numeric.
#'
#' @return A `data.frame` with one row per cluster and columns:
#' \describe{
#'   \item{cluster}{Cluster identifier.}
#'   \item{x}{Mean x-coordinate for the cluster.}
#'   \item{y}{Mean y-coordinate for the cluster.}
#' }
#'
#' @details
#' This is a thin wrapper around \code{\link[stats:aggregate]{aggregate}} using
#' \code{FUN = mean}. Missing values are handled according to `mean()`'s default
#' behavior (i.e., `NA`s will propagate unless you pre-handle them).
#'
#' @examples
#' df <- data.frame(
#'   cluster = c(1, 1, 2, 2),
#'   x = c(0, 2, 10, 12),
#'   y = c(1, 3, 5, 7)
#' )
#' compute_centroids(df)
#'
#' @importFrom stats aggregate
#' @export
compute_centroids <- function(df) {
    stats::aggregate(cbind(x, y) ~ cluster, data = df, FUN = mean)
}

#' Rasterize a line between two points using Bresenham's algorithm
#'
#' Generates integer grid coordinates along the line segment connecting two
#' points
#' using Bresenham's line algorithm.
#'
#' @param p0 A `data.frame` with at least columns `x` and `y`. If it contains
#'   multiple rows, only the first row is used.
#' @param p1 A `data.frame` with at least columns `x` and `y`. If it contains
#'   multiple rows, only the first row is used.
#' @param snap Logical; if `TRUE` (default), `p0` and `p1` coordinates are
#' rounded
#' to the nearest integer before running the algorithm. If `FALSE`, coordinates
#'   are truncated via `as.integer()`.
#'
#' @return A `data.frame` with columns `x` and `y` giving the integer grid
#' points visited by the line, including both endpoints, in traversal order from `p0`
#' to `p1`.
#'
#' @details
#' The function:
#' \itemize{
#' \item Validates inputs are data frames with `x`/`y` columns and at least one
#' row.
#' \item Converts endpoints to integer coordinates (optionally rounding first).
#'   \item Applies the classic Bresenham algorithm to enumerate grid points.
#' }
#'
#' This is useful for tracing discrete paths on an integer lattice (e.g.,
#' pixels, tile grids, spatial transcriptomics spot grids).
#'
#' @examples
#' p0 <- data.frame(x = 1.2, y = 2.7)
#' p1 <- data.frame(x = 7.9, y = 5.1)
#'
#' # Default (snap = TRUE): rounds endpoints first
#' bresenham_line(p0, p1)
#'
#' # No rounding (snap = FALSE): integer coercion truncates toward zero
#' bresenham_line(p0, p1, snap = FALSE)
#'
#' @export
bresenham_line <- function(p0, p1, snap = TRUE) {
stopifnot(is.data.frame(p0), is.data.frame(p1))
stopifnot(all(c("x","y") %in% names(p0)), all(c("x","y") %in% names(p1)))
stopifnot(nrow(p0) >= 1, nrow(p1) >= 1)

x0 <- p0$x[1]; y0 <- p0$y[1]
x1 <- p1$x[1]; y1 <- p1$y[1]

if (snap) {
    x0 <- round(x0); y0 <- round(y0)
    x1 <- round(x1); y1 <- round(y1)
}

x0 <- as.integer(x0); y0 <- as.integer(y0)
x1 <- as.integer(x1); y1 <- as.integer(y1)

dx <- abs(x1 - x0); dy <- abs(y1 - y0)
sx <- if (x0 < x1) 1L else -1L
sy <- if (y0 < y1) 1L else -1L
err <- dx - dy

outX <- integer(0); outY <- integer(0)

repeat {
    outX <- c(outX, x0)
    outY <- c(outY, y0)
    if (x0 == x1 && y0 == y1) break

    e2 <- 2L * err
    if (e2 > -dy) { err <- err - dy; x0 <- x0 + sx }
    if (e2 <  dx) { err <- err + dx; y0 <- y0 + sy }
}

data.frame(x = outX, y = outY)
}

#' Point-to-segment distance (vectorized)
#'
#' Computes the Euclidean distance from one or more points `(px, py)` to the
#' line segment defined by endpoints `A(ax, ay)` and `B(bx, by)`. The projection
#' is clamped to the segment (i.e., `pos_on_seg` is restricted to `[0, 1]`).
#'
#' @param px,py Numeric vectors (or scalars) giving the x/y coordinates of the
#'   query point(s).
#' @param ax,ay Numeric scalars giving the x/y coordinates of endpoint A.
#' @param bx,by Numeric scalars giving the x/y coordinates of endpoint B.
#'
#' @return A numeric vector of distances, with length equal to
#'   `max(length(px), length(py))` (after R's usual vector recycling rules).
#'
#' @details
#' Let `v = B - A` and `pos_on_seg = ((P - A) · v) / ||v||^2`. The closest point
#' on the
#' infinite line is `A + pos_on_seg v`; clamping `pos_on_seg` to `[0, 1]` yields
#' the closest point
#' on the segment.
#'
#' If `A` and `B` are identical (`||v||^2 == 0`), the function returns the
#' distance from `P` to `A`.
#'
#' @examples
#' # Single point to a horizontal segment
#' point_segment_distance_vec(px = 1, py = 2, ax = 0, ay = 0, bx = 3, by = 0)
#'
#' # Multiple points (vectorized)
#' px <- c(0, 1, 2, 3)
#' py <- c(1, 1, 1, 1)
#' point_segment_distance_vec(px, py, ax = 0, ay = 0, bx = 3, by = 0)
#'
#' # Degenerate segment (A == B)
#' point_segment_distance_vec(px = c(0, 1), py = c(0, 1), ax = 0, ay = 0, bx =
#' 0, by = 0)
#'
#' @export
point_segment_distance_vec <- function(px, py, ax, ay, bx, by) {
vx <- bx - ax
vy <- by - ay
denom <- vx*vx + vy*vy

if (denom == 0) {
    return(sqrt((px - ax)^2 + (py - ay)^2))
}

pos_on_seg <- ((px - ax) * vx + (py - ay) * vy) / denom
pos_on_seg <- pmax(0, pmin(1, pos_on_seg))

projx <- ax + pos_on_seg * vx
projy <- ay + pos_on_seg * vy

sqrt((px - projx)^2 + (py - projy)^2)
}

#' Select spots near a segment and order them along the segment
#'
#' Computes the distance from each spot in `df` to the segment `[p0, p1]`,
#' optionally filters by a maximum distance, keeps the `top_n` closest spots,
#' and finally orders the retained spots by their projection position
#' `pos_on_seg`
#' along the segment (from `p0` to `p1`). This is the foundation function for
#' building
#' trajectory lines across spatial spots.
#'
#' @param df A data frame of spots containing at least columns `x` and `y`.
#'   Additional columns are preserved in the output.
#' @param p0 A data frame with columns `x` and `y` defining the first endpoint
#'   of the segment. If it has multiple rows, only the first row is used.
#' @param p1 A data frame with columns `x` and `y` defining the second endpoint
#'   of the segment. If it has multiple rows, only the first row is used.
#' @param top_n Integer; number of closest spots to keep (default `100`).
#'   If `NULL`, no top-N truncation is applied.
#' @param max_dist Numeric; if not `NULL`, only spots with distance to the
#'   segment `<= max_dist` are kept.
#'
#' @return A data frame containing the selected spots, with two extra columns:
#' \describe{
#'   \item{dist_to_seg}{Euclidean distance from the spot to the segment.}
#' \item{pos_on_seg}{Clamped projection parameter in `[0, 1]` indicating
#' position
#'     along the segment (0 at `p0`, 1 at `p1`).}
#' }
#' The rows are ordered by `pos_on_seg` (i.e., from `p0` to `p1`).
#'
#' @details
#' Distances are computed using \code{point_segment_distance_vec()}.
#' Ordering uses the projection parameter
#' \eqn{pos_on_seg = ((P-A)\cdot(B-A))/||B-A||^2} clamped to `[0, 1]`.
#'
#' This function uses \code{dplyr} verbs (`mutate`, `filter`, `arrange`,
#' `slice_head`) via the `|>` pipe.
#'
#' @examples
#' # Minimal example with base data.frame inputs
#' df <- data.frame(
#'   x = c(0, 1, 2, 3, 4),
#'   y = c(0, 1, 1, 2, 4),
#'   id = letters[1:5]
#' )
#' p0 <- data.frame(x = 0, y = 0)
#' p1 <- data.frame(x = 4, y = 0)
#'
#' # Keep 3 closest spots (no distance threshold)
#' res <- build_one_trajectory(df, p0, p1, top_n = 3)
#' head(res)
#' @importFrom dplyr mutate filter arrange slice_head
#' @export
build_one_trajectory <- function(df, p0, p1, top_n = 100, max_dist = NULL) {
dist_to_seg <- NULL
stopifnot(all(c("x","y") %in% names(df)))
stopifnot(all(c("x","y") %in% names(p0)), all(c("x","y") %in% names(p1)))

ax <- p0$x[1]; ay <- p0$y[1]
bx <- p1$x[1]; by <- p1$y[1]

# distances
d <- point_segment_distance_vec(df$x, df$y, ax, ay, bx, by)

# Projection parameter pos_on_seg along the segment, used for ordering.
vx <- bx - ax
vy <- by - ay
denom <- vx*vx + vy*vy
if (denom == 0) stop("Invalid segment: p0 and p1 coincide (zero-length segment).")
pos_on_seg <- ((df$x - ax) * vx + (df$y - ay) * vy) / denom
pos_on_seg <- pmax(0, pmin(1, pos_on_seg))

out <- df |> dplyr::mutate(dist_to_seg = d, pos_on_seg = pos_on_seg)

if (!is.null(max_dist)) out <- out |> dplyr::filter(dist_to_seg <= max_dist)
if (!is.null(top_n))   out <- out |> dplyr::arrange(dist_to_seg) |> dplyr::slice_head(n = top_n)
out <- out |> dplyr::arrange(pos_on_seg)

out |> dplyr::mutate(trajectory_id = "main")
}


# unique key for removing (x, y) points
.xy_key <- function(x, y) paste0(round(x), "_", round(y))

#' Remove rows whose rounded (x, y) coordinates have already been used
#'
#' Filters `df` to drop any rows whose rounded coordinate key (as produced by
#' `.xy_key(x, y)`) appears in `used_df`.
#'
#' @param df A data frame containing at least columns `x` and `y`.
#' @param used_df A data frame containing at least columns `x` and `y` defining
#'   the set of already-used points to exclude.
#'
#' @return A filtered data frame with the same columns as `df`, excluding rows
#'   whose rounded `(x, y)` matches any row in `used_df`.
#'
#' @details
#' Matching is performed on rounded coordinates, not exact floating-point
#' values.
#' Internally, keys are computed as `paste0(round(x), "_", round(y))`.
#'
#' This function uses `dplyr::filter()` (and the `|>` pipe), so `dplyr` must be
#' installed and a pipe operator must be available.
#'
#' @examples
#' df <- data.frame(x = c(0.1, 1.2, 1.6, 2.0),
#'                  y = c(0.2, 3.4, 3.4, 9.5),
#'                  id = 1:4)
#' used_df <- data.frame(x = c(1.49, 2.1),
#'                       y = c(3.49, 9.49))
#'
#' # Removes rows rounding to (1,3) and (2,9)
#' remove_used_points(df, used_df)
#'
#' @importFrom dplyr mutate filter
#' @export
remove_used_points <- function(df, used_df) {
x <- y <- NULL 
used_keys <- .xy_key(used_df$x, used_df$y)
df |> dplyr::filter(!(.xy_key(x, y) %in% used_keys))
}

#' Find the closest spot to a target point
#'
#' Returns the row of `df` corresponding to the spot nearest to the target
#' coordinates `(tx, ty)` using squared Euclidean distance.
#'
#' @param df A data frame containing at least numeric columns `x` and `y`
#'   representing spot coordinates.
#' @param tx Numeric scalar; target x-coordinate.
#' @param ty Numeric scalar; target y-coordinate.
#'
#' @return A one-row data frame (same columns as `df`) corresponding to the
#' closest spot. The result is always a data frame (`drop = FALSE`).
#'
#' @details
#' The function minimizes `(x - tx)^2 + (y - ty)^2`. Squared distances are used
#' for efficiency; taking the square root is unnecessary for argmin.
#'
#' If multiple spots are tied for the minimum distance, the first occurrence is
#' returned (as per `which.min()`).
#'
#' @examples
#' df <- data.frame(x = c(0, 2, 5), y = c(0, 2, 1), id = c("a", "b", "c"))
#' closest_spot(df, tx = 1.9, ty = 2.1)
#'
#' @export
closest_spot <- function(df, tx, ty) {
d2 <- (df$x - tx)^2 + (df$y - ty)^2
df[which.min(d2), , drop = FALSE]
}

#' Pick a point adjacent to a selected endpoint, on a given side of a segment
#'
#' Given an endpoint (typically the start or end of a previously selected path),
#' this function computes a target point located one perpendicular step away
#' from the endpoint, on the "left" or "right" side relative to the directed
#' segment \eqn{A \rightarrow B}. It then returns the closest spot to that
#' target
#' among the remaining candidates `df_rest`.
#'
#' @param df_rest A data frame of candidate spots containing at least columns
#'   `x` and `y`.
#' @param endpoint A one-row data frame with columns `x` and `y` representing
#'   the endpoint from which to step sideways. If it has multiple rows, only the
#'   first row is used.
#' @param A A data frame with columns `x` and `y` representing point A defining
#'   the direction of the segment. Only the first row is used.
#' @param B A data frame with columns `x` and `y` representing point B defining
#'   the direction of the segment. Only the first row is used.
#' @param spacing Numeric scalar; step size (in the same coordinate units as
#'   `x`/`y`) used to move perpendicularly from `endpoint`.
#' @param side Character; which side to step to relative to the vector
#'   \eqn{A \rightarrow B}. Must be `"left"` or `"right"`.
#'
#' @return A one-row data frame (same columns as `df_rest`) corresponding to the
#' spot closest to the computed perpendicular target point.
#'
#' @details
#' The unit direction vector is computed from \eqn{A \rightarrow B} and a unit
#' left normal \eqn{(-vy, vx)/||v||} is derived. The target point is:
#' \deqn{T = endpoint + sign \cdot spacing \cdot n}
#' where `sign = +1` for `"left"` and `-1` for `"right"`.
#'
#' The closest spot is selected with \code{closest_spot()} using squared
#' Euclidean distance.
#'
#' @examples
#' df <- data.frame(x = c(0, 1, 1, 2, 2), y = c(0, 0, 1, 0, 1), id = 1:5)
#' A <- data.frame(x = 0, y = 0)
#' B <- data.frame(x = 2, y = 0)
#' endpoint <- data.frame(x = 1, y = 0)
#'
#' # Step "left" of A->B (here, toward positive y)
#' adjacent_endpoint(df, endpoint, A, B, spacing = 1, side = "left")
#'
#' # Step "right" of A->B (here, toward negative y)
#' adjacent_endpoint(df, endpoint, A, B, spacing = 1, side = "right")
#'
#' @export
adjacent_endpoint <- function(df_rest, endpoint, A, B, spacing, side = c("left",
"right")) {side <- match.arg(side)

# direction A->B
vx <- B$x[1] - A$x[1]
vy <- B$y[1] - A$y[1]
vnorm <- sqrt(vx*vx + vy*vy)
if (vnorm == 0) stop("A & B are identical.")

# normal left
nx <- -vy / vnorm
ny <-  vx / vnorm

sign <- if (side == "left") 1 else -1

# cible = endpoint + 1 pas perpendiculaire
tx <- endpoint$x[1] + sign * spacing * nx
ty <- endpoint$y[1] + sign * spacing * ny

closest_spot(df_rest, tx, ty)
}


#' Estimate spot spacing from nearest-neighbor distances
#'
#' Estimates the typical spacing between spots by computing each spot's
#' nearest-neighbor distance (excluding itself) and returning the median of
#' those distances. For large datasets, the estimate is computed on a random
#' subsample for speed.
#'
#' @param df A data frame containing at least columns `x` and `y` (numeric)
#'   representing spot coordinates.
#' @param sample_n Integer; maximum number of spots to sample (default `1000`).
#'   If `nrow(df) > sample_n`, a random subset of size `sample_n` is used;
#'   otherwise all spots are used.
#'
#' @return A numeric scalar: the median nearest-neighbor distance (in the same
#' units as `x` and `y`).
#'
#' @details
#' Nearest neighbors are computed with \code{RANN::nn2()} using `k = 2`. The
#' first neighbor is the point itself (distance 0), so the function uses the
#' second neighbor distance `nn$nn.dists[, 2]` as the true nearest-neighbor
#' distance.
#'
#' Missing values are handled via `median(..., na.rm = TRUE)`.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   x = rep(1:5, each = 5),
#'   y = rep(1:5, times = 5)
#' )
#' estimate_spot_spacing(df)
#'
#' @importFrom RANN nn2
#' @importFrom stats median
#' @export
estimate_spot_spacing <- function(df, sample_n = 1000) {
stopifnot(all(c("x","y") %in% names(df)))
n <- nrow(df)

if (n > sample_n) {
    idx <- sample.int(n, sample_n)
    xy <- as.matrix(df[idx, c("x","y")])
} else {
    xy <- as.matrix(df[, c("x","y")])
}

nn <- RANN::nn2(xy, xy, k = 2)

return(stats::median(nn$nn.dists[, 2], na.rm = TRUE))

}



#' Compute the left unit normal of the directed segment A->B
#'
#' Returns the unit-length normal vector pointing to the left of the directed
#' segment \eqn{A \rightarrow B}. For a direction vector \eqn{v = (vx, vy)},
#' the left normal is \eqn{(-vy, vx)} normalized by \eqn{||v||}.
#'
#' @param A A data frame with columns `x` and `y` representing point A.
#'   If it contains multiple rows, only the first row is used.
#' @param B A data frame with columns `x` and `y` representing point B.
#'   If it contains multiple rows, only the first row is used.
#'
#' @return A numeric vector of length 2 with named components:
#' \describe{
#'   \item{nx}{x-component of the left unit normal.}
#'   \item{ny}{y-component of the left unit normal.}
#' }
#'
#' @examples
#' # Horizontal segment to the right: (0,0) -> (1,0)
#' #    Left normal points upward: (0, 1)
#' A <- data.frame(x = 0, y = 0)
#' B <- data.frame(x = 1, y = 0)
#' unit_normal_left(A, B)
#' @export
unit_normal_left <- function(A, B) {
vx <- B$x[1] - A$x[1]
vy <- B$y[1] - A$y[1]
vnorm <- sqrt(vx*vx + vy*vy)
if (vnorm == 0) stop("A and B are identical")
c(nx = -vy / vnorm, ny = vx / vnorm)
}


#' Shift a point along a given direction vector
#'
#' Translates a point `P` by an amount `offset` in the direction `(nx, ny)`.
#' The returned point is:
#' \deqn{P' = (P_x + offset \cdot nx,\; P_y + offset \cdot ny).}
#'
#' @param P A data frame with columns `x` and `y` representing the point to
#'   shift. If it contains multiple rows, only the first row is used.
#' @param nx Numeric scalar; x-component of the direction vector.
#' @param ny Numeric scalar; y-component of the direction vector.
#' @param offset Numeric scalar; translation magnitude (in the same units as
#'   `x`/`y`). Positive values move in the `(nx, ny)` direction; negative values
#'   move in the opposite direction.
#'
#' @return A one-row data frame with columns `x` and `y` giving the shifted
#' point.
#'
#' @examples
#' # Shift the point (1, 2) by 3 units in the direction (0, 1)
#' P <- data.frame(x = 1, y = 2)
#' shift_point(P, nx = 0, ny = 1, offset = 3)
#' @export
shift_point <- function(P, nx, ny, offset) {
data.frame(x = P$x[1] + offset * nx,
            y = P$y[1] + offset * ny)
}


#' Build a single line of spots along a segment
#'
#' Selects spots near the segment defined by endpoints `A` and `B` using
#' `build_one_trajectory()`, then orders the selected spots by their projection
#' parameter (from `A` to `B`). The result is returned as a data frame.
#'
#' @param df_rest A data frame of candidate spots containing at least columns
#'   `x` and `y`.
#' @param A A data frame with columns `x` and `y` defining endpoint A of the
#'   segment. If it contains multiple rows, only the first row is used.
#' @param B A data frame with columns `x` and `y` defining endpoint B of the
#'   segment. If it contains multiple rows, only the first row is used.
#' @param top_n Integer; number of closest spots to keep (default `19`).
#'   If `NULL`, no top-N truncation is applied.
#' @param max_dist Numeric; if not `NULL`, only spots with distance to the
#'   segment `<= max_dist` are kept.
#'
#' @return A data frame of selected spots (subset of `df_rest`) ordered along
#' the
#' segment (in increasing `pos_on_seg`). The output includes the extra columns
#' added by
#' `build_one_trajectory()` (typically `dist_to_seg` and `pos_on_seg`).
#'
#' @details
#' This function is a thin wrapper around `build_one_trajectory()` that returns
#' only the ordered `selected` data frame (and not the segment metadata).
#'
#' It uses the base R pipe `|>` and `dplyr::arrange()`.
#'
#' @examples
#' df <- data.frame(
#'   x = c(0, 1, 2, 3, 4, 2),
#'   y = c(0, 0, 0, 0, 0, 1),
#'   id = 1:6
#' )
#' A <- data.frame(x = 0, y = 0)
#' B <- data.frame(x = 4, y = 0)
#'
#' build_one_line(df, A, B, top_n = 5)
#'
#' @export
build_one_line <- function(df_rest, A, B, top_n = 19, max_dist = NULL) {

res <- build_one_trajectory(df_rest, A, B, top_n = top_n, max_dist = max_dist)
res |> dplyr::arrange(pos_on_seg)
}

#' @noRd
.build_at_offset <- function(df_rest, A, B, nx, ny, offset, label,
                            top_n = 19, max_dist = NULL) {
Ash <- shift_point(A, nx, ny, offset)
Bsh <- shift_point(B, nx, ny, offset)

build_one_line(df_rest, Ash, Bsh, top_n = top_n, max_dist = max_dist) |>
    dplyr::mutate(trajectory_id = label, offset = offset)
}

#' Build parallel spot lines around a central segment
#'
#' Constructs a central line of spots near the segment defined by endpoints `A`
#' and `B`, then builds additional parallel lines on the left, right, or both
#' sides by translating the segment along its left unit normal. After each line
#' is built, its spots are removed from the remaining pool to avoid reuse.
#'
#' @param df A data frame of candidate spots containing at least columns `x` and
#'   `y`. Additional columns are preserved.
#' @param A A data frame with columns `x` and `y` defining endpoint A of the
#'   central segment. If it contains multiple rows, only the first row is used.
#' @param B A data frame with columns `x` and `y` defining endpoint B of the
#'   central segment. If it contains multiple rows, only the first row is used.
#' @param top_n Integer; number of closest spots to keep per line (default
#' `19`).
#'   If `NULL`, no top-N truncation is applied in the underlying selection.
#' @param n_extra Integer; number of additional lines to build on each requested
#'   side (default `2`). For example, `n_extra = 2` with `side = "both"` yields
#'   1 center line + 2 left + 2 right.
#' @param side Character; which side(s) to build relative to the vector
#'   \eqn{A \rightarrow B}. One of `"left"`, `"right"`, or `"both"`.
#' @param lane_width_factor Numeric scalar; multiplier applied to the estimated
#'   spot spacing to obtain the inter-line offset (default `1.15`).
#' @param max_dist Numeric; if not `NULL`, only spots with distance to each
#'   (translated) segment `<= max_dist` are considered when building each line.
#'
#' @return A data frame containing all selected spots from all lines,
#' stacked together, with additional columns `trajectory_id` (e.g., `"main"`,
#' `"left_1"`, `"right_1"`) and `offset` (signed offset used for that line).
#'
#' @details
#' The inter-line distance is computed as:
#' \deqn{w = lane\_width\_factor \cdot spacing}
#' where `spacing` is the median nearest-neighbor distance estimated from `df`.
#'
#' Lines are generated by shifting both endpoints `A` and `B` by `offset * n`,
#' where `n` is the left unit normal of \eqn{A \rightarrow B}. Positive offsets
#' correspond to the left side; negative offsets correspond to the right side.
#'
#' This function relies on helpers such as `estimate_spot_spacing()`,
#' `unit_normal_left()`, `shift_point()`, `build_one_line()`, and
#' `remove_used_points()`.
#'
#' The implementation uses the base R pipe `|>` and `dplyr::mutate()` /
#' `dplyr::bind_rows()`. Ensure `dplyr` is installed.
#'
#' @examples
#' df <- data.frame(
#'   x = rep(1:10, each = 3),
#'   y = rep(1:3, times = 10),
#'   id = seq_len(30)
#' )
#' A <- data.frame(x = 1, y = 2)
#' B <- data.frame(x = 10, y = 2)
#'
#' out <- build_similar_trajectories(df, A, B, top_n = 5, n_extra = 1, side =
#' "both")
#' head(out)
#' unique(out$trajectory_id)
#'
#' @export
#' @importFrom dplyr mutate arrange bind_rows
build_similar_trajectories <- function(df, A, B,
                                top_n = 19,
                                n_extra = 2,
                                side = c("left","right","both"),
                                lane_width_factor = 1.15,
                                max_dist = NULL) {

side <- match.arg(side)
spacing <- estimate_spot_spacing(df)
w <- lane_width_factor * spacing

nvec <- unit_normal_left(A, B)
nx <- nvec["nx"]; ny <- nvec["ny"]

lines <- list()

# center line (main trajectory)
sel0 <- build_one_line(df, A, B, top_n = top_n, max_dist = max_dist) |>
    dplyr::mutate(trajectory_id = "main", offset = 0)

lines[[1]] <- sel0
df_rest <- remove_used_points(df, sel0)

# build successively (and remove used points each time)
if (side %in% c("left","both")) {
    for (k in seq_len(n_extra)) {
    offset <-  k * w
    selk <- .build_at_offset(df_rest, A, B, nx, ny, offset,
                            label = paste0("left_", k),
                            top_n = top_n, max_dist = max_dist)
    lines[[length(lines) + 1]] <- selk
    df_rest <- remove_used_points(df_rest, selk)
    }
}

if (side %in% c("right","both")) {
    for (k in seq_len(n_extra)) {
    offset <- -k * w
    selk <- .build_at_offset(df_rest, A, B, nx, ny, offset,
                            label = paste0("right_", k),
                            top_n = top_n, max_dist = max_dist)
    lines[[length(lines) + 1]] <- selk
    df_rest <- remove_used_points(df_rest, selk)
    }
}

dplyr::bind_rows(lines)
}

#' Filter lines by endpoint cluster membership
#'
#' Filters trajectory data by checking the cluster labels at each trajectory's
#' endpoints.
#' For every `trajectory_id`, the start endpoint is defined as the spot with the
#' smallest projection parameter `pos_on_seg`, and the end endpoint as the spot
#' with the largest `pos_on_seg`. Only trajectories whose start cluster is in
#' `allowed_start_clusters`
#' *and* whose end cluster is in `allowed_end_clusters` are kept.
#'
#' @param out A data frame containing selected spots with columns
#' `trajectory_id`,
#' `cluster`, and `pos_on_seg`, typically the output of
#' `build_similar_trajectories()`.
#' @param allowed_start_clusters Vector of allowed cluster labels for the start
#'   endpoint.
#' @param allowed_end_clusters Vector of allowed cluster labels for the end
#'   endpoint.
#'
#' @return The same data frame as `out`, filtered to keep only the
#' lines matching the allowed endpoint cluster constraints.
#'
#' @details
#' The endpoint clusters are computed per `trajectory_id`:
#' \itemize{
#'   \item `start_cluster = cluster[which.min(pos_on_seg)]`
#'   \item `end_cluster   = cluster[which.max(pos_on_seg)]`
#' }
#' If multiple spots share the same minimum/maximum `pos_on_seg`, the first is
#' taken (as per `which.min()` / `which.max()`).
#'
#' This function uses `dplyr` (`group_by`, `summarise`, `filter`, `pull`) and
#' the  base R pipe `|>`.
#'
#' @examples
#' # Minimal example
#' out <- data.frame(
#'   trajectory_id = c("L1","L1","L2","L2"),
#'   pos_on_seg = c(0.0, 1.0, 0.0, 1.0),
#'   cluster  = c("A", "B", "A", "C"),
#'   x = 1:4, y = 1:4
#' )
#'
#' # Keep only trajectories starting in A and ending in B
#' filter_out_by_endpoint_clusters(out, allowed_start_clusters = "A",
#' allowed_end_clusters = "B")
#'
#' @export
#' @importFrom dplyr group_by summarise filter
filter_out_by_endpoint_clusters <- function(out,
                                        allowed_start_clusters,
                                        allowed_end_clusters) {
trajectory_id <- cluster <- start_cluster <- end_cluster <- NULL                                       
stopifnot(is.data.frame(out))
stopifnot("trajectory_id" %in% names(out))
stopifnot("cluster" %in% names(out))
stopifnot("pos_on_seg" %in% names(out))  # provided by build_one_trajectory() and related helpers
ends <- out |>
    dplyr::group_by(trajectory_id) |>
    dplyr::summarise(
    start_cluster = cluster[which.min(pos_on_seg)],
    end_cluster   = cluster[which.max(pos_on_seg)],
    .groups = "drop"
    )

keep_ids <- ends |>
    dplyr::filter(
    start_cluster %in% allowed_start_clusters,
    end_cluster   %in% allowed_end_clusters
    ) |>
    dplyr::pull(trajectory_id)

out |>
    dplyr::filter(trajectory_id %in% keep_ids)
}