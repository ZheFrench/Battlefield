#' Detect the grid type (square vs hexagonal) from spatial coordinates
#'
#' This function inspects k-nearest-neighbor distances in a 2D coordinate set to
#' infer whether the points lie on a **square** or **hexagonal** lattice.
#'
#' The method:
#' \enumerate{
#'   \item Computes the nearest-neighbor distance for each point and uses the
#'   median of those distances as an estimate of the grid step size.
#'   \item Normalizes all neighbor distances by this step size.
#'   \item Looks for a characteristic secondary distance peak:
#'   \itemize{
#'     \item square grid: ratio \eqn{\approx \sqrt{2}}
#'     \item hexagonal grid: ratio \eqn{\approx \sqrt{3}}
#'   }
#' }
#'
#' @param df A data.frame containing spatial coordinates.
#' @param coords Character vector of length 2 giving the coordinate column
#' names.
#'   Default is `c("x","y")`.
#' @param k Integer. Number of neighbors to consider (internally uses `k + 1`
#'   to include the point itself, then removes it). Default is 12.
#' @param tolerance Numeric. Tolerance for matching the characteristic ratio to
#'   \eqn{\sqrt{2}} or \eqn{\sqrt{3}}. Default is 0.1.
#' @param verbose Logical. If `TRUE`, prints diagnostic messages. Default is
#' `TRUE`.
#'
#' @details
#' The grid step size is estimated as the median of the per-point
#' nearest-neighbor
#' distance. The "median distance ratio" is computed from normalized neighbor
#' distances restricted to the interval (1.2, 2.0), which tends to capture the
#' second-shell distances on regular grids.
#'
#' @return A list with:
#' \describe{
#'   \item{grid_type}{One of `"square"`, `"hexagonal"`, or `"unknown"`.}
#'   \item{step}{Estimated grid step size (median nearest-neighbor distance).}
#' \item{median_ratio}{Median normalized distance ratio in the (1.2, 2.0)
#' range.}
#' }
#'
#' @examples
#' # --- Example 1: square grid ---
#' sq <- expand.grid(x = 0:9, y = 0:9)
#' res_sq <- detect_grid_type(sq, coords = c("x","y"), k = 12, tolerance = 0.1)
#' res_sq$grid_type
#'
#' # --- Example 2: hexagonal grid (pointy-top axial-like layout) ---
#' nx <- 10; ny <- 10
#' hex <- expand.grid(i = 0:(nx-1), j = 0:(ny-1))
#' hex$x <- hex$i + 0.5 * (hex$j %% 2)
#' hex$y <- (sqrt(3)/2) * hex$j
#' hex <- hex[, c("x","y")]
#' res_hex <- detect_grid_type(hex, coords = c("x","y"), k = 12, tolerance =
#' 0.1)
#' res_hex$grid_type
#'
#' @importFrom RANN nn2
#' @export
detect_grid_type <- function(df,
                            coords = c("x", "y"),
                            k = 12,
                            tolerance = 0.1,
                            verbose = TRUE) {
stopifnot(all(coords %in% colnames(df)))

X <- as.matrix(df[, coords])

nn <- RANN::nn2(X, X, k = k + 1)
dists <- nn$nn.dists[, -1, drop = FALSE]

# Flatten distances and keep strictly positive ones
v <- as.vector(dists)
v <- v[v > 0]

# Per-point minimum distance (nearest neighbor) => estimate grid step
d0 <- apply(dists, 1, min)
step <- median(d0)

# Normalize distances by step size
ratios <- v / step

# Look for the median ratio near sqrt(2) or sqrt(3)
med_ratio <- median(ratios[ratios > 1.2 & ratios < 2.0])

grid_type <- if (abs(med_ratio - sqrt(2)) < tolerance) {
    "square"
} else if (abs(med_ratio - sqrt(3)) < tolerance) {
    "hexagonal"
} else {
    "unknown"
}

if (verbose) {
    message("---- Grid type detection ----")
    message("Estimated grid step: ", signif(step, 4))
    message("Median distance ratio: ", signif(med_ratio, 4))
    message("Detected grid type: ", grid_type)
    message("--------------------------------") 

}

invisible(list(
    grid_type = grid_type,
    step = step,
    median_ratio = med_ratio
))
}

#' Suggest neighborhood parameters based on detected grid type
#'
#' This helper detects the spatial grid type (square vs hexagonal) and returns
#' recommended neighborhood parameters for building local adjacency / neighbor
#' queries (e.g., for border detection or spatial graphs).
#'
#' It relies on [detect_grid_type()] to estimate the grid step size and grid
#' geometry, then chooses typical parameters:
#' \itemize{
#'   \item **Hexagonal (Visium)**: 6-neighborhood, radius ~ step
#'   \item **Square (Visium HD)**: 4- or 8-neighborhood, radius ~ step (4) or
#'   ~ sqrt(2)*step (8)
#' }
#'
#' @param df A data.frame containing spatial coordinates.
#' @param coords Character vector of length 2 giving the coordinate column
#' names.
#'   Default is `c("x","y")`.
#' @param square_connectivity Integer-like choice for square grids: `4` (Von
#' Neumann)
#'   or `8` (Moore). Default is `c(4, 8)` which selects the first value (`4`).
#' @param tolerance Numeric. Passed to [detect_grid_type()] to match
#' characteristic
#'   ratios (sqrt(2) / sqrt(3)). Default is 0.1.
#' @param verbose Logical. If `TRUE`, prints diagnostic messages. Default is
#' `TRUE`.
#'
#' @details
#' Returned values:
#' \describe{
#' \item{grid_type}{Detected grid type: `"hexagonal"`, `"square"`, or
#' `"unknown"`.}
#' \item{connectivity}{Nominal grid connectivity (6 for hex, 4 or 8 for
#' square).}
#'   \item{r}{Recommended distance threshold (radius) to define adjacency.}
#'   \item{k_max}{A suggested upper bound for kNN queries (used as a safe cap).}
#'   \item{comment}{Human-readable description of the choice.}
#' }
#'
#' Notes:
#' \itemize{
#' \item `radius` is set to `1.01 * step` (or `1.01 * sqrt(2) * step` for
#' #' 8-connectivity)
#'   to be slightly permissive under small numerical noise.
#'   \item `k` is deliberately larger than the nominal connectivity to make sure
#'   enough candidates are retrieved before applying distance-based filtering.
#' }
#'
#' @return A named list with neighborhood parameters: `grid_type`,
#'`connectivity`,
#' `radius`, `k`, and `comment`.
#'
#' @examples
#' # Square grid example
#' sq <- expand.grid(x = 0:9, y = 0:9)
#' get_neighborhood_params(sq, square_connectivity = 4, verbose = FALSE)
#'
#' # Hexagonal grid example
#' nx <- 10; ny <- 10
#' hex <- expand.grid(i = 0:(nx-1), j = 0:(ny-1))
#' hex$x <- hex$i + 0.5 * (hex$j %% 2)
#' hex$y <- (sqrt(3)/2) * hex$j
#' hex <- hex[, c("x","y")]
#' get_neighborhood_params(hex, verbose = FALSE)
#'
#' @export
get_neighborhood_params <- function(df,
                                    coords = c("x", "y"),
                                    square_connectivity = c(4, 8),
                                    tolerance = 0.1,
                                    verbose = TRUE) {

square_connectivity <- match.arg(as.character(square_connectivity),
                                choices = c("4", "8"))
square_connectivity <- as.integer(square_connectivity)

# 1) Detect grid type
det <- detect_grid_type(df,
                        coords = coords,
                        tolerance = tolerance,
                        verbose = verbose)

step <- det$step
grid <- det$grid_type

# 2) Default parameters based on grid geometry
if (grid == "hexagonal") {

    params <- list(
    grid_type    = "hexagonal",
    connectivity = 6,
    radius       = 1.01 * step,
    k           = 6,
    comment      = "Standard Visium: 6 hexagonal neighbors"
    )

} else if (grid == "square") {

    if (square_connectivity == 4) {
    radius <- 1.01 * step
    k <- 16
    conn <- 4
    } else {
    radius <- 1.01 * sqrt(2) * step
    k <- 32
    conn <- 8
    }

    params <- list(
    grid_type    = "square",
    connectivity = conn,
    radius       = radius,
    k            = k,
    comment      = paste0("Visium HD: ", conn, "-connectivity (square grid)")
    )

} else {

    # Conservative fallback
    params <- list(
    grid_type    = "unknown",
    connectivity = NA,
    radius       = 1.01 * step,
    k        = 32,
    comment      = "Grid not detected: conservative parameters"
    )
}

if (verbose) {
    message("---- Neighborhood parameters ----")
    for (n in names(params)) {
    message(n, ": ", params[[n]])
    }
    message("--------------------------------") 
}

invisible(params)
}