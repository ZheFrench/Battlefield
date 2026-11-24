# Local quadratic smoother matrix (tricube weights), as loess(degree = 2).
# The window always spans at least 3 distinct positions.
.trend_smoother <- function(pos, grid, span) {
    q <- max(1L, floor(span * length(pos)))
    distinct <- unique(round(pos, 8))
    t(vapply(grid, function(g) {
        u <- pos - g
        h <- max(sort(abs(u))[q], sort(abs(distinct - g))[5])
        w <- pmax(1 - abs(u / h)^3, 0)^3
        X <- cbind(1, u / h, (u / h)^2)
        solve(crossprod(X, w * X), t(w * X))[1, ]
    }, numeric(length(pos))))
}

# Slope of each row of Y along pos, adjusted for line, by generalized least
# squares. Residual covariance is exponential plus nugget; range and nugget are
# chosen per gene by REML over a grid shared by all genes.
.trend_gls <- function(Y, pos, line, coords) {
    X <- model.matrix(if (nlevels(line) > 1L) ~ line + pos else ~pos)
    q <- ncol(X)
    dof <- length(pos) - q
    D <- as.matrix(dist(coords))
    if (!any(D > 0)) stop("All trajectory spots share the same coordinates.")
    ranges <- exp(seq(log(min(D[D > 0]) / 2), log(max(D)), length.out = 8))
    nuggets <- c(0.02, 0.1, 0.25, 0.5, 0.75, 0.95)
    best <- rep(-Inf, nrow(Y))
    slope <- se <- rep(NA_real_, nrow(Y))
    for (r in ranges) {
        for (nugget in nuggets) {
            U <- chol((1 - nugget) * exp(-D / r) + diag(nugget, length(pos)))
            Xs <- backsolve(U, X, transpose = TRUE)
            Ys <- backsolve(U, t(Y), transpose = TRUE)
            R <- chol(crossprod(Xs))
            XtX_inv <- chol2inv(R)
            beta <- XtX_inv %*% crossprod(Xs, Ys)
            rss <- pmax(colSums((Ys - Xs %*% beta)^2), .Machine$double.eps)
            reml <- -sum(log(diag(U))) - sum(log(diag(R))) - dof / 2 * log(rss)
            k <- reml > best
            best[k] <- reml[k]
            slope[k] <- beta[q, k]
            se[k] <- sqrt(rss[k] / dof * XtX_inv[q, q])
        }
    }
    t_stat <- slope / se
    list(slope = slope, p_down = pt(t_stat, dof),
        p_up = pt(t_stat, dof, lower.tail = FALSE))
}

#' Compare gene trends along a fixed trajectory
#'
#' Tests whether genes decrease or increase along trajectories already stored
#' by [add_trajectories_to_spe()], and reports how similar their smoothed
#' profile is to that of a reference gene. Spots are
#' fixed; only the genes change.
#'
#' @param spe A SpatialExperiment with `trajectory_id` and `pos_on_seg` in its
#'   colData.
#' @param genes Character vector of genes to test.
#' @param reference Reference gene. It is reported but never tested.
#' @param assay_name Assay holding normalized expression.
#' @param feature_names Optional rowData column with gene names. Row names are
#'   used if `NULL`.
#' @param trajectory_ids Trajectories to pool. All are used if `NULL`.
#' @param min_detected Minimal fraction of trajectory spots expressing a gene.
#' @param span Fraction of spots used by the local smoother, in (0, 1].
#' @param fdr_cutoff FDR threshold used for `trend`.
#'
#' @details
#' `slope` is the estimated change in expression from the start to the end of
#' the trajectory (`pos_on_seg` from 0 to 1). It comes from a generalized least
#' squares model adjusted for the trajectory line. Residuals of neighbouring
#' spots are allowed to be correlated (exponential covariance on the spatial
#' coordinates, plus a nugget), with parameters chosen per gene by REML.
#' P-values are one-sided and adjusted by Benjamini-Hochberg over tested genes.
#' Smoothed profiles pool the spots of all lines. `r_curve` is the Pearson
#' correlation of each profile with that of the reference; it is descriptive
#' and has no p-value.
#' A spatial pattern about as long as the trajectory cannot be told apart from
#' a trend, so the FDR remain optimistic in that case: use them to rank genes
#' within a sample, and replicate across samples to conclude.
#' `trend` is "dec", "inc", "unknown" or "ref" (`NA` if not tested).
#' Do not use the gene that guided the choice of trajectory as evidence.
#'
#' @return A data.frame with one row per gene: the reference first, then the
#'   tested genes by increasing slope, then the genes that could not be tested
#'   (statistics are `NA`), with columns `gene`, `slope`, `r_curve`, `fdr_down`,
#'   `fdr_up`, `trend` and `status`. `status` is "reference", "tested", "not
#'   found",
#'   "duplicated", "missing values" or "too sparse or constant". Smoothed
#'   profiles, for plotting, are in `attr(result, "curves")`.
#'
#' @examples
#' data("visiumHD_16um_simulated_spe", package = "Battlefield")
#' spe <- visiumHD_16um_simulated_spe
#' df <- data.frame(
#'   spot_id = colnames(spe),
#'   x = SpatialExperiment::spatialCoords(spe)[, 1],
#'   y = SpatialExperiment::spatialCoords(spe)[, 2],
#'   cluster = SummarizedExperiment::colData(spe)$cluster
#' )
#' centroids <- compute_centroids(df)
#' A <- centroids[centroids$cluster == 3, c("x", "y")]
#' B <- centroids[centroids$cluster == 4, c("x", "y")]
#' trajs <- build_similar_trajectories(df, A, B, top_n = 12, n_extra = 2,
#'   side = "both")
#'
#' # The example data hold a single gene: simulate genes along the trajectory
#' set.seed(1)
#' pos <- rep(0.5, nrow(df))
#' pos[match(trajs$spot_id, df$spot_id)] <- trajs$pos_on_seg
#' sim <- function(slope) 3 + slope * pos + stats::rnorm(length(pos), sd = 0.3)
#' logcounts <- rbind(REF = sim(-2), DOWN = sim(-1.5), UP = sim(1.5),
#'   FLAT = sim(0))
#' colnames(logcounts) <- df$spot_id
#' toy <- SpatialExperiment::SpatialExperiment(
#'   assays = list(logcounts = logcounts),
#'   spatialCoords = as.matrix(df[, c("x", "y")])
#' )
#' toy <- add_trajectories_to_spe(toy, trajectory = trajs)
#'
#' compare_gene_trends(toy, genes = c("DOWN", "UP", "FLAT"),
#'   reference = "REF")
#'
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom stats cor dist model.matrix p.adjust pt
#' @export
compare_gene_trends <- function(spe, genes, reference,
    assay_name = "logcounts", feature_names = NULL, trajectory_ids = NULL,
    min_detected = 0.1, span = 0.8, fdr_cutoff = 0.05) {
    stopifnot(is(spe, "SpatialExperiment"), is.character(genes),
        length(reference) == 1L, span > 0, span <= 1)

    # Trajectory spots
    cd <- colData(spe)
    if (!all(c("trajectory_id", "pos_on_seg") %in% colnames(cd))) {
        stop("Run add_trajectories_to_spe() first.")
    }
    keep <- !is.na(cd$trajectory_id) & is.finite(cd$pos_on_seg)
    if (!is.null(trajectory_ids)) {
        keep <- keep & cd$trajectory_id %in% trajectory_ids
    }
    pos <- cd$pos_on_seg[keep]
    line <- factor(cd$trajectory_id[keep])
    if (length(pos) < 20L) stop("Fewer than 20 trajectory spots.")
    if (length(unique(round(pos, 8))) < 5L) {
        stop("Fewer than 5 distinct positions along the trajectory.")
    }

    # Expression of the genes matching exactly one row
    if (is.null(feature_names)) {
        row_genes <- rownames(spe)
    } else if (feature_names %in% colnames(rowData(spe))) {
        row_genes <- as.character(rowData(spe)[[feature_names]])
    } else {
        stop("'feature_names' is not a rowData column.")
    }
    genes <- unique(c(reference, genes))
    hits <- vapply(genes, function(g) sum(row_genes == g, na.rm = TRUE), 0L)
    if (hits[[1]] != 1L) stop("'reference' must match exactly one gene.")
    found <- genes[hits == 1L]
    Y <- as.matrix(assay(spe, assay_name)[
        match(found, row_genes), keep, drop = FALSE])
    rownames(Y) <- found

    # Status of each gene
    has_na <- rowSums(!is.finite(Y)) > 0
    constant <- apply(Y, 1, function(y) length(unique(y)) == 1L)
    if (has_na[[1]]) stop("'reference' has missing values.")
    if (constant[[1]]) stop("'reference' is constant.")
    sparse <- rowMeans(Y > 0) < min_detected
    status <- ifelse(hits == 0L, "not found", "duplicated")
    status[found] <- ifelse(has_na, "missing values",
        ifelse(sparse | constant, "too sparse or constant", "tested"))
    status[1] <- "reference"
    kept <- genes[status %in% c("reference", "tested")]
    if (length(kept) < 2L) stop("No testable gene.")
    Y <- Y[kept, , drop = FALSE]

    # Slopes and smoothed profiles
    fit <- .trend_gls(Y, pos, line,
        SpatialExperiment::spatialCoords(spe)[keep, , drop = FALSE])
    grid <- seq(min(pos), max(pos), length.out = 40)
    curve <- tcrossprod(.trend_smoother(pos, grid, span), Y)

    # The reference is not tested: its p-values are left out of the adjustment
    tested <- data.frame(gene = kept, slope = fit$slope,
        r_curve = drop(cor(curve, curve[, 1])),
        fdr_down = p.adjust(replace(fit$p_down, 1, NA), "BH"),
        fdr_up = p.adjust(replace(fit$p_up, 1, NA), "BH"))
    tested$trend <- ifelse(tested$fdr_down < fdr_cutoff, "dec",
        ifelse(tested$fdr_up < fdr_cutoff, "inc", "unknown"))
    tested$trend[1] <- "ref"

    res <- tested[match(genes, kept), ]
    res$gene <- genes
    res$status <- unname(status)
    res <- res[order(res$status != "reference", res$slope), ]
    rownames(res) <- NULL
    attr(res, "curves") <- data.frame(gene = rep(kept, each = length(grid)),
        position = rep(grid, length(kept)),
        mean = as.vector(curve))
    res
}
