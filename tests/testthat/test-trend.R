trend_spe <- function(m, pos, line) {
    colnames(m) <- paste0("s", seq_along(pos))
    spe <- SpatialExperiment::SpatialExperiment(
        assays = list(logcounts = m),
        spatialCoords = cbind(x = pos * 39, y = as.numeric(factor(line))))
    spe$trajectory_id <- line
    spe$pos_on_seg <- pos
    spe
}
pos <- rep(seq(0, 1, length.out = 40), 3)
line <- rep(c("main", "left_1", "right_1"), each = 40)
sim <- function(slope, at = pos) 5 + slope * at + rnorm(length(at), sd = 0.3)

test_that("compare_gene_trends recovers direction and reports bad genes", {
    set.seed(1)
    m <- rbind(REF = sim(-2), DOWN = sim(-1.5), UP = sim(1.5), FLAT = sim(0),
        RARE = c(1, rep(0, 119)), GAP = c(NA, sim(0)[-1]), DUP = sim(0),
        DUP = sim(0))
    genes <- c("FLAT", "UP", "DOWN", "RARE", "NOPE", "GAP", "DUP")
    res <- compare_gene_trends(trend_spe(m, pos, line), genes, "REF")

    expect_equal(colnames(res), c("gene", "slope", "r_curve", "fdr_down",
        "fdr_up", "trend", "status"))
    expect_equal(res$gene,
        c("REF", "DOWN", "FLAT", "UP", "RARE", "NOPE", "GAP", "DUP"))
    expect_equal(res$trend, c("ref", "dec", "unknown", "inc", rep(NA, 4)))
    expect_equal(res$status, c("reference", rep("tested", 3),
        "too sparse or constant", "not found", "missing values", "duplicated"))
    expect_true(all(is.na(res$slope[5:8])))
    expect_gt(res$r_curve[2], 0.8)
    expect_lt(res$r_curve[4], 0)

    curves <- attr(res, "curves")
    expect_equal(colnames(curves), c("gene", "position", "mean"))
    expect_equal(nrow(curves), 4 * 40)
})

test_that("compare_gene_trends rejects invalid input", {
    set.seed(1)
    spe <- trend_spe(rbind(REF = sim(-2), DOWN = sim(-2)), pos, line)
    expect_error(compare_gene_trends(spe, "DOWN", "NOPE"), "reference")
    expect_error(compare_gene_trends(spe, "DOWN", "REF", span = 1.5), "span")
    expect_error(compare_gene_trends(spe, "DOWN", "REF",
        feature_names = "zz"), "feature_names")
    expect_error(compare_gene_trends(spe, "DOWN", "REF",
        trajectory_ids = "zz"), "20 trajectory spots")
    expect_error(compare_gene_trends(spe, "NOPE", "REF"), "testable")
})

test_that("compare_gene_trends supports one line and tied positions", {
    set.seed(1)
    tied <- c(rep(0, 25), seq(0, 1, length.out = 10), rep(1, 25))
    spe <- trend_spe(rbind(REF = sim(-2, tied), DOWN = sim(-2, tied)), tied,
        rep("main", 60))
    res <- compare_gene_trends(spe, "DOWN", "REF")
    expect_true(all(is.finite(res$slope)) && all(is.finite(res$r_curve)))
})

test_that("trend smoother handles few symmetric positions", {
    x <- rep(seq(0, 1, length.out = 5), 4)
    S <- .trend_smoother(x, seq(0, 1, length.out = 9), 0.2)
    expect_equal(drop(S %*% x^2), seq(0, 1, length.out = 9)^2)
})

test_that("trend GLS agrees with least squares on independent noise", {
    set.seed(4)
    Y <- matrix(rnorm(300 * 120), 300)
    coords <- cbind(pos * 39, rep(1:3, each = 40))
    fit <- .trend_gls(Y, pos, factor(line), coords)
    ols <- apply(Y, 1, function(y) coef(lm(y ~ line + pos))[["pos"]])
    expect_gt(cor(fit$slope, ols), 0.95)
    expect_lt(mean(fit$p_down < 0.05), 0.1)
})
