viewer_spe <- function(img = FALSE) {
    n <- 6
    spe <- SpatialExperiment::SpatialExperiment(
        assays = list(counts = matrix(1, 1, n,
            dimnames = list("g", paste0("s", seq_len(n))))),
        spatialCoords = cbind(x = seq_len(n) * 10.123456789, y = n:1),
        sample_id = "S1")
    spe$cluster <- factor(c("a", "a", "b", "b", NA, "solo"))
    spe$score <- c(1:5 / 3, NA)
    if (img) {
        f <- tempfile(fileext = ".png")
        magick::image_write(magick::image_blank(40, 20, "pink"), f)
        spe <- SpatialExperiment::addImg(spe, sample_id = "S1",
            image_id = "lowres", imageSource = f, scaleFactor = 0.5,
            load = TRUE)
    }
    spe
}
skip_viewer <- function() {
    skip_if_not_installed("plotly")
    skip_if_not_installed("htmlwidgets")
}
# traces and layout as sent to the browser (first trace is empty)
built <- function(p) plotly::plotly_build(p)$x
ids_of <- function(p) p$jsHooks$render[[1]]$data

test_that("view_spots keeps ids aligned with trace points", {
    skip_viewer()
    spe <- viewer_spe()
    for (col in list(NULL, "cluster", "score")) {
        p <- view_spots(spe, color_by = col)
        ids <- ids_of(p)
        tr <- utils::tail(built(p)$data, length(ids))
        expect_setequal(unlist(ids), colnames(spe))
        xy <- SpatialExperiment::spatialCoords(spe)
        rownames(xy) <- colnames(spe)
        for (i in seq_along(ids)) {
            expect_equal(as.numeric(tr[[i]]$x),
                unname(xy[as.character(ids[[i]]), "x"]))
            expect_equal(as.numeric(tr[[i]]$y),
                unname(xy[as.character(ids[[i]]), "y"]))
        }
    }
    expect_length(ids_of(view_spots(spe)), 1)
    expect_length(ids_of(view_spots(spe, "cluster")), 4)
    expect_length(ids_of(view_spots(spe, "score")), 2)
})

test_that("view_spots places the slide image in spatialCoords units", {
    skip_viewer()
    skip_if_not_installed("magick")
    skip_if_not_installed("base64enc")
    spe <- viewer_spe(img = TRUE)
    lay <- built(view_spots(spe))$layout
    expect_equal(lay$images[[1]]$sizex, 80)
    expect_equal(lay$images[[1]]$sizey, 40)
    expect_equal(lay$yaxis$range, c(40, 0))
    expect_equal(built(view_spots(spe, scale_factor = 1))$layout$images[[1]]$sizex, 40)
    expect_length(built(view_spots(spe, image = FALSE))$layout$images, 0)
    expect_error(view_spots(spe, flip_y = FALSE), "misalign")
    expect_error(view_spots(spe, image_id = "nope"), "No image")
    expect_error(view_spots(spe, image_file = "nope.png"), "not found")

    # an unreadable stored image is not fatal
    f <- tempfile(fileext = ".png")
    file.create(f)
    bad <- SpatialExperiment::addImg(viewer_spe(), sample_id = "S1",
        image_id = "x", imageSource = f, scaleFactor = 1, load = FALSE)
    expect_warning(p <- view_spots(bad), "not drawn")
    expect_length(built(p)$layout$images, 0)
})

test_that("view_spots checks its inputs", {
    skip_viewer()
    spe <- viewer_spe()
    expect_error(view_spots(data.frame()))
    expect_error(view_spots(spe, color_by = "nope"), "color_by")
    expect_error(view_spots(spe, sample_id = "nope"), "sample_id")
    two <- SpatialExperiment::SpatialExperiment(
        assays = list(counts = matrix(1, 1, 6,
            dimnames = list("g", paste0("s", 1:6)))),
        spatialCoords = cbind(x = 1:6, y = 6:1),
        sample_id = rep(c("S1", "S2"), each = 3))
    expect_error(view_spots(two), "Several samples")
    expect_setequal(unlist(ids_of(view_spots(two, sample_id = "S2"))),
        c("s4", "s5", "s6"))
})
