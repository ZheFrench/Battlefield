# JavaScript run on the plotly widget: click / lasso selection and the panel
# holding the R code to copy. `ids` holds one vector of spot ids per trace.
.viewer_js <- "
function(el, x, ids) {
  var picked = [], n0 = ids.length, first = el.data.length - n0;
  var box = document.createElement('div');
  box.style.cssText = 'position:fixed;left:10px;bottom:10px;z-index:10;' +
    'width:420px;max-width:90vw;box-sizing:border-box;background:#fff;' +
    'color:#222;border:1px solid #999;padding:8px;font:12px sans-serif';
  box.innerHTML = '<b class=\"n\"></b> <button class=\"c\">Copy</button> ' +
    '<button class=\"u\">Undo last</button> ' +
    '<button class=\"z\">Clear</button>' +
    '<textarea readonly style=\"width:100%;height:90px;' +
    'box-sizing:border-box;margin-top:6px;font:11px monospace\"></textarea>';
  document.body.appendChild(box);
  var ta = box.querySelector('textarea');
  function find(id) {
    return picked.findIndex(function(q) { return q.id === id; });
  }
  function get(pt) {
    var k = pt.curveNumber - first;
    if (k < 0 || k >= n0) return null;
    return {id: ids[k][pt.pointNumber], x: pt.x, y: pt.y};
  }
  var ov = el.data.length;
  Plotly.addTraces(el, {x: [], y: [], type: el.data[first].type,
    mode: 'markers+text',
    text: [], textposition: 'top center', hoverinfo: 'skip',
    showlegend: false, textfont: {size: 16, color: 'black'},
    marker: {size: 14, color: 'rgba(255,255,0,0.6)',
      line: {color: 'black', width: 2}}});
  function draw() {
    Plotly.restyle(el, {
      x: [picked.map(function(q) { return q.x; })],
      y: [picked.map(function(q) { return q.y; })],
      text: [picked.map(function(q, i) {
        if (picked.length <= 2) return ['A', 'B'][i];
        return picked.length <= 20 ? String(i + 1) : '';
      })]
    }, [ov]);
    box.querySelector('.n').textContent = picked.length + ' spot(s)';
    var s = '';
    if (picked.length) {
      s = 'spots <- c(' + picked.map(function(q) {
        return JSON.stringify(q.id); }).join(', ') + ')\\n';
      if (picked.length <= 2) picked.forEach(function(q, i) {
        s += ['A', 'B'][i] + ' <- data.frame(x = ' + q.x +
          ', y = ' + q.y + ')\\n'; });
    }
    ta.value = s;
  }
  el.on('plotly_click', function(e) {
    var q = get(e.points[0]);
    if (!q) return;
    var k = find(q.id);
    if (k >= 0) picked.splice(k, 1); else picked.push(q);
    draw();
  });
  el.on('plotly_selected', function(e) {
    if (!e || !e.points || !e.points.length) return;
    e.points.forEach(function(pt) {
      var q = get(pt);
      if (q && find(q.id) < 0) picked.push(q);
    });
    draw();
  });
  function copy() { ta.select(); document.execCommand('copy'); }
  box.querySelector('.u').onclick = function() { picked.pop(); draw(); };
  box.querySelector('.z').onclick = function() { picked = []; draw(); };
  box.querySelector('.c').onclick = function() {
    if (navigator.clipboard) navigator.clipboard.writeText(ta.value).catch(copy);
    else copy();
  };
  draw();
}"

# Slide image as a data URI, with its size in spatialCoords units.
# Returns NULL when the object holds no image for this sample.
.viewer_image <- function(spe, sample_id, image_id, image_file, scale_factor,
                        max_px) {
    img_data <- SpatialExperiment::imgData(spe)
    if (nrow(img_data) > 0) {
        keep <- img_data$sample_id == sample_id
        if (!is.null(image_id)) keep <- keep & img_data$image_id == image_id
        img_data <- img_data[keep, , drop = FALSE]
    }
    found <- nrow(img_data) > 0
    if (!found && !is.null(image_id)) {
        stop("No image '", image_id, "' for sample '", sample_id, "'.")
    }
    if (is.null(scale_factor) && found) {
        scale_factor <- unname(img_data$scaleFactor[1])
    }
    if (is.null(image_file)) {
        if (!found) return(NULL)
        img <- img_data$data[[1]]
    } else {
        if (!file.exists(image_file)) stop("File not found: ", image_file)
        img <- SpatialExperiment::SpatialImage(image_file, is.url = FALSE)
    }
    if (is.null(scale_factor) || is.na(scale_factor)) {
        stop("'scale_factor' is required: the object does not provide it.")
    }
    img <- magick::image_read(SpatialExperiment::imgRaster(img))
    info <- magick::image_info(img)
    # "WxH>" only shrinks images larger than max_px
    img <- magick::image_resize(img, paste0(max_px, "x", max_px, ">"))
    img <- magick::image_background(img, "white", flatten = TRUE)
    raw <- magick::image_write(img, format = "jpeg", quality = 85)
    list(uri = base64enc::dataURI(raw, mime = "image/jpeg"),
        width = info$width / scale_factor,
        height = info$height / scale_factor)
}

# Split spots into traces: list(rows, name, marker) according to the colouring.
.viewer_traces <- function(value, n, size, alpha, has_image) {
    if (is.null(value)) {
        marker <- list(size = size,
            color = if (has_image) "rgba(0,0,0,0)" else "rgba(150,150,150,0.5)",
            line = list(color = "rgba(60,60,60,0.5)", width = 0.7))
        return(list(list(rows = seq_len(n), name = "spots", marker = marker)))
    }
    grey <- list(size = size, color = "#BBBBBB", opacity = alpha)
    if (is.numeric(value)) {
        ok <- is.finite(value)
        marker <- list(size = size, color = value[ok], opacity = alpha,
            colorscale = "Viridis", showscale = TRUE)
        res <- list(list(rows = which(ok), name = "value", marker = marker),
            list(rows = which(!ok), name = "NA", marker = grey))
    } else {
        value <- droplevels(as.factor(value))
        lev <- levels(value)
        if (length(lev) > 100) {
            stop("'color_by' has more than 100 distinct values.")
        }
        cols <- grDevices::hcl.colors(max(length(lev), 2), "Dark 3")
        res <- lapply(seq_along(lev), function(i) {
            list(rows = which(value == lev[i]), name = lev[i],
                marker = list(size = size, color = cols[i], opacity = alpha))
        })
        res <- c(res, list(list(rows = which(is.na(value)), name = "NA",
            marker = grey)))
    }
    Filter(function(tr) length(tr$rows) > 0, res)
}

#' Interactive viewer to pick spots on the tissue
#'
#' Draws the spots of a SpatialExperiment object on an interactive plotly map,
#' on top of the slide image when there is one, to pick spots and retrieve
#' their identifiers and coordinates.
#'
#' @param spe A SpatialExperiment object.
#' @param color_by Optional name of a `colData` column used to colour the
#'   spots. Numeric columns use a continuous scale, other columns one colour
#'   per value. If `NULL` (default), spots are not coloured: only their outline
#'   is drawn over the slide image.
#' @param image Logical. Draw the slide image? Default is `TRUE`. If the image
#'   stored in the object is missing or cannot be read, the map is drawn
#'   without it.
#' @param sample_id Sample to display. Required if `spe` holds several samples.
#' @param image_id Optional `image_id` in `imgData(spe)`. By default, the first
#'   image of the sample is used.
#' @param image_file Optional path to a slide image, used instead of the one
#'   stored in the object (e.g., when the stored path is no longer valid).
#' @param scale_factor Optional numeric. Image pixels per `spatialCoords` unit.
#'   By default, taken from the selected row of `imgData(spe)`. Set it if
#'   the image has another resolution (1 for a full-resolution image).
#' @param point_size Numeric. Marker size in screen pixels. Default is 6.
#' @param alpha Numeric between 0 and 1. Opacity of coloured spots. Default is
#'   0.8.
#' @param flip_y Logical. Draw the y axis downwards, as in image coordinates?
#'   By default, `TRUE` if an image is drawn and `FALSE` otherwise. This only
#'   affects the display, never the coordinates returned. `FALSE` is not
#'   allowed with an image.
#' @param max_px Integer. Images larger than this are shrunk to keep the page
#'   light. Default is 2000.
#' @param file Optional path of a self-contained HTML page to write.
#'
#' @details
#' Hovering a spot shows its identifier, coordinates and `color_by` value.
#' Clicking a spot adds it to the selection, and clicking it again removes it.
#' The box and lasso tools add all enclosed spots. A panel gives
#' the R code to copy: a `spots` vector of identifiers and, for one or two
#' spots, data frames `A` and `B` that can be passed to
#' [build_one_trajectory()] or [build_similar_trajectories()].
#'
#' Spot identifiers are taken from the `spot_id` column of `colData` if
#' present, otherwise from `colnames(spe)`.
#'
#' This function needs the suggested packages \pkg{SpatialExperiment},
#' \pkg{plotly} and \pkg{htmlwidgets}, and \pkg{magick} and \pkg{base64enc}
#' to draw the slide image.
#'
#' @return A plotly htmlwidget.
#'
#' @examples
#' data("visium_simulated_spe", package = "Battlefield")
#' if (requireNamespace("plotly", quietly = TRUE) &&
#'     requireNamespace("htmlwidgets", quietly = TRUE)) {
#'     # No colour
#'     view_spots(visium_simulated_spe)
#'     # Coloured by a colData column
#'     view_spots(visium_simulated_spe, color_by = "cluster")
#' }
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom methods is
#' @export
view_spots <- function(spe,
                    color_by = NULL,
                    image = TRUE,
                    sample_id = NULL,
                    image_id = NULL,
                    image_file = NULL,
                    scale_factor = NULL,
                    point_size = 6,
                    alpha = 0.8,
                    flip_y = NULL,
                    max_px = 2000,
                    file = NULL) {
    stopifnot(methods::is(spe, "SpatialExperiment"))
    for (pkg in c("SpatialExperiment", "plotly", "htmlwidgets")) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop("Package '", pkg, "' is required.")
        }
    }
    cd <- colData(spe)
    if (!is.null(color_by) &&
        !(length(color_by) == 1 && color_by %in% colnames(cd))) {
        stop("'color_by' must be the name of one colData column.")
    }
    samples <- unique(as.character(cd$sample_id))
    if (is.null(sample_id)) {
        if (length(samples) > 1) {
            stop("Several samples found, set 'sample_id' to one of: ",
                paste(samples, collapse = ", "))
        }
        sample_id <- samples
    }
    if (!(length(sample_id) == 1 && sample_id %in% samples)) {
        stop("'sample_id' not found in the object.")
    }
    keep <- which(cd$sample_id == sample_id)
    ids <- if ("spot_id" %in% colnames(cd)) cd$spot_id else colnames(spe)
    if (is.null(ids)) stop("Spots have no identifier: set colnames(spe).")
    ids <- as.character(ids)[keep]
    if (anyDuplicated(ids)) stop("Spot identifiers are not unique.")
    xy <- SpatialExperiment::spatialCoords(spe)[keep, , drop = FALSE]
    value <- if (!is.null(color_by)) cd[[color_by]][keep]

    img <- NULL
    if (image) {
        load_image <- function() {
            for (pkg in c("magick", "base64enc")) {
                if (!requireNamespace(pkg, quietly = TRUE)) {
                    stop("Package '", pkg, "' is required to draw the image.")
                }
            }
            .viewer_image(spe, sample_id, image_id, image_file, scale_factor,
                max_px)
        }
        has_img <- nrow(SpatialExperiment::imgData(spe)) > 0
        if (!is.null(image_file) || !is.null(image_id)) {
            img <- load_image()
        } else if (has_img) {
            img <- tryCatch(load_image(), error = function(e) {
                warning("Slide image not drawn (see 'image_file'): ",
                    conditionMessage(e), call. = FALSE)
                NULL
            })
        }
    }
    if (is.null(flip_y)) flip_y <- !is.null(img)
    if (!flip_y && !is.null(img)) {
        stop("'flip_y = FALSE' would misalign the image: set image = FALSE.")
    }

    label <- sprintf("<b>%s</b><br>x = %s<br>y = %s", ids,
        format(xy[, 1], trim = TRUE), format(xy[, 2], trim = TRUE))
    if (!is.null(value)) {
        label <- paste0(label, "<br>", color_by, " = ",
            if (is.numeric(value)) signif(value, 4) else as.character(value))
    }

    traces <- .viewer_traces(value, length(ids), point_size, alpha,
        !is.null(img))
    p <- plotly::plot_ly(
        type = if (length(ids) > 5000) "scattergl" else "scatter",
        mode = "markers")
    for (tr in traces) {
        marker <- tr$marker
        if (isTRUE(marker$showscale)) {
            marker$colorbar <- list(title = list(text = color_by))
        }
        p <- plotly::add_trace(p, x = unname(xy[tr$rows, 1]),
            y = unname(xy[tr$rows, 2]), name = tr$name,
            text = label[tr$rows], hoverinfo = "text", marker = marker)
    }
    images <- list()
    xaxis <- list(title = "", zeroline = FALSE, showgrid = FALSE)
    yaxis <- c(xaxis, list(scaleanchor = "x",
        autorange = if (flip_y) "reversed" else TRUE))
    if (!is.null(img)) {
        images <- list(list(source = img$uri, xref = "x", yref = "y",
            x = 0, y = 0, sizex = img$width, sizey = img$height,
            xanchor = "left", yanchor = "top",
            sizing = "stretch", layer = "below"))
        # show the whole image, not only the area covered by spots
        xaxis$range <- c(0, img$width)
        yaxis$autorange <- FALSE
        yaxis$range <- c(img$height, 0)
    }
    p <- plotly::layout(p, dragmode = "pan", hovermode = "closest",
        margin = list(b = 190, t = 30), images = images,
        showlegend = !is.null(value) && !is.numeric(value),
        legend = list(title = list(text = color_by)),
        xaxis = xaxis, yaxis = yaxis)
    p <- plotly::config(p, scrollZoom = TRUE, displaylogo = FALSE,
        modeBarButtonsToAdd = list("select2d", "lasso2d"))
    p <- htmlwidgets::onRender(p, .viewer_js,
        data = lapply(traces, function(tr) I(ids[tr$rows])))
    if (!is.null(file)) htmlwidgets::saveWidget(p, file, selfcontained = TRUE)
    p
}
