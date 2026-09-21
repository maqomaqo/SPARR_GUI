# Helpers for turning an uploaded data.frame into a spatstat ppp/owin,
# in either planar (as-is) or geographic (lat/lon -> projected) mode.

#' Read a CSV robustly (handles BOM, guesses types) and return a data.frame.
read_uploaded_csv <- function(path) {
  df <- as.data.frame(readr::read_csv(path, show_col_types = FALSE, progress = FALSE))
  if (nrow(df) == 0) stop("The uploaded CSV has no data rows.")
  df
}

#' Pick a planar (metric) CRS for a set of lon/lat points: the UTM zone
#' containing their centroid. Good enough for a study-region-sized dataset.
utm_epsg_for_lonlat <- function(lon, lat) {
  zone <- floor((mean(lon, na.rm = TRUE) + 180) / 6) %% 60 + 1
  north <- mean(lat, na.rm = TRUE) >= 0
  base <- if (north) 32600 else 32700
  base + zone
}

#' Expand a window outward by a fixed distance (morphological dilation), to
#' absorb points that fall just outside the true boundary purely because of
#' coordinate recording imprecision. `dist_val` is interpreted as kilometres
#' when `is_geo` is TRUE (coordinates have already been reprojected to a
#' metric CRS by this point, i.e. metres) and otherwise as the window's own
#' map units directly.
apply_boundary_buffer <- function(window, dist_val, is_geo) {
  r <- if (isTRUE(is_geo)) dist_val * 1000 else dist_val
  spatstat.geom::dilation(window, r)
}

#' Build an owin from a bounding box around x,y with `pad_pct` percent padding.
owin_from_bbox <- function(x, y, pad_pct = 5) {
  xr <- range(x); yr <- range(y)
  xpad <- diff(xr) * pad_pct / 100
  ypad <- diff(yr) * pad_pct / 100
  if (xpad == 0) xpad <- 1
  if (ypad == 0) ypad <- 1
  spatstat.geom::owin(xrange = xr + c(-xpad, xpad), yrange = yr + c(-ypad, ypad))
}

#' Build an owin as the convex hull of x,y.
owin_from_hull <- function(x, y) {
  spatstat.geom::convexhull.xy(x, y)
}

#' Signed area of a polygon ring via the shoelace formula: positive for
#' anticlockwise winding, negative for clockwise.
ring_signed_area <- function(x, y) {
  0.5 * sum(x * c(y[-1], y[1]) - c(x[-1], x[1]) * y)
}

#' Ensure a ring winds anticlockwise (spatstat's required orientation for an
#' exterior boundary), reversing it with spatstat.utils::reverse.xypolygon()
#' if it doesn't. GIS tools disagree on winding convention for exported
#' polygons -- Esri/ArcGIS shapefiles traditionally use clockwise exterior
#' rings, the opposite of the GeoJSON/OGC convention -- so this can't be
#' assumed either way and has to be checked per ring.
normalize_ring_orientation <- function(ring) {
  if (ring_signed_area(ring$x, ring$y) < 0) spatstat.utils::reverse.xypolygon(ring) else ring
}

#' Build a single-ring owin, normalizing orientation first (see above).
owin_from_ring <- function(x, y) {
  r <- normalize_ring_orientation(list(x = x, y = y))
  spatstat.geom::owin(poly = list(x = r$x, y = r$y))
}

#' Build an owin from a polygon boundary given as a vertex data.frame with
#' x,y columns (in order around the boundary), and optionally a 3rd column
#' naming which ring/part each vertex belongs to (for a boundary with
#' multiple disjoint pieces, e.g. a mainland coastline plus its islands).
#'
#' If there's no explicit ring column, multi-part boundaries are common
#' anyway when a CSV was exported by flattening a multi-polygon shapefile
#' without keeping its ring/part index (e.g. Reduce()-ing over st_coordinates()
#' output columns) -- the vertex list is just several separate closed rings
#' concatenated with no separator, sometimes with a ring fully *nested*
#' inside a run of another ring's vertices (traced, then resumed) rather
#' than simply appended after it. Naively treating the whole file as one
#' ring connects unrelated, distant parts of the boundary with a straight
#' line -- both at each such seam and at the final implicit closing edge.
#' This is detected and repaired here: anomalously large jumps between
#' consecutive vertices mark ring boundaries, and a chunk between two such
#' jumps is either its own closed ring (first vertex ~= last vertex, e.g. an
#' inserted island) or an open fragment that must be re-joined with the
#' fragments before/after it to reconstitute the ring it was cut out of.
owin_from_polygon <- function(poly_df) {
  if (ncol(poly_df) >= 3 && tolower(names(poly_df)[3]) %in%
      c("part", "ring", "group", "id", "piece", "polygon", "island")) {
    rings <- split(poly_df[1:2], poly_df[[3]])
    ring_list <- lapply(rings, function(r) list(x = r[[1]], y = r[[2]]))
    ring_list <- Filter(function(r) length(r$x) >= 3, ring_list)
    if (length(ring_list) == 0) stop("No valid polygon rings found in the boundary file.")
    if (length(ring_list) == 1) return(owin_from_ring(ring_list[[1]]$x, ring_list[[1]]$y))
    ring_list <- lapply(ring_list, normalize_ring_orientation)
    return(spatstat.geom::owin(poly = ring_list))
  }

  px <- poly_df[[1]]; py <- poly_df[[2]]
  n <- length(px)
  if (n < 3) stop("Boundary file needs at least 3 vertices.")

  gaps <- sqrt(diff(px)^2 + diff(py)^2)
  med <- stats::median(gaps[gaps > 0], na.rm = TRUE)
  breaks <- if (is.finite(med) && med > 0) which(gaps > med * 25) else integer(0)

  if (length(breaks) == 0) return(owin_from_ring(px, py))

  starts <- c(1, breaks + 1)
  ends <- c(breaks, n)
  eps <- max(med * 5, .Machine$double.eps * 100)
  close_enough <- function(i, j) sqrt((px[i] - px[j])^2 + (py[i] - py[j])^2) < eps

  ring_list <- list()
  buf_idx <- integer(0)
  for (k in seq_along(starts)) {
    s <- starts[k]; e <- ends[k]
    if ((e - s) >= 2 && close_enough(s, e)) {
      ring_list[[length(ring_list) + 1]] <- list(x = px[s:e], y = py[s:e])
    } else {
      buf_idx <- c(buf_idx, s:e)
      if (length(buf_idx) >= 3 && close_enough(buf_idx[1], buf_idx[length(buf_idx)])) {
        ring_list[[length(ring_list) + 1]] <- list(x = px[buf_idx], y = py[buf_idx])
        buf_idx <- integer(0)
      }
    }
  }
  if (length(buf_idx) >= 3) ring_list[[length(ring_list) + 1]] <- list(x = px[buf_idx], y = py[buf_idx])

  ring_list <- Filter(function(r) length(r$x) >= 3, ring_list)
  if (length(ring_list) == 0) stop("No valid polygon rings found in the boundary file.")
  if (length(ring_list) == 1) return(owin_from_ring(ring_list[[1]]$x, ring_list[[1]]$y))
  ring_list <- lapply(ring_list, normalize_ring_orientation)
  spatstat.geom::owin(poly = ring_list)
}

#' Build a planar ppp from a data.frame, dropping points outside `window`.
#' `timecol`, if given, is tracked as a plain numeric vector aligned to the
#' final point order (NOT stored as ppp marks -- kept separate from
#' `markcol`'s case/control factor since sparr's spatiotemporal functions
#' take time via an explicit `tt=` argument). NA time values are kept as-is
#' (e.g. a static control group legitimately has no event time), not treated
#' as invalid data to drop.
#' Returns list(pp = <ppp>, n_dropped = <int>).
build_ppp <- function(df, xcol, ycol, window, markcol = NULL, weightcol = NULL, timecol = NULL) {
  x <- suppressWarnings(as.numeric(df[[xcol]]))
  y <- suppressWarnings(as.numeric(df[[ycol]]))
  bad_xy <- is.na(x) | is.na(y)
  if (all(bad_xy)) stop("Coordinate columns '", xcol, "'/'", ycol, "' contain no valid numbers.")

  keep <- !bad_xy
  x <- x[keep]; y <- y[keep]
  marks_vec <- if (!is.null(markcol)) df[[markcol]][keep] else NULL
  weights_vec <- if (!is.null(weightcol)) suppressWarnings(as.numeric(df[[weightcol]][keep])) else NULL
  time_vec <- if (!is.null(timecol)) suppressWarnings(as.numeric(df[[timecol]][keep])) else NULL

  inside <- spatstat.geom::inside.owin(x, y, window)
  n_dropped <- sum(bad_xy) + sum(!inside)

  x <- x[inside]; y <- y[inside]
  if (length(x) < 2) stop("Fewer than 2 points fall inside the chosen study window.")

  if (!is.null(marks_vec)) marks_vec <- marks_vec[inside]
  if (!is.null(weights_vec)) weights_vec <- weights_vec[inside]
  if (!is.null(time_vec)) time_vec <- time_vec[inside]

  pp <- if (!is.null(marks_vec)) {
    spatstat.geom::ppp(x, y, window = window, marks = factor(marks_vec), checkdup = FALSE)
  } else {
    spatstat.geom::ppp(x, y, window = window, checkdup = FALSE)
  }

  list(pp = pp, weights = weights_vec, tt = time_vec, n_dropped = n_dropped)
}

#' Reproject a lon/lat data.frame's coordinate + optional polygon vertices to
#' a metric CRS (UTM zone by default, or an explicit EPSG code), returning
#' plain numeric x/y (already updated in-place on copies of the inputs) plus
#' the EPSG code used. Requires the `sf` package.
reproject_lonlat <- function(df, loncol, latcol, poly_df = NULL, epsg = NULL) {
  lon <- suppressWarnings(as.numeric(df[[loncol]]))
  lat <- suppressWarnings(as.numeric(df[[latcol]]))
  if (is.null(epsg)) epsg <- utm_epsg_for_lonlat(lon, lat)

  pts_sf <- sf::st_as_sf(data.frame(lon = lon, lat = lat), coords = c("lon", "lat"), crs = 4326)
  pts_proj <- sf::st_transform(pts_sf, epsg)
  coords <- sf::st_coordinates(pts_proj)

  df[[loncol]] <- coords[, 1]
  df[[latcol]] <- coords[, 2]

  poly_out <- NULL
  if (!is.null(poly_df)) {
    poly_sf <- sf::st_as_sf(
      data.frame(lon = poly_df[[1]], lat = poly_df[[2]]),
      coords = c("lon", "lat"), crs = 4326
    )
    poly_proj <- sf::st_transform(poly_sf, epsg)
    pc <- sf::st_coordinates(poly_proj)
    poly_out <- data.frame(x = pc[, 1], y = pc[, 2])
    # Carry an optional 3rd ring/part-id column through reprojection, under
    # its original name (owin_from_polygon() looks for it by name) --
    # st_coordinates() drops all non-geometry columns, so reattach it by
    # position (row order is preserved by st_transform/st_coordinates).
    if (ncol(poly_df) >= 3) poly_out[[names(poly_df)[3]]] <- poly_df[[3]]
  }

  list(df = df, poly_df = poly_out, epsg = epsg)
}

#' Validate that a mark column is usable as a case/control factor for risk().
#' Returns the 2 level names (case first) or throws an informative error.
validate_case_control_marks <- function(marks_vec) {
  lv <- levels(factor(marks_vec))
  if (length(lv) != 2) {
    stop(
      "The group column must have exactly 2 distinct values for case/control ",
      "risk estimation (found ", length(lv), ": ", paste(lv, collapse = ", "), ")."
    )
  }
  lv
}
