# Plot/raster helpers shared by the density and risk result tabs.
# sparr's own S3 plot methods (plot.bivden, plot.rrs) already produce the
# right image+contour output via spatstat.geom::plot.im -- these helpers just
# wrap consistent styling (palette, tolerance contour args) around them, and
# convert the underlying pixel image to formats needed for the Leaflet map
# tab and the CSV grid export.

#' A colourblind-safe ramp function suitable for spatstat's `col` argument.
sparr_colour_ramp <- function(option = "viridis") {
  function(n) viridisLite::viridis(n, option = option)
}

#' Overlay contours for the opposite tail of a tolerance p-value surface --
#' i.e. areas of significantly *low* risk, where sparr's own tol.args only
#' ever contours the significantly *high* side (tol.type = "upper"). Mirrors
#' the matrix orientation sparr's plot.rrs/plot.rrst use internally (transpose
#' `P$v` before passing to `contour`), but takes the complement (1 - p) and
#' draws dotted (lty = 3) so it reads as visually distinct from the solid
#' upper-tail contours already on the plot.
add_lower_tol_contours <- function(P, tol_levels) {
  if (is.null(P)) return(invisible(NULL))
  ps <- 1 - t(as.matrix(P))
  suppressWarnings(
    graphics::contour(x = P$xcol, y = P$yrow, z = ps, levels = tol_levels,
                       lty = 3, drawlabels = TRUE, add = TRUE)
  )
  invisible(NULL)
}

#' Render a `bivden` (density) or `rrs` (risk) sparr object using its own
#' plot method, with consistent styling. `tol_levels` and `show_lower` only
#' apply to `rrs` objects that were computed with tolerate = TRUE.
render_sparr_plot <- function(x, tol_levels = c(0.05, 0.01), show_points = FALSE, show_lower = FALSE) {
  ramp <- sparr_colour_ramp()
  if (inherits(x, "rrs")) {
    has_p <- !is.null(x$P)
    plot(
      x,
      col = ramp,
      tol.show = has_p,
      tol.args = list(levels = tol_levels, lty = seq_along(tol_levels), drawlabels = TRUE)
    )
    if (has_p && show_lower) add_lower_tol_contours(x$P, tol_levels)
  } else {
    plot(x, col = ramp, add.pts = show_points)
  }
  invisible(NULL)
}

#' Extract the primary surface (`im` pixel image) from a bivden or rrs object.
extract_surface_im <- function(x) {
  if (inherits(x, "rrs")) x$rr else x$z
}

#' Long-format data.frame (x, y, z) of a spatstat pixel image, for CSV export.
im_to_grid_df <- function(imobj) {
  # as.data.frame.im is an S3 method registered by spatstat.geom (not itself
  # exported as a standalone function), so dispatch through the base generic.
  xy <- as.data.frame(imobj)
  names(xy) <- c("x", "y", "z")
  xy
}

#' The discretised time grid of a stden or rrst object. `stden` stores this
#' at the top level ($tgrid); `rrst` (sparr::spattemp.risk()'s output) does
#' not repeat it there -- it only lives on the case density nested at $f.
spattemp_tgrid <- function(x) {
  if (!is.null(x$tgrid)) x$tgrid else x$f$tgrid
}

#' Render a `stden` (spatiotemporal density) or `rrst` (spatiotemporal risk)
#' sparr object at a chosen time `t_sel`, using its own plot method. `type` is
#' "joint" (unconditional) or "conditional". Neither class supports add.pts.
#' `p_slice`, if supplied, is the already-interpolated P surface for `t_sel`
#' (e.g. a caller's own cached `spattemp_slice_im()$P`) -- lets callers that
#' already sliced P for this exact time/type skip doing it again here.
render_spattemp_plot <- function(x, t_sel, type = "joint", tol_levels = c(0.05, 0.01), show_lower = FALSE, p_slice = NULL) {
  ramp <- sparr_colour_ramp()
  if (inherits(x, "rrst")) {
    has_p <- !is.null(x$P)
    draw_lower <- has_p && show_lower
    if (draw_lower) {
      # plot.rrst's default override.par=TRUE captures the *pre-plot* par()
      # (including par("usr"), the coordinate system) and restores it via
      # on.exit as soon as the call returns -- so a contour(add=TRUE) call
      # made after plot() would draw into the wrong coordinate space and be
      # invisible. Reproduce its override.par=TRUE margins ourselves and pass
      # override.par=FALSE so the coordinate system survives until our own
      # on.exit restores it once we've finished drawing the extra contour.
      old_par <- par(mfrow = c(1, 1), mar = rep(2, 4))
      on.exit(par(old_par))
    }
    plot(
      x, tselect = t_sel, type = type, col = ramp,
      tol.show = has_p,
      tol.args = list(levels = tol_levels, lty = seq_along(tol_levels), drawlabels = TRUE),
      override.par = !draw_lower
    )
    if (draw_lower) {
      # Same slice() path the Map/Table tabs use (spattemp_slice_im), so the
      # overlaid contour lines up with whatever plot.rrst just drew above --
      # it interpolates for a t_sel that isn't exactly on x's own tgrid.
      P <- if (!is.null(p_slice)) p_slice else spattemp_slice_im(x, t_sel, type)$P
      add_lower_tol_contours(P, tol_levels)
    }
  } else {
    plot(x, tselect = t_sel, type = type, col = ramp)
  }
  invisible(NULL)
}

#' Extract the pixel image (and, for risk, its p-value surface) at one time
#' from a stden/rrst object, for the Table/Download tabs. `spattemp.slice()`
#' interpolates for a time that isn't exactly on the object's own tgrid.
spattemp_slice_im <- function(x, t_sel, type = "joint") {
  sl <- sparr::spattemp.slice(x, tt = t_sel)
  if (inherits(x, "rrst")) {
    z <- if (identical(type, "conditional")) sl$rr.cond[[1]] else sl$rr[[1]]
    P <- if (!is.null(x$P)) {
      if (identical(type, "conditional")) sl$P.cond[[1]] else sl$P[[1]]
    } else NULL
    list(z = z, P = P)
  } else {
    list(z = if (identical(type, "conditional")) sl$z.cond[[1]] else sl$z[[1]], P = NULL)
  }
}

#' Convert a spatstat pixel image (planar, in a projected CRS `epsg`) to a
#' `terra` SpatRaster reprojected to WGS84 (EPSG:4326), for Leaflet display.
im_to_wgs84_raster <- function(imobj, epsg) {
  # spatstat's im$v has rows sorted by increasing y (bottom row = min y);
  # terra::rast() on a plain matrix expects the opposite (top row = max y).
  v_flipped <- imobj$v[nrow(imobj$v):1, , drop = FALSE]
  r <- terra::rast(
    v_flipped,
    extent = terra::ext(imobj$xrange[1], imobj$xrange[2], imobj$yrange[1], imobj$yrange[2]),
    crs = paste0("EPSG:", epsg)
  )
  terra::project(r, "EPSG:4326")
}
