# Download handlers: plot image export, raster grid CSV export, and a
# generated R script reproducing the analysis just run outside the app.

write_plot_png <- function(path, result, tol_levels, show_lower = FALSE, width = 1000, height = 850) {
  grDevices::png(path, width = width, height = height, res = 120)
  on.exit(grDevices::dev.off())
  render_sparr_plot(result, tol_levels = tol_levels, show_lower = show_lower)
}

write_plot_pdf <- function(path, result, tol_levels, show_lower = FALSE, width = 9, height = 7.5) {
  grDevices::pdf(path, width = width, height = height)
  on.exit(grDevices::dev.off())
  render_sparr_plot(result, tol_levels = tol_levels, show_lower = show_lower)
}

write_grid_csv <- function(path, imobj) {
  readr::write_csv(im_to_grid_df(imobj), path)
}

write_spattemp_plot_png <- function(path, result, t_sel, type, tol_levels, show_lower = FALSE, width = 1000, height = 850) {
  grDevices::png(path, width = width, height = height, res = 120)
  on.exit(grDevices::dev.off())
  render_spattemp_plot(result, t_sel, type = type, tol_levels = tol_levels, show_lower = show_lower)
}

write_spattemp_plot_pdf <- function(path, result, t_sel, type, tol_levels, show_lower = FALSE, width = 9, height = 7.5) {
  grDevices::pdf(path, width = width, height = height)
  on.exit(grDevices::dev.off())
  render_spattemp_plot(result, t_sel, type = type, tol_levels = tol_levels, show_lower = show_lower)
}

#' Build a standalone R script string that reproduces the analysis the app
#' just ran, so the user can graduate to scripting sparr/spatstat directly.
#' `p` is a named list of the parameters used (see mod_density.R / mod_risk.R).
generate_repro_script <- function(p) {
  lines <- c(
    "# Reproduces the analysis run in the sparr GUI app.",
    "library(sparr)",
    "library(spatstat.geom)",
    "",
    "df <- read.csv(\"your_data.csv\")",
    ""
  )

  window_code <- switch(p$window_method,
    bbox = sprintf(
      "w <- owin(xrange = range(df$%s) + c(-1,1)*diff(range(df$%s))*%s/100, yrange = range(df$%s) + c(-1,1)*diff(range(df$%s))*%s/100)",
      p$xcol, p$xcol, p$pad_pct, p$ycol, p$ycol, p$pad_pct
    ),
    hull = sprintf("w <- convexhull.xy(df$%s, df$%s)", p$xcol, p$ycol),
    polygon = "w <- owin(poly = list(x = boundary$x, y = boundary$y))  # boundary <- read.csv(\"boundary.csv\")"
  )
  lines <- c(lines, window_code)

  if (isTRUE(p$boundary_buffer > 0)) {
    r_expr <- if (isTRUE(p$is_geo)) {
      sprintf("%s * 1000  # km -> metres; assumes df$%s/df$%s are already in a projected (metric) CRS here",
              p$boundary_buffer, p$xcol, p$ycol)
    } else {
      as.character(p$boundary_buffer)
    }
    lines <- c(lines, sprintf("w <- dilation(w, %s)  # relax boundary to absorb near-boundary points", r_expr))
  }
  lines <- c(lines, "")

  is_case_control <- identical(p$analysis, "risk") || identical(p$analysis, "spattemp_risk")
  if (!is.null(p$markcol) && is_case_control && !is.null(p$case_level)) {
    # sparr::risk()/spattemp.risk() treat the FIRST factor level as the case
    # group, so the level order here must match whichever group was chosen
    # as "case" in the app -- not just the data's natural/alphabetical order.
    lines <- c(lines,
      sprintf(
        "pp <- ppp(df$%s, df$%s, window = w, marks = factor(df$%s, levels = c(%s, %s)))  # %s = case, %s = control",
        p$xcol, p$ycol, p$markcol, deparse(p$case_level), deparse(p$control_level),
        deparse(p$case_level), deparse(p$control_level)
      )
    )
  } else if (!is.null(p$markcol)) {
    lines <- c(lines,
      sprintf("pp <- ppp(df$%s, df$%s, window = w, marks = factor(df$%s))", p$xcol, p$ycol, p$markcol)
    )
  } else {
    lines <- c(lines, sprintf("pp <- ppp(df$%s, df$%s, window = w)", p$xcol, p$ycol))
  }

  is_spattemp <- identical(p$analysis, "spattemp_density") || identical(p$analysis, "spattemp_risk")
  if (is_spattemp && !is.null(p$timecol)) {
    if (isTRUE(p$time_is_date)) {
      lines <- c(lines, sprintf("tt <- as.numeric(as.Date(df$%s))", p$timecol))
    } else {
      lines <- c(lines, sprintf("tt <- df$%s", p$timecol))
    }
  }
  lines <- c(lines, "")

  has_subset <- (identical(p$analysis, "density") || identical(p$analysis, "spattemp_density")) &&
    !is.null(p$markcol) && !is.null(p$subset_level) && !identical(p$subset_level, "__all__")
  weights_expr <- if (!is.null(p$weightcol)) sprintf("df$%s", p$weightcol) else NULL
  if (has_subset) {
    lines <- c(lines, sprintf("keep <- marks(pp) == %s  # restrict to the chosen population", deparse(p$subset_level)))
    if (identical(p$analysis, "spattemp_density")) {
      lines <- c(lines, "pp <- pp[keep]; tt <- tt[keep]")
    } else {
      lines <- c(lines, "pp <- pp[keep]")
    }
    lines <- c(lines, "")
    if (!is.null(weights_expr)) weights_expr <- paste0(weights_expr, "[keep]")
  }
  if (is_spattemp && !is.null(p$timecol)) {
    lines <- c(lines,
      "has_time <- !is.na(tt)  # spattemp functions need a time value for every point",
      "pp <- pp[has_time]; tt <- tt[has_time]", ""
    )
  }

  if (identical(p$analysis, "density")) {
    args <- c(
      sprintf("h0 = %s", p$h0),
      if (isTRUE(p$adapt)) "adapt = TRUE",
      if (isTRUE(p$adapt) && !is.null(p$hp)) sprintf("hp = %s", p$hp),
      sprintf("resolution = %s", p$resolution),
      sprintf("edge = \"%s\"", p$edge),
      if (!is.null(weights_expr)) sprintf("weights = %s", weights_expr)
    )
    lines <- c(lines,
      sprintf("result <- bivariate.density(pp, %s)", paste(args, collapse = ", ")),
      "plot(result)"
    )
  } else if (identical(p$analysis, "risk")) {
    args <- c(
      sprintf("h0 = %s", p$h0),
      if (isTRUE(p$adapt)) "adapt = TRUE",
      sprintf("resolution = %s", p$resolution),
      sprintf("tolerate = %s", if (isTRUE(p$tolerate)) "TRUE" else "FALSE")
    )
    lines <- c(lines,
      sprintf("result <- risk(pp, %s)", paste(args, collapse = ", ")),
      "plot(result)"
    )
  } else if (identical(p$analysis, "spattemp_density")) {
    args <- c(
      sprintf("h = %s", p$h), sprintf("lambda = %s", p$lambda), "tt = tt",
      sprintf("sedge = \"%s\"", p$edge), sprintf("tedge = \"%s\"", p$edge),
      sprintf("sres = %s", p$sres), sprintf("tres = %s", p$tres)
    )
    lines <- c(lines,
      sprintf("result <- spattemp.density(pp, %s)", paste(args, collapse = ", ")),
      "plot(result, tselect = mean(result$tlim))  # pick any time within result$tlim"
    )
  } else if (identical(p$analysis, "spattemp_risk")) {
    lines <- c(lines,
      sprintf("case_pp <- pp[marks(pp) == %s]; case_tt <- tt[marks(pp) == %s]", deparse(p$case_level), deparse(p$case_level)),
      sprintf("control_pp <- pp[marks(pp) == %s]", deparse(p$control_level)),
      "",
      sprintf(
        "f <- spattemp.density(case_pp, tt = case_tt, h = %s, lambda = %s, sres = %s, tres = %s)",
        p$h, p$lambda, p$resolution, p$tres
      ),
      sprintf("g <- bivariate.density(control_pp, h0 = %s, resolution = %s)  # static control", p$h0, p$resolution),
      sprintf("result <- spattemp.risk(f, g, tolerate = %s)", if (isTRUE(p$tolerate)) "TRUE" else "FALSE"),
      "plot(result, tselect = mean(result$tlim))  # pick any time within result$tlim"
    )
  }

  paste(lines, collapse = "\n")
}
