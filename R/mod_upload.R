# Upload module: CSV -> column mapping -> observation window -> ppp object.
# Shared by the Density and Risk tabs (both need the same point pattern).

upload_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    # The lead decision -- what data, and whether its x/y are lon/lat (which
    # changes how every downstream field, including boundary buffer's units,
    # is interpreted) -- gets real visual weight. Study region below defaults
    # to a sensible bounding box + 5% padding that most runs won't need to
    # touch, so it starts collapsed rather than competing for attention here.
    shiny::div(
      class = "ic-hero-choice",
      shiny::tags$h5("1. Choose your data"),
      shiny::selectInput(
        ns("example"), "Use an example dataset",
        choices = c("(none — upload my own)" = "", "Burkitt's lymphoma, Uganda (case-control)" = "burk",
                    "Foot-and-mouth disease outbreak" = "fmd", "Primary biliary cirrhosis" = "pbc")
      ),
      shiny::fileInput(ns("csv"), "Or upload a CSV of points", accept = ".csv"),
      shiny::uiOutput(ns("column_map")),
      shiny::checkboxInput(ns("is_geo"), "My x/y columns are longitude/latitude (decimal degrees)", value = FALSE),
      shiny::conditionalPanel(
        condition = "input.is_geo == true", ns = ns,
        shiny::numericInput(ns("epsg"), "Projected CRS to use (EPSG code, blank = auto UTM zone)", value = NA)
      )
    ),
    shiny::div(
      class = "ic-group",
      bslib::accordion(
        open = FALSE,
        bslib::accordion_panel(
          "Study region settings (default: bounding box, 5% padding)",
          shiny::radioButtons(
            ns("window_method"), "Study region (observation window)",
            choices = c("Bounding box around points" = "bbox", "Convex hull of points" = "hull",
                        "Upload a boundary polygon" = "polygon"),
            selected = "bbox"
          ),
          shiny::numericInput(ns("boundary_buffer"), "Boundary buffer / relax boundary", value = 0, min = 0, step = 0.1),
          shiny::helpText(
            "Expands the study window outward by this amount, so points just outside the boundary aren't ",
            "dropped purely because of coordinate recording imprecision (e.g. GPS accuracy of a few hundred ",
            "metres, or a coastline traced at lower resolution than your point data). Applies on top of ",
            "whichever window method is selected below. Units: kilometres if 'longitude/latitude' is checked ",
            "above, otherwise the same units as your x/y columns."
          ),
          shiny::conditionalPanel(
            condition = "input.window_method == 'bbox'", ns = ns,
            shiny::sliderInput(ns("pad_pct"), "Padding around points (%)", min = 0, max = 50, value = 5)
          ),
          shiny::conditionalPanel(
            condition = "input.window_method == 'polygon'", ns = ns,
            shiny::fileInput(
              ns("poly_csv"),
              "Boundary vertex CSV (x, y columns in order around the boundary; an optional 3rd 'part'/'ring'/'island' column names which piece each vertex belongs to, for a boundary with multiple disjoint pieces e.g. a mainland plus its islands)",
              accept = ".csv"
            ),
            shiny::helpText(
              "Vertex winding direction (clockwise vs anticlockwise) is detected and corrected automatically ",
              "per piece, so exports from any GIS tool work regardless of that tool's convention (e.g. Esri/",
              "ArcGIS shapefiles traditionally export clockwise, opposite to the GeoJSON/OGC convention)."
            )
          )
        )
      )
    ),
    shiny::div(
      class = "ic-group",
      shiny::actionButton(ns("build"), "Build point pattern", class = "btn-primary"),
      shiny::hr(),
      shiny::verbatimTextOutput(ns("summary"))
    )
  )
}

upload_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    raw_df <- shiny::reactive({
      if (!is.null(input$csv)) {
        read_uploaded_csv(input$csv$datapath)
      } else if (nzchar(input$example)) {
        utils::read.csv(file.path("data", "examples", paste0(input$example, ".csv")))
      } else {
        NULL
      }
    })

    output$column_map <- shiny::renderUI({
      df <- raw_df()
      shiny::validate(shiny::need(df, "Choose an example dataset or upload a CSV to continue."))
      cols <- names(df)
      shiny::tagList(
        shiny::selectInput(ns("xcol"), "X / longitude column", cols, selected = guess_col(cols, c("x", "lon", "longitude", "easting"))),
        shiny::selectInput(ns("ycol"), "Y / latitude column", cols, selected = guess_col(cols, c("y", "lat", "latitude", "northing"))),
        shiny::selectInput(ns("markcol"), "Group column (optional — needed for case/control risk)",
                            c("(none)" = "", cols), selected = guess_col(cols, c("group", "mark", "type", "case"))),
        shiny::selectInput(ns("weightcol"), "Weight column (optional)", c("(none)" = "", cols)),
        # Defaults to "(none)", unlike the other guessed columns: a wrong
        # x/y/group guess is obvious and easy to fix, but a wrong time-column
        # guess silently makes the Spatiotemporal tabs available with the
        # wrong column selected -- safer to make the user opt in explicitly.
        shiny::selectInput(ns("timecol"), "Time column (optional — needed for the Spatiotemporal tabs)",
                            c("(none)" = "", cols), selected = ""),
        shiny::conditionalPanel(
          condition = "input.timecol != ''", ns = ns,
          shiny::checkboxInput(ns("time_is_date"), "Time column holds dates (e.g. 2020-03-14) — convert to days automatically", value = FALSE)
        )
      )
    })

    result <- shiny::eventReactive(input$build, {
      df <- raw_df()
      shiny::req(df, input$xcol, input$ycol)
      shiny::validate(shiny::need(input$xcol != input$ycol, "X and Y columns must be different."))

      markcol <- if (nzchar(input$markcol %||% "")) input$markcol else NULL
      weightcol <- if (nzchar(input$weightcol %||% "")) input$weightcol else NULL
      timecol <- if (nzchar(input$timecol %||% "")) input$timecol else NULL

      if (!is.null(timecol) && isTRUE(input$time_is_date)) {
        df[[timecol]] <- as.numeric(as.Date(as.character(df[[timecol]])))
      }

      epsg <- NA_integer_
      poly_df <- NULL
      if (identical(input$window_method, "polygon")) {
        shiny::validate(shiny::need(input$poly_csv, "Upload a boundary vertex CSV, or choose a different window method."))
        poly_raw <- read_uploaded_csv(input$poly_csv$datapath)
        poly_df <- poly_raw[, seq_len(min(3, ncol(poly_raw))), drop = FALSE]
      }

      if (isTRUE(input$is_geo)) {
        rp <- reproject_lonlat(
          df, input$xcol, input$ycol, poly_df = poly_df,
          epsg = if (is.na(input$epsg)) NULL else input$epsg
        )
        df <- rp$df
        poly_df <- rp$poly_df
        epsg <- rp$epsg
      }

      x <- suppressWarnings(as.numeric(df[[input$xcol]]))
      y <- suppressWarnings(as.numeric(df[[input$ycol]]))
      ok <- !is.na(x) & !is.na(y)

      window <- switch(input$window_method,
        bbox = owin_from_bbox(x[ok], y[ok], input$pad_pct),
        hull = owin_from_hull(x[ok], y[ok]),
        polygon = owin_from_polygon(poly_df)
      )

      buffer_val <- if (is.null(input$boundary_buffer) || is.na(input$boundary_buffer)) 0 else input$boundary_buffer
      if (buffer_val > 0) {
        window <- apply_boundary_buffer(window, buffer_val, isTRUE(input$is_geo))
      }

      built <- build_ppp(df, input$xcol, input$ycol, window, markcol = markcol, weightcol = weightcol, timecol = timecol)

      list(
        pp = built$pp, weights = built$weights, tt = built$tt, n_dropped = built$n_dropped,
        is_geo = isTRUE(input$is_geo), epsg = epsg,
        xcol = input$xcol, ycol = input$ycol, markcol = markcol, weightcol = weightcol,
        timecol = timecol, time_is_date = isTRUE(input$time_is_date),
        window_method = input$window_method, pad_pct = input$pad_pct, boundary_buffer = buffer_val
      )
    })

    output$summary <- shiny::renderPrint({
      r <- result()
      cat("Points in window:", spatstat.geom::npoints(r$pp), "\n")
      cat("Points dropped (invalid or outside window):", r$n_dropped, "\n")
      cat("Window area:", format(spatstat.geom::area(r$pp$window), digits = 4), "\n")
      if (isTRUE(r$boundary_buffer > 0)) {
        cat("Boundary buffer applied:", r$boundary_buffer, if (isTRUE(r$is_geo)) "km\n" else "map units\n")
      }
      if (!is.null(r$markcol)) {
        cat("\nGroup counts:\n")
        print(table(spatstat.geom::marks(r$pp)))
      }
      if (!is.null(r$timecol)) {
        n_timed <- sum(!is.na(r$tt))
        cat("\nTime column '", r$timecol, "': ", n_timed, " of ", length(r$tt),
            " points have a time value", sep = "")
        if (n_timed > 0) cat(", range ", paste(signif(range(r$tt, na.rm = TRUE), 6), collapse = " to "), sep = "")
        cat("\n")
      }
      if (isTRUE(r$is_geo)) cat("\nReprojected to EPSG:", r$epsg, "\n")
    })

    result
  })
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

guess_col <- function(cols, candidates) {
  hit <- cols[tolower(cols) %in% candidates]
  if (length(hit) > 0) hit[1] else cols[1]
}
