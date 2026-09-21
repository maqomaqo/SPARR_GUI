# Generic results display: Plot / Map / Table / Download tabs. Shared by the
# Density and Risk modules — both hand it a `bivden` or `rrs` sparr object.
# `result` is a reactive returning the sparr object (or NULL before a run).
# `meta` is a reactive returning list(analysis, is_geo, epsg, repro_params).

results_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fluidRow(
      shiny::column(4, shiny::checkboxInput(ns("show_points"), "Show observed points on plot", value = FALSE)),
      shiny::column(8, shiny::uiOutput(ns("tol_ui")))
    ),
    shiny::tabsetPanel(
      # ic-results-frame caps Plot/Table to a sensible viewing/reading width
      # instead of stretching a fixed-size plot image or a short summary
      # table across whatever width sidebarLayout's mainPanel leaves (which
      # grows with the browser window). Map is left unconstrained -- more
      # width is genuinely useful for a pannable basemap.
      shiny::tabPanel("Plot", shiny::div(
        class = "ic-result-pane ic-results-frame",
        shinycssloaders::withSpinner(shiny::plotOutput(ns("plot"), height = 550))
      )),
      shiny::tabPanel("Map", shiny::div(
        class = "ic-result-pane",
        shinycssloaders::withSpinner(leaflet::leafletOutput(ns("map"), height = 550))
      )),
      shiny::tabPanel("Table", shiny::div(
        class = "ic-result-pane ic-results-frame",
        DT::DTOutput(ns("table"))
      )),
      shiny::tabPanel(
        "Download",
        shiny::br(),
        shiny::downloadButton(ns("dl_png"), "Plot (PNG)"),
        shiny::downloadButton(ns("dl_pdf"), "Plot (PDF)"),
        shiny::downloadButton(ns("dl_grid"), "Surface grid (CSV)"),
        shiny::downloadButton(ns("dl_script"), "Reproducible R script"),
      )
    )
  )
}

results_server <- function(id, result, meta) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$tol_ui <- shiny::renderUI({
      r <- result()
      if (is.null(r) || !inherits(r, "rrs") || is.null(r$P)) return(NULL)
      shiny::tagList(
        shiny::sliderInput(ns("tol_levels"), "Tolerance contour significance level(s)",
                            min = 0.001, max = 0.2, value = c(0.05), step = 0.001),
        shiny::checkboxInput(ns("show_lower"), "Also outline areas of lower significance (dotted)", value = FALSE)
      )
    })

    tol_levels <- shiny::reactive({
      if (!is.null(input$tol_levels)) input$tol_levels else c(0.05, 0.01)
    })
    show_lower <- shiny::reactive(isTRUE(input$show_lower))

    output$plot <- shiny::renderPlot({
      r <- result()
      shiny::validate(shiny::need(r, "Run an analysis to see the plot."))
      render_sparr_plot(r, tol_levels = tol_levels(), show_points = isTRUE(input$show_points), show_lower = show_lower())
      # Fixed width: on a tab pane's first-ever paint the client can briefly
      # report 0 for the container's measured width (before its resize
      # observer has fired), which otherwise fails the PNG device outright.
      # A fixed width sidesteps that negotiation entirely.
    }, width = 900, height = 550)

    output$map <- leaflet::renderLeaflet({
      r <- result()
      shiny::validate(shiny::need(r, "Run an analysis to see the map."))
      m <- meta()
      if (!isTRUE(m$is_geo)) {
        return(
          leaflet::leaflet() |>
            leaflet::addProviderTiles(leaflet::providers$CartoDB.Positron) |>
            leaflet::setView(0, 0, zoom = 1)
        )
      }
      imobj <- extract_surface_im(r)
      ras <- im_to_wgs84_raster(imobj, m$epsg)
      pal <- leaflet::colorNumeric(sparr_colour_ramp()(256), terra::values(ras), na.color = "transparent")
      leaflet::leaflet() |>
        leaflet::addProviderTiles(leaflet::providers$CartoDB.Positron) |>
        leaflet::addRasterImage(ras, colors = pal, opacity = 0.75) |>
        leaflet::addLegend(pal = pal, values = terra::values(ras), title = "Surface value")
    })

    output$table <- DT::renderDT({
      r <- result()
      shiny::validate(shiny::need(r, "Run an analysis to see summary statistics."))
      imobj <- extract_surface_im(r)
      all_vals <- imobj$v[!is.na(imobj$v)]
      # Log relative risk is legitimately +-Inf where the control (or case)
      # density estimate rounds to zero at a pixel; exclude those from the
      # summary stats (min/max/mean would otherwise just read Inf/-Inf/NaN)
      # but report how many pixels that affected.
      vals <- all_vals[is.finite(all_vals)]
      n_infinite <- sum(!is.finite(all_vals))
      stats <- data.frame(
        Statistic = c("Grid resolution", "Min", "Max", "Mean", "Median"),
        Value = c(
          paste0(imobj$dim[2], " x ", imobj$dim[1]),
          signif(min(vals), 4), signif(max(vals), 4),
          signif(mean(vals), 4), signif(stats::median(vals), 4)
        )
      )
      if (n_infinite > 0) {
        stats <- rbind(stats, data.frame(
          Statistic = "Non-finite pixels excluded above (±Inf)", Value = n_infinite
        ))
      }
      if (inherits(r, "rrs") && !is.null(r$P)) {
        for (lev in tol_levels()) {
          stats <- rbind(stats, data.frame(
            Statistic = paste0("Pixels significant at p<", lev),
            Value = sum(r$P$v < lev, na.rm = TRUE)
          ))
        }
      }
      # An explicit width, not "100%" of .ic-results-frame: this is a
      # 2-column, half-dozen-row summary table, not a data grid -- letting
      # it fill the frame would just spread "Statistic" and "Value" apart
      # with a wide gap of nothing between them.
      DT::datatable(stats, options = list(dom = "t", paging = FALSE), rownames = FALSE, width = "480px")
    })
    for (nm in c("plot", "map", "table")) {
      shiny::outputOptions(output, nm, suspendWhenHidden = FALSE)
    }

    output$dl_png <- shiny::downloadHandler(
      filename = function() "sparr_plot.png",
      content = function(file) write_plot_png(file, result(), tol_levels(), show_lower())
    )
    output$dl_pdf <- shiny::downloadHandler(
      filename = function() "sparr_plot.pdf",
      content = function(file) write_plot_pdf(file, result(), tol_levels(), show_lower())
    )
    output$dl_grid <- shiny::downloadHandler(
      filename = function() "sparr_surface_grid.csv",
      content = function(file) write_grid_csv(file, extract_surface_im(result()))
    )
    output$dl_script <- shiny::downloadHandler(
      filename = function() "sparr_analysis.R",
      content = function(file) writeLines(generate_repro_script(meta()$repro_params), file)
    )
    # Download links live in a tab pane that isn't shown by default; leaving
    # suspendWhenHidden on its default can leave the href empty (disabled)
    # even after the pane becomes visible. See mod_risk.R for the same fix.
    for (nm in c("dl_png", "dl_pdf", "dl_grid", "dl_script")) {
      shiny::outputOptions(output, nm, suspendWhenHidden = FALSE)
    }
  })
}
