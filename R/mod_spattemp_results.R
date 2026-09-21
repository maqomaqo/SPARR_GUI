# Spatiotemporal results display: time slider (+ play/pause animation),
# joint/conditional toggle, then Plot / Map / Table / Download tabs. Shared
# by the Spatiotemporal Density and Risk modules -- both hand it a `stden`
# or `rrst` sparr object. Mirrors mod_results.R's structure but every tab
# needs a chosen time first (plot(tselect=) for the Plot tab,
# spattemp.slice() to get the actual pixel image for Table/Download).

spattemp_results_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::uiOutput(ns("time_ui")),
    shiny::fluidRow(
      shiny::column(6,
        shiny::radioButtons(ns("type"), "Surface", inline = TRUE,
                             choices = c("Joint (unconditional)" = "joint", "Conditional on time" = "conditional"))
      ),
      shiny::column(6, shiny::uiOutput(ns("tol_ui")))
    ),
    shiny::tabsetPanel(
      # ic-results-frame caps Plot/Table to a sensible viewing/reading width
      # instead of stretching a fixed-size plot image or a short summary
      # table across whatever width sidebarLayout's mainPanel leaves (which
      # grows with the browser window). Map is left unconstrained -- more
      # width is genuinely useful for a pannable basemap.
      shiny::tabPanel("Plot", shiny::div(
        class = "ic-result-pane ic-results-frame",
        shinycssloaders::withSpinner(shiny::imageOutput(ns("plot"), height = 550))
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
        shiny::downloadButton(ns("dl_png"), "Plot at this time (PNG)"),
        shiny::downloadButton(ns("dl_pdf"), "Plot at this time (PDF)"),
        shiny::downloadButton(ns("dl_grid"), "Surface grid at this time (CSV)"),
        shiny::downloadButton(ns("dl_script"), "Reproducible R script"),
      )
    )
  )
}

spattemp_results_server <- function(id, result, meta) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$time_ui <- shiny::renderUI({
      r <- result()
      shiny::validate(shiny::need(r, "Run an analysis to see the time slider."))
      # Deliberately bounded to range(tgrid), not the wider $tlim: verified
      # that spattemp.slice() (used for Table/Download) warns and can return
      # unusable output right at the exact $tlim edges, since sparr pads its
      # internal grid inward from the requested tlim for kernel edge
      # handling -- tgrid's own extent is the actually-safe sliceable range.
      tg <- spattemp_tgrid(r)
      shiny::fluidRow(
        shiny::column(
          9, shiny::sliderInput(ns("time_sel"), "Time", min = min(tg), max = max(tg), value = min(tg), step = diff(range(tg)) / 200, width = "100%")
        ),
        shiny::column(3, shiny::br(), shiny::actionButton(ns("play_pause"), "Play"))
      )
    })
    shiny::outputOptions(output, "time_ui", suspendWhenHidden = FALSE)

    output$tol_ui <- shiny::renderUI({
      r <- result()
      if (is.null(r) || !inherits(r, "rrst") || is.null(r$P)) return(NULL)
      shiny::tagList(
        shiny::sliderInput(ns("tol_levels"), "Tolerance contour significance level(s)",
                            min = 0.001, max = 0.2, value = c(0.05), step = 0.001),
        shiny::checkboxInput(ns("show_lower"), "Also outline areas of lower significance (dotted)", value = FALSE)
      )
    })
    shiny::outputOptions(output, "tol_ui", suspendWhenHidden = FALSE)

    tol_levels <- shiny::reactive(if (!is.null(input$tol_levels)) input$tol_levels else c(0.05, 0.01))
    show_lower <- shiny::reactive(isTRUE(input$show_lower))
    surf_type <- shiny::reactive(input$type %||% "joint")

    # The slider fires an input update on every pixel of drag, not just on
    # release -- and Plot/Map/Table (see the suspendWhenHidden note below)
    # all three re-render on *every* one of those, even for hidden tabs.
    # Undebounced, that queues far more renders than the app can keep up
    # with while dragging, so playback visibly lags behind the handle.
    # Debouncing collapses a fast drag down to one render per pause.
    time_sel_raw <- shiny::reactive({
      shiny::req(input$time_sel)
      input$time_sel
    })
    time_sel <- shiny::debounce(time_sel_raw, millis = 150)

    # Shared slice so Map and Table (and the Plot tab's lower-significance
    # overlay) don't each separately re-interpolate the same time point.
    current_slice <- shiny::reactive({
      r <- result()
      shiny::req(r, time_sel())
      spattemp_slice_im(r, time_sel(), type = surf_type())
    })

    # Frame cache for Play: every tgrid time point's Plot-tab PNG, pre-rendered
    # to a temp file up front so that during playback each advance is just an
    # instant file swap, not a fresh spatstat/contour render racing the
    # invalidateLater timer below (which is what used to let the slider move
    # on to the next time before the previous frame had actually drawn).
    # Keyed on whatever affects the pixels (a fresh analysis run, surface
    # type, tolerance levels, lower-significance toggle) so a stale cache is
    # never served after any of those change.
    cache_gen <- shiny::reactiveVal(0)
    shiny::observeEvent(result(), cache_gen(cache_gen() + 1), ignoreNULL = FALSE)
    frame_cache <- shiny::reactiveVal(NULL)
    cache_key <- shiny::reactive(list(gen = cache_gen(), type = surf_type(), tol = tol_levels(), lower = show_lower()))
    last_live_tmp <- NULL

    session$onSessionEnded(function() {
      fc <- shiny::isolate(frame_cache())
      if (!is.null(fc)) unlink(unlist(fc$paths))
      if (!is.null(last_live_tmp)) unlink(last_live_tmp)
    })

    # Play/Pause: auto-advance the time slider through the object's own
    # computed time grid every 400ms while playing, satisfying "animation"
    # without rendering to a GIF. req(playing()) halts the invalidateLater
    # chain as soon as the button is toggled back off. Starting playback
    # first (synchronously, with a progress bar) renders every frame the
    # animation will need, so the 400ms loop below never has to wait on one.
    playing <- shiny::reactiveVal(FALSE)
    shiny::observeEvent(input$play_pause, {
      if (!playing()) {
        r <- result()
        shiny::req(r)
        key <- cache_key()
        fc <- frame_cache()
        if (is.null(fc) || !identical(fc$key, key)) {
          if (!is.null(fc)) unlink(unlist(fc$paths))
          tg <- spattemp_tgrid(r)
          paths <- stats::setNames(vector("list", length(tg)), as.character(tg))
          shiny::withProgress(message = "Preparing animation frames…", value = 0, {
            for (i in seq_along(tg)) {
              tmp <- tempfile(fileext = ".png")
              write_spattemp_plot_png(tmp, r, tg[i], surf_type(), tol_levels(), show_lower(), width = 900, height = 550)
              paths[[i]] <- tmp
              shiny::incProgress(1 / length(tg), detail = sprintf("Frame %d of %d", i, length(tg)))
            }
          })
          frame_cache(list(key = key, paths = paths))
        }
      }
      playing(!playing())
    })
    shiny::observe(shiny::updateActionButton(session, "play_pause", label = if (playing()) "Pause" else "Play"))
    shiny::observe({
      shiny::req(playing())
      r <- result()
      shiny::req(r)
      shiny::invalidateLater(400, session)
      # isolate(): input$time_sel must NOT be a reactive dependency of this
      # observer. It used to be (via a direct read below), which meant every
      # updateSliderInput() call here -- once the client echoed it back as a
      # fresh input$time_sel -- re-triggered this same observer immediately,
      # on top of the 400ms invalidateLater timer. The two triggers raced,
      # each computing "next" from whatever input$time_sel happened to be at
      # that instant, so the slider jumped around non-monotonically far
      # faster than once per 400ms instead of stepping through the grid in
      # order -- which is what actually made playback outrun rendering.
      cur <- shiny::isolate(input$time_sel)
      shiny::req(cur)
      tg <- spattemp_tgrid(r)
      nxt_idx <- which(tg > cur + .Machine$double.eps^0.5)[1]
      nxt <- if (is.na(nxt_idx)) tg[1] else tg[nxt_idx]
      shiny::updateSliderInput(session, "time_sel", value = nxt)
    })

    output$plot <- shiny::renderImage({
      r <- result()
      shiny::validate(shiny::need(r, "Run an analysis to see the plot."))
      ts <- time_sel()
      shiny::validate(shiny::need(!is.null(ts), "Pick a time."))

      # Serve the pre-rendered frame if Play already built one for this exact
      # time/settings combo (also lets a manual drag that happens to land on
      # a cached tgrid tick render instantly); otherwise fall back to a
      # one-off live render, replacing whichever one-off file preceded it.
      fc <- frame_cache()
      cached_path <- NULL
      if (!is.null(fc) && identical(fc$key, cache_key())) {
        nums <- as.numeric(names(fc$paths))
        idx <- which.min(abs(nums - ts))
        if (length(idx) && abs(nums[idx] - ts) < 1e-6) cached_path <- fc$paths[[idx]]
      }
      if (!is.null(cached_path)) {
        list(src = cached_path, contentType = "image/png", width = 900, height = 550)
      } else {
        if (!is.null(last_live_tmp)) unlink(last_live_tmp)
        tmp <- tempfile(fileext = ".png")
        write_spattemp_plot_png(tmp, r, ts, surf_type(), tol_levels(), show_lower(), width = 900, height = 550)
        last_live_tmp <<- tmp
        list(src = tmp, contentType = "image/png", width = 900, height = 550)
      }
    }, deleteFile = FALSE)

    output$map <- leaflet::renderLeaflet({
      r <- result()
      shiny::validate(shiny::need(r, "Run an analysis to see the map."))
      shiny::validate(shiny::need(!is.null(time_sel()), ""))
      m <- meta()
      if (!isTRUE(m$is_geo)) {
        return(
          leaflet::leaflet() |>
            leaflet::addProviderTiles(leaflet::providers$CartoDB.Positron) |>
            leaflet::setView(0, 0, zoom = 1)
        )
      }
      sl <- current_slice()
      ras <- im_to_wgs84_raster(sl$z, m$epsg)
      pal <- leaflet::colorNumeric(sparr_colour_ramp()(256), terra::values(ras), na.color = "transparent")
      leaflet::leaflet() |>
        leaflet::addProviderTiles(leaflet::providers$CartoDB.Positron) |>
        leaflet::addRasterImage(ras, colors = pal, opacity = 0.75) |>
        leaflet::addLegend(pal = pal, values = terra::values(ras), title = "Surface value")
    })

    output$table <- DT::renderDT({
      r <- result()
      shiny::validate(shiny::need(r, "Run an analysis to see summary statistics."))
      shiny::validate(shiny::need(!is.null(time_sel()), ""))
      sl <- current_slice()
      imobj <- sl$z
      all_vals <- imobj$v[!is.na(imobj$v)]
      vals <- all_vals[is.finite(all_vals)]
      n_infinite <- sum(!is.finite(all_vals))
      stats <- data.frame(
        Statistic = c("Time", "Surface", "Grid resolution", "Min", "Max", "Mean", "Median"),
        Value = c(
          signif(time_sel(), 6), surf_type(), paste0(imobj$dim[2], " x ", imobj$dim[1]),
          signif(min(vals), 4), signif(max(vals), 4), signif(mean(vals), 4), signif(stats::median(vals), 4)
        )
      )
      if (n_infinite > 0) {
        stats <- rbind(stats, data.frame(Statistic = "Non-finite pixels excluded above (±Inf)", Value = n_infinite))
      }
      if (!is.null(sl$P)) {
        for (lev in tol_levels()) {
          stats <- rbind(stats, data.frame(
            Statistic = paste0("Pixels significant at p<", lev),
            Value = sum(sl$P$v < lev, na.rm = TRUE)
          ))
        }
      }
      # An explicit width, not "100%" of .ic-results-frame: this is a
      # 2-column, half-dozen-row summary table, not a data grid -- letting
      # it fill the frame would just spread "Statistic" and "Value" apart
      # with a wide gap of nothing between them.
      DT::datatable(stats, options = list(dom = "t", paging = FALSE), rownames = FALSE, width = "480px")
    })
    # See mod_results.R for why these stay at suspendWhenHidden=FALSE: the
    # default (TRUE) hits a visibility-detection race in a hidden sub-tab
    # where the output can stay stuck "recalculating" forever.
    for (nm in c("plot", "map", "table")) shiny::outputOptions(output, nm, suspendWhenHidden = FALSE)

    output$dl_png <- shiny::downloadHandler(
      filename = function() paste0("sparr_spattemp_t", round(input$time_sel %||% 0), ".png"),
      content = function(file) write_spattemp_plot_png(file, result(), input$time_sel, surf_type(), tol_levels(), show_lower())
    )
    output$dl_pdf <- shiny::downloadHandler(
      filename = function() paste0("sparr_spattemp_t", round(input$time_sel %||% 0), ".pdf"),
      content = function(file) write_spattemp_plot_pdf(file, result(), input$time_sel, surf_type(), tol_levels(), show_lower())
    )
    output$dl_grid <- shiny::downloadHandler(
      filename = function() paste0("sparr_spattemp_grid_t", round(input$time_sel %||% 0), ".csv"),
      content = function(file) write_grid_csv(file, spattemp_slice_im(result(), input$time_sel, surf_type())$z)
    )
    output$dl_script <- shiny::downloadHandler(
      filename = function() "sparr_spattemp_analysis.R",
      content = function(file) writeLines(generate_repro_script(meta()$repro_params), file)
    )
    for (nm in c("dl_png", "dl_pdf", "dl_grid", "dl_script")) shiny::outputOptions(output, nm, suspendWhenHidden = FALSE)
  })
}
