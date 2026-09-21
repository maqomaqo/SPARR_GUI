# Density module: UI/server for sparr::bivariate.density() on a single
# population. `ppp_data` is the reactive returned by upload_server().

density_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::div(
      class = "ic-group",
      shiny::uiOutput(ns("subset_ui")),
      shiny::radioButtons(
        ns("h0_method"), "Global bandwidth (h0)",
        choices = c("Oversmoothing (OS)" = "OS", "Normal scale (NS)" = "NS",
                    "Least-squares CV (LSCV)" = "LSCV", "Likelihood CV (LIK)" = "LIK",
                    "Fixed value" = "fixed"),
        selected = "OS"
      ),
      shiny::conditionalPanel(
        condition = "input.h0_method == 'fixed'", ns = ns,
        shiny::numericInput(ns("h0_fixed"), "h0 value", value = 1, min = 0.0001)
      )
    ),
    shiny::div(
      class = "ic-group",
      shiny::checkboxInput(ns("adapt"), "Adaptive smoothing", value = FALSE),
      shiny::conditionalPanel(
        condition = "input.adapt == true", ns = ns,
        shiny::numericInput(ns("hp"), "Pilot bandwidth (hp, blank = same as h0)", value = NA)
      ),
      shiny::selectInput(ns("edge"), "Edge correction",
                          choices = c("Uniform" = "uniform", "Diggle" = "diggle", "None" = "none")),
      shiny::numericInput(ns("resolution"), "Grid resolution", value = 128, min = 16, max = 512, step = 16)
    ),
    shiny::div(
      class = "ic-group",
      shiny::helpText("LSCV/LIK and adaptive smoothing can take from several seconds to a few minutes, ",
                       "especially above resolution 128."),
      shiny::actionButton(ns("run"), "Run density estimation", class = "btn-primary"),
      shiny::hr(),
      shiny::verbatimTextOutput(ns("status"))
    )
  )
}

density_server <- function(id, ppp_data) {
  shiny::moduleServer(id, function(input, output, session) {

    output$subset_ui <- shiny::renderUI({
      d <- ppp_data()
      if (is.null(d) || is.null(d$markcol)) return(NULL)
      lv <- levels(factor(spatstat.geom::marks(d$pp)))
      shiny::selectInput(
        session$ns("subset_level"), "Population to smooth",
        choices = c("All groups pooled" = "__all__", stats::setNames(lv, lv)),
        selected = "__all__"
      )
    })
    shiny::outputOptions(output, "subset_ui", suspendWhenHidden = FALSE)

    computed <- shiny::eventReactive(input$run, {
      d <- shiny::req(ppp_data())
      pp <- d$pp
      weights <- d$weights
      # input$subset_level can retain a stale value from a previous upload
      # that had a group column even after this one doesn't (the selectize
      # widget disappears, but Shiny keeps the last value it sent) -- only
      # trust it when there's actually a group column to subset by.
      subset_level <- if (!is.null(d$markcol)) (input$subset_level %||% "__all__") else "__all__"
      if (!is.null(d$markcol) && !identical(subset_level, "__all__")) {
        # Subset by logical index (rather than split.ppp, which reorders by
        # group) so `weights` -- aligned to `pp`'s original point order --
        # stays in lockstep with which points survive the filter.
        keep_idx <- spatstat.geom::marks(pp) == subset_level
        pp <- pp[keep_idx]
        if (!is.null(weights)) weights <- weights[keep_idx]
      }

      shiny::withProgress(message = "Running density estimation…", value = 0.2, {
        h0 <- switch(input$h0_method,
          fixed = input$h0_fixed,
          OS    = sparr::OS(pp),
          NS    = sparr::NS(pp),
          LSCV  = sparr::LSCV.density(pp, verbose = FALSE),
          LIK   = sparr::LIK.density(pp, verbose = FALSE)
        )
        shiny::incProgress(0.4, detail = "Estimating surface")
        hp <- if (isTRUE(input$adapt) && !is.na(input$hp)) input$hp else NULL
        res <- sparr::bivariate.density(
          pp, h0 = h0, hp = hp, adapt = isTRUE(input$adapt),
          resolution = input$resolution, edge = input$edge,
          weights = weights, verbose = FALSE
        )
        shiny::incProgress(0.4)
      })
      list(result = res, h0 = h0, subset_level = subset_level)
    })

    output$status <- shiny::renderPrint({
      tryCatch({
        c <- computed()
        cat("Population smoothed:", if (identical(c$subset_level, "__all__")) "all groups pooled" else paste0("'", c$subset_level, "'"),
            " (", spatstat.geom::npoints(c$result$pp), " points)\n", sep = "")
        cat("h0 used:", signif(c$h0, 4), "\n")
        if (isTRUE(input$adapt)) cat("hp used:", signif(c$result$hp, 4), "\n")
      }, error = function(e) cat("Error:", conditionMessage(e)))
    })
    # See comment in mod_risk.R: avoid this output getting stuck waiting on
    # a visibility signal that's already true.
    shiny::outputOptions(output, "status", suspendWhenHidden = FALSE)

    result <- shiny::reactive(tryCatch(computed()$result, error = function(e) NULL))

    meta <- shiny::reactive({
      d <- ppp_data()
      shiny::req(d)
      list(
        analysis = "density", is_geo = isTRUE(d$is_geo), epsg = d$epsg,
        repro_params = list(
          analysis = "density", xcol = d$xcol, ycol = d$ycol, markcol = d$markcol,
          subset_level = tryCatch(computed()$subset_level, error = function(e) "__all__"),
          window_method = d$window_method, pad_pct = d$pad_pct,
          boundary_buffer = d$boundary_buffer, is_geo = isTRUE(d$is_geo),
          h0 = tryCatch(signif(computed()$h0, 4), error = function(e) "NA  # run failed"),
          adapt = isTRUE(input$adapt),
          hp = if (isTRUE(input$adapt) && !is.na(input$hp)) input$hp else NULL,
          resolution = input$resolution, edge = input$edge, weightcol = d$weightcol
        )
      )
    })

    list(result = result, meta = meta)
  })
}
