# Spatiotemporal density module: UI/server for sparr::spattemp.density() on a
# single population evolving over time. `ppp_data` is the reactive returned
# by upload_server(); needs a time column mapped in the Upload tab.

spattemp_density_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::div(
      class = "ic-group",
      shiny::uiOutput(ns("subset_ui")),
      shiny::radioButtons(
        ns("bw_method"), "Bandwidths (spatial h, temporal lambda)",
        choices = c("Automatic (oversmoothing + Sheather-Jones)" = "auto",
                    "Normal scale (NS.spattemp)" = "NS",
                    "Least-squares CV (LSCV.spattemp)" = "LSCV",
                    "Likelihood CV (LIK.spattemp)" = "LIK",
                    "Fixed values" = "fixed"),
        selected = "auto"
      ),
      shiny::conditionalPanel(
        condition = "input.bw_method == 'fixed'", ns = ns,
        shiny::numericInput(ns("h_fixed"), "Spatial bandwidth (h)", value = 1, min = 0.0001),
        shiny::numericInput(ns("lambda_fixed"), "Temporal bandwidth (lambda)", value = 1, min = 0.0001)
      )
    ),
    shiny::div(
      class = "ic-group",
      shiny::selectInput(ns("edge"), "Edge correction (spatial & temporal)",
                          choices = c("Uniform" = "uniform", "None" = "none")),
      shiny::numericInput(ns("sres"), "Spatial grid resolution", value = 64, min = 16, max = 256, step = 16),
      shiny::numericInput(ns("tres"), "Time grid resolution (number of time steps)", value = 50, min = 5, max = 300, step = 5),
      shiny::helpText(
        "Time steps are spread evenly across the data's time range by default. A dataset spanning years ",
        "at daily resolution can have thousands of possible time points -- tres caps how many are actually ",
        "computed, since each one is a full spatial density surface. 50 is a good starting point; raise it ",
        "for a smoother scrub/animation at the cost of longer computation."
      )
    ),
    shiny::div(
      class = "ic-group",
      shiny::helpText("Cross-validation and fine grids can take from several seconds to a minute or more."),
      shiny::actionButton(ns("run"), "Run spatiotemporal density estimation", class = "btn-primary"),
      shiny::hr(),
      shiny::verbatimTextOutput(ns("status"))
    )
  )
}

spattemp_density_server <- function(id, ppp_data) {
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
      shiny::validate(shiny::need(d$timecol, "Map a time column in the Upload tab, then rebuild the point pattern."))
      pp <- d$pp
      tt <- d$tt

      subset_level <- if (!is.null(d$markcol)) (input$subset_level %||% "__all__") else "__all__"
      if (!is.null(d$markcol) && !identical(subset_level, "__all__")) {
        keep_idx <- spatstat.geom::marks(pp) == subset_level
        pp <- pp[keep_idx]
        tt <- tt[keep_idx]
      }

      has_time <- !is.na(tt)
      n_missing_time <- sum(!has_time)
      shiny::validate(shiny::need(sum(has_time) >= 2, "Fewer than 2 points in this population have a time value."))
      pp <- pp[has_time]
      tt <- tt[has_time]

      shiny::withProgress(message = "Running spatiotemporal density estimation…", value = 0.15, {
        hlam <- switch(input$bw_method,
          auto  = c(h = NA_real_, lambda = NA_real_),
          fixed = c(h = input$h_fixed, lambda = input$lambda_fixed),
          NS    = sparr::NS.spattemp(pp, tt = tt),
          LSCV  = sparr::LSCV.spattemp(pp, tt = tt, sres = input$sres, tres = input$sres, verbose = FALSE),
          LIK   = sparr::LIK.spattemp(pp, tt = tt, verbose = FALSE)
        )
        shiny::incProgress(0.3, detail = "Estimating spatiotemporal surface")
        res <- sparr::spattemp.density(
          pp, tt = tt,
          h = if (is.na(hlam[["h"]])) NULL else unname(hlam[["h"]]),
          lambda = if (is.na(hlam[["lambda"]])) NULL else unname(hlam[["lambda"]]),
          sedge = input$edge, tedge = input$edge,
          sres = input$sres, tres = input$tres, verbose = FALSE
        )
        shiny::incProgress(0.4)
      })
      list(result = res, h = res$h, lambda = res$lambda, subset_level = subset_level, n_missing_time = n_missing_time)
    })

    output$status <- shiny::renderPrint({
      tryCatch({
        c <- computed()
        cat("Population smoothed:", if (identical(c$subset_level, "__all__")) "all groups pooled" else paste0("'", c$subset_level, "'"),
            " (", spatstat.geom::npoints(c$result$pp), " points with a time value)\n", sep = "")
        if (c$n_missing_time > 0) cat(c$n_missing_time, "point(s) excluded for missing a time value\n")
        cat("h (spatial) used:", signif(c$h, 4), "\n")
        cat("lambda (temporal) used:", signif(c$lambda, 4), "\n")
        cat("Time range:", paste(signif(c$result$tlim, 6), collapse = " to "), "\n")
        cat("Time steps computed:", length(c$result$tgrid), "\n")
      }, error = function(e) cat("Error:", conditionMessage(e)))
    })
    shiny::outputOptions(output, "status", suspendWhenHidden = FALSE)

    result <- shiny::reactive(tryCatch(computed()$result, error = function(e) NULL))

    meta <- shiny::reactive({
      d <- ppp_data()
      shiny::req(d)
      list(
        analysis = "spattemp_density", is_geo = isTRUE(d$is_geo), epsg = d$epsg,
        repro_params = list(
          analysis = "spattemp_density", xcol = d$xcol, ycol = d$ycol, markcol = d$markcol,
          timecol = d$timecol, time_is_date = isTRUE(d$time_is_date),
          subset_level = tryCatch(computed()$subset_level, error = function(e) "__all__"),
          window_method = d$window_method, pad_pct = d$pad_pct,
          boundary_buffer = d$boundary_buffer, is_geo = isTRUE(d$is_geo),
          h = tryCatch(signif(computed()$h, 4), error = function(e) "NA  # run failed"),
          lambda = tryCatch(signif(computed()$lambda, 4), error = function(e) "NA  # run failed"),
          edge = input$edge, sres = input$sres, tres = input$tres, weightcol = d$weightcol
        )
      )
    })

    list(result = result, meta = meta)
  })
}
