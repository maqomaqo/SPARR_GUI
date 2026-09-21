# Spatiotemporal risk module: UI/server for sparr::spattemp.risk() -- a
# time-varying case density (spattemp.density) against a control group.
# v1 treats the control as static over time (bivariate.density), matching
# sparr's own documented pattern (used in both its burk and fmd examples,
# since those datasets' control/at-risk populations don't carry meaningful
# event times) -- a time-varying control is a possible future enhancement.

spattemp_risk_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::div(
      class = "ic-group",
      shiny::uiOutput(ns("case_ui")),
      shiny::radioButtons(
        ns("case_bw_method"), "Case bandwidths (spatial h, temporal lambda)",
        choices = c("Automatic (oversmoothing + Sheather-Jones)" = "auto",
                    "Normal scale (NS.spattemp)" = "NS",
                    "Least-squares CV (LSCV.spattemp)" = "LSCV",
                    "Likelihood CV (LIK.spattemp)" = "LIK",
                    "Fixed values" = "fixed"),
        selected = "auto"
      ),
      shiny::conditionalPanel(
        condition = "input.case_bw_method == 'fixed'", ns = ns,
        shiny::numericInput(ns("h_fixed"), "Case spatial bandwidth (h)", value = 1, min = 0.0001),
        shiny::numericInput(ns("lambda_fixed"), "Case temporal bandwidth (lambda)", value = 1, min = 0.0001)
      )
    ),
    shiny::div(
      class = "ic-group",
      shiny::radioButtons(
        ns("control_bw_method"), "Control bandwidth (static density, h0)",
        choices = c("Oversmoothing (OS)" = "OS", "Normal scale (NS)" = "NS", "Fixed value" = "fixed"),
        selected = "OS"
      ),
      shiny::conditionalPanel(
        condition = "input.control_bw_method == 'fixed'", ns = ns,
        shiny::numericInput(ns("h0_fixed"), "Control h0 value", value = 1, min = 0.0001)
      )
    ),
    shiny::div(
      class = "ic-group",
      shiny::numericInput(ns("resolution"), "Grid resolution (shared by case and control -- required to match)",
                           value = 64, min = 16, max = 256, step = 16),
      shiny::numericInput(ns("tres"), "Time grid resolution (number of time steps)", value = 50, min = 5, max = 300, step = 5),
      shiny::checkboxInput(ns("tolerate"), "Compute tolerance contours (asymptotic p-value surfaces)", value = TRUE)
    ),
    shiny::div(
      class = "ic-group",
      shiny::helpText(
        "The control group is treated as static over time (its own spatial density, unchanging), matching ",
        "sparr's documented usage -- appropriate when only the case/event group has meaningful timestamps ",
        "(e.g. infection dates), which is the common case for this kind of data."
      ),
      shiny::helpText("Cross-validation and fine grids can take from several seconds to a minute or more."),
      shiny::actionButton(ns("run"), "Run spatiotemporal risk estimation", class = "btn-primary"),
      shiny::hr(),
      shiny::verbatimTextOutput(ns("status"))
    )
  )
}

spattemp_risk_server <- function(id, ppp_data) {
  shiny::moduleServer(id, function(input, output, session) {

    output$case_ui <- shiny::renderUI({
      d <- ppp_data()
      if (is.null(d) || is.null(d$markcol)) {
        return(shiny::helpText("Select a group column with exactly 2 levels in the Upload tab first."))
      }
      lv <- tryCatch(validate_case_control_marks(spatstat.geom::marks(d$pp)), error = function(e) NULL)
      if (is.null(lv)) return(shiny::helpText("Group column must have exactly 2 distinct values."))
      shiny::selectInput(
        session$ns("case_level"), "Which group is the 'case' (time-varying numerator)?",
        choices = lv, selected = lv[1]
      )
    })
    shiny::outputOptions(output, "case_ui", suspendWhenHidden = FALSE)

    computed <- shiny::eventReactive(input$run, {
      d <- shiny::req(ppp_data())
      shiny::validate(shiny::need(d$timecol, "Map a time column in the Upload tab, then rebuild the point pattern."))
      shiny::validate(shiny::need(
        d$markcol, "Select a group column with exactly 2 levels in the Upload tab, then rebuild the point pattern."
      ))
      lv <- validate_case_control_marks(spatstat.geom::marks(d$pp))
      case_level <- if (!is.null(input$case_level) && input$case_level %in% lv) input$case_level else lv[1]
      control_level <- setdiff(lv, case_level)

      marks_vec <- spatstat.geom::marks(d$pp)
      case_idx <- marks_vec == case_level
      control_idx <- marks_vec == control_level

      case_pp <- d$pp[case_idx]
      case_tt <- d$tt[case_idx]
      has_time <- !is.na(case_tt)
      n_missing_time <- sum(!has_time)
      shiny::validate(shiny::need(sum(has_time) >= 2, "Fewer than 2 case points have a time value."))
      case_pp <- case_pp[has_time]
      case_tt <- case_tt[has_time]

      control_pp <- d$pp[control_idx]
      shiny::validate(shiny::need(spatstat.geom::npoints(control_pp) >= 2, "Fewer than 2 control points."))

      shiny::withProgress(message = "Running spatiotemporal risk estimation…", value = 0.1, {
        hlam <- switch(input$case_bw_method,
          auto  = c(h = NA_real_, lambda = NA_real_),
          fixed = c(h = input$h_fixed, lambda = input$lambda_fixed),
          NS    = sparr::NS.spattemp(case_pp, tt = case_tt),
          LSCV  = sparr::LSCV.spattemp(case_pp, tt = case_tt, sres = input$resolution, tres = input$resolution, verbose = FALSE),
          LIK   = sparr::LIK.spattemp(case_pp, tt = case_tt, verbose = FALSE)
        )
        shiny::incProgress(0.2, detail = "Estimating case spatiotemporal density")
        f <- sparr::spattemp.density(
          case_pp, tt = case_tt,
          h = if (is.na(hlam[["h"]])) NULL else unname(hlam[["h"]]),
          lambda = if (is.na(hlam[["lambda"]])) NULL else unname(hlam[["lambda"]]),
          sres = input$resolution, tres = input$tres, verbose = FALSE
        )

        h0 <- switch(input$control_bw_method,
          fixed = input$h0_fixed,
          OS    = sparr::OS(control_pp),
          NS    = sparr::NS(control_pp)
        )
        shiny::incProgress(0.3, detail = "Estimating control density")
        g <- sparr::bivariate.density(control_pp, h0 = h0, resolution = input$resolution, verbose = FALSE)

        shiny::incProgress(0.2, detail = "Computing relative risk")
        res <- sparr::spattemp.risk(f, g, tolerate = isTRUE(input$tolerate), verbose = FALSE)
        shiny::incProgress(0.2)
      })
      list(
        result = res, h = f$h, lambda = f$lambda, h0 = h0,
        case_level = case_level, control_level = control_level, n_missing_time = n_missing_time
      )
    })

    output$status <- shiny::renderPrint({
      tryCatch({
        c <- computed()
        cat("Case group: '", c$case_level, "' (", spatstat.geom::npoints(c$result$f$pp), " points with a time value)\n", sep = "")
        if (c$n_missing_time > 0) cat(c$n_missing_time, "case point(s) excluded for missing a time value\n")
        cat("Control group: '", c$control_level, "' (static, ", spatstat.geom::npoints(c$result$g$pp), " points)\n", sep = "")
        cat("Case h (spatial):", signif(c$h, 4), " | Case lambda (temporal):", signif(c$lambda, 4), "\n")
        cat("Control h0:", signif(c$h0, 4), "\n")
        cat("Time range:", paste(signif(c$result$tlim, 6), collapse = " to "), "\n")
      }, error = function(e) cat("Error:", conditionMessage(e)))
    })
    shiny::outputOptions(output, "status", suspendWhenHidden = FALSE)

    result <- shiny::reactive(tryCatch(computed()$result, error = function(e) NULL))

    meta <- shiny::reactive({
      d <- ppp_data()
      shiny::req(d)
      list(
        analysis = "spattemp_risk", is_geo = isTRUE(d$is_geo), epsg = d$epsg,
        repro_params = list(
          analysis = "spattemp_risk", xcol = d$xcol, ycol = d$ycol, markcol = d$markcol,
          timecol = d$timecol, time_is_date = isTRUE(d$time_is_date),
          case_level = tryCatch(computed()$case_level, error = function(e) NULL),
          control_level = tryCatch(computed()$control_level, error = function(e) NULL),
          window_method = d$window_method, pad_pct = d$pad_pct,
          boundary_buffer = d$boundary_buffer, is_geo = isTRUE(d$is_geo),
          h = tryCatch(signif(computed()$h, 4), error = function(e) "NA  # run failed"),
          lambda = tryCatch(signif(computed()$lambda, 4), error = function(e) "NA  # run failed"),
          h0 = tryCatch(signif(computed()$h0, 4), error = function(e) "NA  # run failed"),
          resolution = input$resolution, tres = input$tres, tolerate = isTRUE(input$tolerate),
          weightcol = d$weightcol
        )
      )
    })

    list(result = result, meta = meta)
  })
}
