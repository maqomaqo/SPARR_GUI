# Risk module: UI/server for sparr::risk() case/control relative risk with
# tolerance contours. Requires the uploaded ppp to have a 2-level group mark.

risk_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::div(
      class = "ic-group",
      shiny::uiOutput(ns("case_ui")),
      shiny::radioButtons(
        ns("h0_method"), "Global bandwidth (h0)",
        choices = c("Oversmoothing, geometric (OS)" = "OS", "Normal scale, geometric (NS)" = "NS",
                    "Jointly-optimal LSCV" = "LSCVrisk", "Fixed value" = "fixed"),
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
        shiny::selectInput(ns("pilot_symmetry"), "Adaptive pilot symmetry",
                            choices = c("None (asymmetric)" = "none", "Case" = "f",
                                        "Control" = "g", "Pooled" = "pooled"))
      ),
      shiny::numericInput(ns("resolution"), "Grid resolution", value = 128, min = 16, max = 512, step = 16),
      shiny::checkboxInput(ns("tolerate"), "Compute tolerance contours (asymptotic p-value surface)", value = TRUE)
    ),
    shiny::div(
      class = "ic-group",
      shiny::helpText("LSCV and adaptive smoothing can take from several seconds to a few minutes, ",
                       "especially above resolution 128."),
      shiny::actionButton(ns("run"), "Run risk estimation", class = "btn-primary"),
      shiny::hr(),
      shiny::verbatimTextOutput(ns("status"))
    )
  )
}

risk_server <- function(id, ppp_data) {
  shiny::moduleServer(id, function(input, output, session) {

    output$case_ui <- shiny::renderUI({
      d <- ppp_data()
      if (is.null(d) || is.null(d$markcol)) {
        return(shiny::helpText("Select a group column with exactly 2 levels in the Upload tab first."))
      }
      lv <- tryCatch(validate_case_control_marks(spatstat.geom::marks(d$pp)), error = function(e) NULL)
      if (is.null(lv)) return(shiny::helpText("Group column must have exactly 2 distinct values."))
      shiny::selectInput(
        session$ns("case_level"), "Which group is the 'case' (numerator)?",
        choices = lv, selected = lv[1]
      )
    })
    shiny::outputOptions(output, "case_ui", suspendWhenHidden = FALSE)

    computed <- shiny::eventReactive(input$run, {
      d <- shiny::req(ppp_data())
      pp <- d$pp
      shiny::validate(shiny::need(
        d$markcol, "Select a group column with exactly 2 levels in the Upload tab, then rebuild the point pattern."
      ))
      lv <- validate_case_control_marks(spatstat.geom::marks(pp))
      case_level <- if (!is.null(input$case_level) && input$case_level %in% lv) input$case_level else lv[1]
      control_level <- setdiff(lv, case_level)
      # sparr::risk() treats the FIRST factor level of a dichotomous-marked
      # ppp as the case group and the second as control -- relevel so the
      # user's choice (not just alphabetical order) decides which is which.
      spatstat.geom::marks(pp) <- factor(spatstat.geom::marks(pp), levels = c(case_level, control_level))

      shiny::withProgress(message = "Running risk estimation…", value = 0.15, {
        h0 <- switch(input$h0_method,
          fixed    = input$h0_fixed,
          OS       = sparr::OS(pp, nstar = "geometric"),
          NS       = sparr::NS(pp, nstar = "geometric"),
          LSCVrisk = sparr::LSCV.risk(pp, verbose = FALSE)
        )
        shiny::incProgress(0.3, detail = "Estimating risk surface")
        res <- sparr::risk(
          pp, h0 = h0, adapt = isTRUE(input$adapt),
          pilot.symmetry = input$pilot_symmetry %||% "none",
          resolution = input$resolution, tolerate = isTRUE(input$tolerate),
          verbose = FALSE
        )
        shiny::incProgress(0.4)
      })
      list(result = res, h0 = h0, case_level = case_level, control_level = control_level)
    })

    output$status <- shiny::renderPrint({
      tryCatch({
        c <- computed()
        cat("Case group: '", c$case_level, "' (", spatstat.geom::npoints(c$result$f$pp), " points)\n", sep = "")
        cat("Control group: '", c$control_level, "' (", spatstat.geom::npoints(c$result$g$pp), " points)\n", sep = "")
        cat("h0 used:", paste(signif(c$h0, 4), collapse = ", "), "\n")
      }, error = function(e) cat("Error:", conditionMessage(e)))
    })
    # Never suspend: the sidebar status text lives outside the visible-tab
    # detection that gates the results tabset, so leaving suspendWhenHidden
    # on its default (TRUE) can leave this stuck waiting for a visibility
    # signal that's already true, and the run never appears to finish.
    shiny::outputOptions(output, "status", suspendWhenHidden = FALSE)

    result <- shiny::reactive(tryCatch(computed()$result, error = function(e) NULL))

    meta <- shiny::reactive({
      d <- ppp_data()
      shiny::req(d)
      list(
        analysis = "risk", is_geo = isTRUE(d$is_geo), epsg = d$epsg,
        repro_params = list(
          analysis = "risk", xcol = d$xcol, ycol = d$ycol, markcol = d$markcol,
          case_level = tryCatch(computed()$case_level, error = function(e) NULL),
          control_level = tryCatch(computed()$control_level, error = function(e) NULL),
          window_method = d$window_method, pad_pct = d$pad_pct,
          boundary_buffer = d$boundary_buffer, is_geo = isTRUE(d$is_geo),
          h0 = tryCatch(paste(signif(computed()$h0, 4), collapse = ", "), error = function(e) "NA  # run failed"),
          adapt = isTRUE(input$adapt), resolution = input$resolution,
          tolerate = isTRUE(input$tolerate), weightcol = d$weightcol
        )
      )
    })

    list(result = result, meta = meta)
  })
}
