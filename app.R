# sparr GUI — entry point. Run with `Rscript run_app.R` or shiny::runApp().

options(shiny.maxRequestSize = 30 * 1024^2)  # 30 MB CSV upload limit

# macOS: the default "quartz" PNG device backend requires a live window-server
# session and fails ("invalid quartz() device size") when R runs headless via
# Rscript, which is how this app is launched. Cairo doesn't have that
# dependency and is the standard fix for headless/server Shiny on macOS.
if (Sys.info()[["sysname"]] == "Darwin" && capabilities("cairo")) {
  options(bitmapType = "cairo")
}

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f, local = FALSE)

ui <- bslib::page_navbar(
  title = "sparr — spatial relative risk explorer",
  theme = bslib::bs_theme(version = 5, primary = "#2C5F8A"),
  header = shiny::tags$head(shiny::tags$link(rel = "stylesheet", type = "text/css", href = "layout.css")),
  # bslib pages default to a fill/flexbox layout; combined with
  # sidebarLayout()/tabsetPanel() nesting that isn't fill-aware, that can
  # collapse plotOutput's measured width to 0 on a freshly-shown tab. Plain
  # block layout avoids it and this app has no need for viewport-filling.
  fillable = FALSE,
  bslib::nav_panel(
    "1. Data & window",
    shiny::fluidRow(
      shiny::column(5, upload_ui("upload")),
      shiny::column(7, shiny::div(class = "ic-intro-copy", shiny::markdown(paste(
        "**How this works:** upload a CSV of point locations (or pick an example), map the",
        "columns, choose how the study region is defined, then click **Build point pattern**.",
        "Once built, go to the **Density** tab to estimate a smoothed intensity surface for a",
        "single set of points, or **Risk** tab (needs a 2-level group column) for a case/control",
        "relative risk surface with tolerance contours."
      ))))
    )
  ),
  bslib::nav_panel(
    "2. Density",
    shiny::sidebarLayout(
      shiny::sidebarPanel(density_ui("density"), width = 3),
      shiny::mainPanel(results_ui("density_results"), width = 9)
    )
  ),
  bslib::nav_panel(
    "3. Risk",
    shiny::sidebarLayout(
      shiny::sidebarPanel(risk_ui("risk"), width = 3),
      shiny::mainPanel(results_ui("risk_results"), width = 9)
    )
  ),
  bslib::nav_panel(
    "4. Spatiotemporal density",
    shiny::sidebarLayout(
      shiny::sidebarPanel(spattemp_density_ui("spattemp_density"), width = 3),
      shiny::mainPanel(spattemp_results_ui("spattemp_density_results"), width = 9)
    )
  ),
  bslib::nav_panel(
    "5. Spatiotemporal risk",
    shiny::sidebarLayout(
      shiny::sidebarPanel(spattemp_risk_ui("spattemp_risk"), width = 3),
      shiny::mainPanel(spattemp_results_ui("spattemp_risk_results"), width = 9)
    )
  )
)

server <- function(input, output, session) {
  ppp_data <- upload_server("upload")

  den <- density_server("density", ppp_data)
  results_server("density_results", den$result, den$meta)

  rk <- risk_server("risk", ppp_data)
  results_server("risk_results", rk$result, rk$meta)

  std <- spattemp_density_server("spattemp_density", ppp_data)
  spattemp_results_server("spattemp_density_results", std$result, std$meta)

  str <- spattemp_risk_server("spattemp_risk", ppp_data)
  spattemp_results_server("spattemp_risk_results", str$result, str$meta)
}

shiny::shinyApp(ui, server)
