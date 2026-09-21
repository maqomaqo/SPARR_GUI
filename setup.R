# Installs every package this app needs. Safe to re-run: already-installed
# packages are skipped. Run with:  Rscript setup.R

if (getRversion() < "4.1.0") {
  stop(
    "R ", getRversion(), " detected. This app needs R >= 4.1.0 ",
    "(4.3+ recommended). Install a current R from https://cran.r-project.org/ ",
    "and re-run this script."
  )
}

# --- macOS preflight: sparr pulls in misc3d, which Imports tcltk. On macOS,
# tcltk needs XQuartz (X11) to load at all -- without it, `library(sparr)`
# fails outright even though sparr's own code never uses tcltk/X11. Check for
# this now with a clear message instead of a cryptic failure mid-install.
if (Sys.info()[["sysname"]] == "Darwin" && !file.exists("/opt/X11/lib/libX11.dylib")) {
  stop(
    "\n\nXQuartz (X11) is not installed.\n\n",
    "sparr depends on the 'misc3d' package, which requires 'tcltk', which on ",
    "macOS requires XQuartz to load -- even though this app never uses any ",
    "tcltk/X11 features directly. Without it, `library(sparr)` fails outright.\n\n",
    "Fix: download and install XQuartz from https://www.xquartz.org/ ",
    "(the .pkg installer), then LOG OUT AND BACK IN (or reboot) so the X11 ",
    "libraries are on the dynamic linker path, then re-run this script.\n"
  )
}

required_pkgs <- c(
  # Shiny app shell
  "shiny", "bslib", "DT", "shinyjs", "shinycssloaders",
  # sparr and its declared dependencies
  "sparr", "spatstat", "spatstat.geom", "spatstat.explore", "spatstat.random",
  "spatstat.utils", "spatstat.univar", "doParallel", "parallel", "foreach", "misc3d",
  # geographic mode (lat/lon basemap overlay)
  "sf", "leaflet", "terra",
  # misc
  "readr", "viridisLite", "zip",
  # optional: faster FFT, sparr falls back to base stats::fft() if missing
  "fftwtools"
)

pkg_loads <- function(pkg) {
  isTRUE(tryCatch({
    suppressPackageStartupMessages(requireNamespace(pkg, quietly = TRUE))
  }, error = function(e) FALSE))
}

missing_pkgs <- required_pkgs[!vapply(required_pkgs, pkg_loads, logical(1))]

if (length(missing_pkgs) == 0) {
  message("All required packages are already installed.")
} else {
  message("Installing ", length(missing_pkgs), " missing package(s): ",
          paste(missing_pkgs, collapse = ", "))
  # Platform-default type (binary on macOS/Windows, source on Linux) —
  # letting R choose avoids "no binary package available" errors on Linux.
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")

  still_missing <- missing_pkgs[!vapply(missing_pkgs, pkg_loads, logical(1))]
  if (length(still_missing) > 0) {
    message("\nThe following packages installed but still fail to load:")
    for (pkg in still_missing) {
      err <- tryCatch({
        suppressPackageStartupMessages(library(pkg, character.only = TRUE))
        NULL
      }, error = function(e) conditionMessage(e))
      message("  - ", pkg, ": ", if (is.null(err)) "(loaded on retry)" else err)
    }
    stop(
      "Failed to load: ", paste(still_missing, collapse = ", "),
      ". See the specific error(s) above."
    )
  }
  message("All packages installed successfully.")
}
