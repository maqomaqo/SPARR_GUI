# Launches the app in your default browser. Run with:  Rscript run_app.R
# (or double-click run_app.command on macOS)

# Make sure the working directory is this script's directory, regardless of
# where it was invoked from, so relative paths (R/, data/) resolve correctly.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
if (length(file_arg) == 1) setwd(dirname(normalizePath(file_arg)))

shiny::runApp(".", launch.browser = TRUE)
