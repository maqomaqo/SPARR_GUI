# One-off script that exports sparr's built-in datasets to CSV for the app's
# example dropdown. Uses only spatstat.geom (not sparr itself) to read the
# .rda data files, since data() doesn't require loading sparr's namespace.
# Run once from the repo root: Rscript data/examples/export_examples.R

library(spatstat.geom)
out_dir <- "data/examples"

# --- pbc: already a single dichotomous case/control ppp ---
data(pbc, package = "sparr")
write.csv(
  data.frame(x = pbc$x, y = pbc$y, group = as.character(marks(pbc))),
  file.path(out_dir, "pbc.csv"), row.names = FALSE
)

# --- burk: $cases (time-marked, day of observation) + $controls (spatial-only,
# static "at-risk population" -- no meaningful event time) -> case/control CSV
# with a `day` column (blank for controls) for the spatiotemporal tabs.
data(burk, package = "sparr")
burk_df <- rbind(
  data.frame(x = burk$cases$x, y = burk$cases$y, group = "case", day = marks(burk$cases)),
  data.frame(x = burk$controls$x, y = burk$controls$y, group = "control", day = NA)
)
write.csv(burk_df, file.path(out_dir, "burk.csv"), row.names = FALSE)

# --- fmd: $cases (day-marked, day of infection) + $controls (uninfected farms,
# no event time) -> case/control CSV with a `day` column (blank for controls).
data(fmd, package = "sparr")
fmd_df <- rbind(
  data.frame(x = fmd$cases$x, y = fmd$cases$y, group = "case", day = marks(fmd$cases)),
  data.frame(x = fmd$controls$x, y = fmd$controls$y, group = "control", day = NA)
)
write.csv(fmd_df, file.path(out_dir, "fmd.csv"), row.names = FALSE)

message("Wrote pbc.csv (", nrow(read.csv(file.path(out_dir, "pbc.csv"))), " rows), ",
        "burk.csv (", nrow(burk_df), " rows), fmd.csv (", nrow(fmd_df), " rows) to ", out_dir)
