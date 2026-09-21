# sparr GUI

A local browser app (R Shiny) for [sparr](https://github.com/cran/sparr), the R package for
kernel-smoothed spatial and spatiotemporal relative risk estimation. Upload a CSV of point
locations, configure the same parameters `sparr` exposes in R, and get plots, an optional map,
summary tables, and downloadable output — no R code required.

## What it covers

- **Density**: single-population smoothed intensity/density (`bivariate.density`), with
  bandwidth selection (OS, NS, LSCV, LIK, or a fixed value), fixed or adaptive smoothing, edge
  correction, and optional per-point weights.
- **Risk**: case/control relative risk (`risk`) with asymptotic tolerance contours
  (significance regions of elevated risk), fixed or adaptive smoothing, and an explicit control
  over which group is treated as the case (numerator) — sparr itself just picks the
  alphabetically-first factor level, which this app surfaces and lets you override.
- **Spatiotemporal density**: `spattemp.density()` for a single population's density evolving
  over time, with the same style of bandwidth selection (Automatic, NS.spattemp, LSCV.spattemp,
  LIK.spattemp, or fixed h/lambda).
- **Spatiotemporal risk**: `spattemp.risk()` — a time-varying case density against a control
  group treated as static over time (matching sparr's own documented usage), with asymptotic
  tolerance contours.
- Both spatiotemporal tabs get a **time slider** (scrub to any point in the data's time range,
  not just the computed grid points — sparr interpolates) with a **Joint/Conditional** toggle and
  a **Play/Pause** button that animates through the computed time steps.
- Plot tab (the native sparr/spatstat image+contour plot), Map tab (Leaflet basemap overlay if
  your data is geographic lat/lon), Table tab (summary statistics), and a Download tab (PNG, PDF,
  the underlying surface as a CSV grid, and a generated R script reproducing the exact analysis —
  for the spatiotemporal tabs, all four reflect whichever time/surface is currently selected).
- **Study region**: bounding box, convex hull, or an uploaded boundary polygon — including
  multi-part boundaries (e.g. a mainland plus its islands) with automatic ring reconstruction and
  winding-direction correction (so exports from any GIS tool work, regardless of that tool's
  clockwise/anticlockwise convention). A **boundary buffer** setting relaxes the window outward
  by a fixed distance, to absorb points that land just outside the true boundary due to
  coordinate recording imprecision.

**Not yet included** (see "Roadmap" below): multi-scale density, bootstrap bandwidths,
point-level tolerance classification, Monte Carlo (as opposed to asymptotic) tolerance contours,
and a time-varying (rather than static) control group for spatiotemporal risk.

## 1. Install R (one-time)

You need R installed to run this app — there's no way around that, since `sparr` is an R
package with no equivalent elsewhere. If you don't have it:

1. Download and install R (≥ 4.1, ideally 4.3+) from <https://cran.r-project.org/>.
2. RStudio is optional — not needed to run this app.

**On this machine, R is already installed** (R 4.6.1), so you can skip straight to step 2.

### macOS only: XQuartz is also required

This is a non-obvious catch worth knowing about: `sparr` imports the `misc3d` package, which in
turn imports `tcltk`. On macOS, loading `tcltk` at all requires **XQuartz** (Apple removed X11
from macOS years ago) — even though neither this app nor `sparr` itself actually uses any
tcltk/X11 features. Without XQuartz, `library(sparr)` fails immediately with an X11-related
error, regardless of how the rest of your R setup looks.

Fix (one-time, ~5 minutes):
```bash
brew install --cask xquartz
```
(or download the installer directly from <https://www.xquartz.org/> if you don't use Homebrew).
**Log out and back in (or reboot)** afterwards so the X11 libraries are on the linker path, then
continue below. `setup.R` checks for this and will tell you clearly if it's still missing.

Windows and Linux users are not affected by this (Windows' `tcltk` ships its own Tcl/Tk; Linux
desktop environments normally already have X11).

## 2. Install R packages (one-time)

From this folder:
```bash
Rscript setup.R
```
This installs `sparr` plus every other package the app needs — all from CRAN, all as
pre-built binaries on macOS/Windows (no compiler toolchain needed). Safe to re-run; already
satisfied packages are skipped.

## 3. Run the app

```bash
Rscript run_app.R
```
or on macOS, double-click **run_app.command** in Finder. Either opens the app in your default
browser at a local address (e.g. `http://127.0.0.1:PORT`) — nothing leaves your machine.

## Usage

1. **Data & window tab**: pick one of the bundled example datasets, or upload your own CSV. Map
   which columns are x/y (or longitude/latitude — check the box if so), optionally a group
   column (needed for the Risk tabs; must have exactly 2 distinct values), a weight column, and
   a time column (needed for the Spatiotemporal tabs — dates like `2020-03-14` are converted to
   days automatically if you check the box). Choose how the study region boundary is defined
   (bounding box, convex hull, or your own uploaded polygon) and an optional boundary buffer,
   then **Build point pattern**.
2. **Density tab**: choose a bandwidth method and other parameters, click **Run density
   estimation**. Results appear in Plot / Map / Table / Download.
3. **Risk tab**: same idea, for case/control relative risk with tolerance contours.
4. **Spatiotemporal density / Spatiotemporal risk tabs**: same idea again, but the group column's
   time-varying member(s) need a time value for every point (a static control group, like the
   `burk`/`fmd` examples' at-risk populations, doesn't). After running, use the **Population to
   smooth** (density) or **case-group** (risk) selector alongside the time slider, Joint/
   Conditional toggle, and Play button to explore the result.

Parameter changes don't auto-recompute — some of these calculations (especially cross-validation
bandwidth selection and adaptive smoothing) can take from seconds to minutes, so you explicitly
click "Run" when ready. Cosmetic changes (contour significance level, showing points on the plot)
update instantly from the last computed result.

Example CSVs live in `data/examples/` (`pbc.csv`, `burk.csv`, `fmd.csv` — all exported from
sparr's own built-in datasets, including each case's event day for the spatiotemporal tabs;
see `data/examples/export_examples.R`).

## Project layout

```
app.R                       Shiny entry point
R/mod_upload.R               CSV upload, column mapping, window construction
R/mod_density.R              bivariate.density() controls
R/mod_risk.R                 risk() + tolerance controls
R/mod_results.R              Plot / Map / Table / Download tabs (shared, spatial)
R/mod_spattemp_density.R     spattemp.density() controls
R/mod_spattemp_risk.R        spattemp.risk() controls (time-varying case vs static control)
R/mod_spattemp_results.R     Plot / Map / Table / Download tabs (shared, spatiotemporal --
                              adds the time slider, Joint/Conditional toggle, Play/Pause)
R/utils_ppp.R                data.frame -> ppp/owin builders, CRS handling, multi-part
                              polygon reconstruction, boundary buffering
R/utils_plot.R               plot styling, raster conversion for the map tab
R/utils_export.R             download handlers, reproducible-script generator
data/examples/                bundled example CSVs
setup.R                       one-time dependency installer
run_app.R / .command          launchers
```

## Roadmap (not built yet)

- Multi-scale density (`multiscale.density`)
- Bootstrap bandwidth selectors (`BOOT.density`, `BOOT.spattemp`)
- Monte Carlo tolerance contours and point-level classification (`tol.classify`)
- Time-varying (rather than static) control group for spatiotemporal risk
- Packaging as a standalone desktop app (no visible browser/R)
