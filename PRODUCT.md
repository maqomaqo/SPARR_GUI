# Product

<!-- impeccable:product-schema 1 -->

## Platform

web

## Users

Mixed / general spatial-statistics audience: GIS analysts, students, and researchers across
domains (public health, ecology, criminology, and similar) who need `sparr`'s kernel-density and
relative-risk estimation methods but don't want to, or can't, write R code directly. Not scoped
to a single discipline — the bundled example datasets happen to be epidemiological (Burkitt's
lymphoma case-control, foot-and-mouth disease outbreak), but the audience is broader than that.

## Product Purpose

A local, browser-based GUI wrapper around the `sparr` R package (kernel-smoothed spatial and
spatiotemporal relative risk / density estimation), which otherwise has no GUI of its own. A user
uploads a CSV of point locations, configures the same parameters `sparr` exposes in R (bandwidth
selection method, edge correction, resolution, tolerance-contour significance levels, etc.), and
gets interactive plots, an optional map, summary tables, and downloadable output — without
writing any R code. Success is a user going from raw point data to a statistically correct
risk/density surface (real tolerance contours, not an ad-hoc heatmap) and either using the
exported output directly or graduating to the generated reproducible R script for further work.

## Positioning

`sparr` itself ships with no GUI — this is the only interactive way to use it. What a generic
"make a heatmap from points" tool could not truthfully copy: this computes actual asymptotic
tolerance contours (statistical significance regions of elevated *or* reduced relative risk, not
just a smoothed density surface implied by color), and it exposes `sparr`'s real bandwidth
machinery (OS / NS / LSCV / LIK / fixed, adaptive smoothing) rather than picking a default for
the user. Results are reproducible science, not a visual approximation.

## Operating Context

Local, single-user, single-session. Launched with `Rscript run_app.R` (or double-clicking
`run_app.command` on macOS); opens in the user's default browser at a local address — "nothing
leaves your machine." No accounts, no auth, no multi-tenant hosting.

Workflow: **Data & window** tab (upload a CSV or pick a bundled example, map x/y or lon/lat
columns, an optional group column for case/control, an optional weight column, an optional time
column for the spatiotemporal tabs, and the study-region boundary method) → **Build point
pattern** → one of four analysis tabs (**Density**, **Risk**, **Spatiotemporal density**,
**Spatiotemporal risk**) → configure parameters → an explicit **Run** (cross-validation bandwidth
selection and adaptive smoothing can take seconds to minutes, so nothing auto-recomputes on
parameter change) → explore results across **Plot / Map / Table / Download** sub-tabs. The two
spatiotemporal tabs add a time slider (sparr interpolates between computed grid points, not just
snapping to them), a Joint/Conditional toggle, and a Play/Pause animation that steps through the
computed time grid.

Requires R ≥ 4.1 installed locally. macOS additionally requires XQuartz (a non-obvious transitive
dependency: `sparr` imports `misc3d`, which imports `tcltk`, which needs XQuartz to load on
macOS at all, even though nothing here uses tcltk/X11 features directly).

## Capabilities and Constraints

**Implemented:**
- `bivariate.density()` — single-population density, with OS/NS/LSCV/LIK/fixed bandwidth
  selection, fixed or adaptive smoothing, edge correction, optional per-point weights.
- `risk()` — case/control relative risk with asymptotic tolerance contours, plus an explicit
  override for which group is the case/numerator (sparr itself just picks the alphabetically
  first factor level).
- `spattemp.density()` / `spattemp.risk()` — spatiotemporal versions; the risk variant treats the
  control group as static over time (time-varying case density vs. a static control), matching
  sparr's own documented usage pattern.
- Study region via bounding box, convex hull, or an uploaded boundary polygon — including
  multi-part polygons (e.g. mainland plus islands) with automatic ring reconstruction and
  winding-direction correction, plus a boundary buffer to absorb near-boundary points.
- Leaflet map overlay when the uploaded data is geographic (lon/lat).
- Downloads: PNG, PDF, the underlying surface as a CSV grid, and a generated reproducible R
  script that reproduces the exact analysis outside the app.

**Not yet built (explicit roadmap):** multi-scale density (`multiscale.density`), bootstrap
bandwidth selectors (`BOOT.density`, `BOOT.spattemp`), Monte Carlo tolerance contours and
point-level classification (`tol.classify`), a time-varying (rather than static) control group
for spatiotemporal risk, and standalone desktop packaging (currently requires a visible browser
window plus a local R process — there is no packaged, browser-free distribution).

**Terminology:**
- *Tolerance contours* — asymptotic p-value significance regions; the app supports contouring
  both areas of significantly elevated risk (solid line) and significantly reduced risk (dotted
  line, an addition beyond sparr's own default single-tail contouring).
- *Case / control* — the two levels of the group column in a risk analysis.
- *tgrid* — the actual time grid points a spatiotemporal analysis computed, as distinct from the
  finer, interpolated range the time slider exposes.

## Brand Commitments

None established beyond the existing app. Header/title text: "sparr — spatial relative risk
explorer". Single accent color (`#2C5F8A`, a blue) set via `bslib::bs_theme(version = 5)`;
otherwise default Bootstrap 5 / bslib styling, no custom theme, no logo. Voice, where it appears
(README, in-app help text), is direct and technically precise rather than promotional.

## Evidence on Hand

Three bundled example CSVs in `data/examples/` (`pbc.csv`, `burk.csv`, `fmd.csv`), exported
directly from `sparr`'s own built-in R datasets (see `data/examples/export_examples.R`) — real
data, but the package's own demonstration datasets, not a user's actual research data. No
testimonials, customer logos, case studies, pricing, or licensing claims exist and none should be
invented; this is a free local tool with no commercial framing.

## Product Principles

1. Faithful to `sparr`, not a reinterpretation — every control mirrors a real parameter `sparr`'s
   R functions expose; nothing here should silently produce different numbers than calling
   `sparr` directly with the same settings.
2. Statistical rigor over visual convenience — tolerance contours are actual computed
   significance regions, not implied by color choice alone.
3. No code required, but no black box either — the reproducible-script download lets a user
   verify or extend results by hand.
4. Explicit compute, not silent recompute — expensive operations always require a deliberate
   "Run," so a multi-minute wait is never a surprise.
5. Local and self-contained — no server, no account, no data leaving the user's machine.
