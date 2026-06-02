# Continuous-G Project Log

Timestamp: 2026-05-19, America/Chicago

This file is the consolidated detailed log for the continuous-G update to the
`fdid` package. It combines the historical project log, completion log,
estimator review verdicts, validation reports, tutorial reorganization notes,
and implementation notes.

Historical markdowns were archived under:

```text
historical_mds_2026-05-19/
```

The other active markdowns are:

- `TODO.md`;
- `CONTINUOUS_G_ESTIMATOR_REFERENCE.md`;
- `HANDOFF.md` (local restart note; excluded from package builds).

Package-standard `NEWS.md` and `LICENSE.md` remain separate.

## 2026-06-01 Tutorial Inference And DML Map Clarification

Scope: StatsClaw-guided tutorial update plus GitHub Pages branch sync.

- Expanded Chapter 3 with a concise reporting-layer map for continuous-G
  results: `summary()` scalar rows, `plot(type = "curve")`,
  `fdid_contrast()`, and `fdid_derivative()`. The chapter now states that
  kernel and `dml_flex` scalar event rows are interval-average level slopes over
  the fitted grid, not substitutes for the stored level/derivative curves.
- Updated Chapter 4 so the kernel tutorial matches the current robust
  implementation. The chapter now distinguishes analytical robust
  stacked-sandwich curve covariance from bootstrap replicate curves,
  percentile intervals, and finite-grid quantile-envelope bands. The reporting
  examples use the analytical robust fit for compact reporting and explicitly
  label the pointwise approximation used when the substantively useful
  full-support `G = 0` boundary contrast does not have usable cross-grid
  covariance entries.
- Expanded Chapter 5's DML tutorial with a target-specific inference table and
  a dedicated `signal_map` table. The text now separates `density_method`
  choices from second-stage curve mapping choices, and labels `blp_spline` as
  the current DML-flex route with stored `theta_vcov` / `delta_vcov`.
- Updated Chapter 6 to clarify that plotting is a display layer. Curve plots
  draw stored curve/inference objects, contrast plots call `fdid_contrast()`,
  and support panels are descriptive rather than inferential.
- Disabled the Quarto HTML source-code tool for the tutorial book. Normal code
  chunks still render, but hidden custom table-construction chunks no longer
  reappear through the page-level source view.
- Added an internal kernel bootstrap performance guard after the first full
  tutorial render attempt revealed that bootstrap replicate curves were also
  recomputing the analytical stacked-sandwich covariance. The main robust curve
  still stores `mu_vcov` / `delta_vcov`; bootstrap replicate fits now skip that
  redundant analytical covariance and compute `mu_boot_vcov` /
  `delta_boot_vcov` from the replicate matrices as intended. A briefly
  considered fixed-`h0` bootstrap shortcut was rejected because bootstrap
  should follow the theory-defined resampling estimator, not tutorial render
  time. The documented kernel bootstrap path continues to reselect the baseline
  bandwidth inside resamples.
- Planned ship scope: push the just-completed continuous-G inference/summary
  code changes and the updated tutorial sources/rendered docs to the existing
  GitHub development branch without adding local planning markdown files.

## 2026-06-01 Continuous-G Inference Clarification And Summary Update

Scope: local R-code and documentation update for continuous-G inference
semantics.

- Added analytical cross-grid covariance storage for the kernel estimator under
  `vartype = "robust"`. The event curve now stores `mu_vcov`, `delta_vcov`,
  and `curve_vcov_method`, computed from stacked local-linear sandwich scores.
  When a prepared cluster column is present, the sandwich meat is cluster-summed.
- Updated kernel scalar summaries and `fdid_contrast()` / `fdid_derivative()` so
  robust kernel fits can use covariance-aware inference instead of the previous
  endpoint-independent fallback.
- Updated `dml_flex` scalar event output from a grid-average level compatibility
  value to the interval-average level slope over the fitted `eval_g` grid.
  When `theta_vcov` is available, the scalar SE uses the covariance-aware
  endpoint formula divided by the grid span.
- Updated `summary.fdid()` to print the scalar event SE and CI method when
  stored, and to describe `dml_flex` scalar output as an interval-average
  theta slope rather than a grid-average theta value.
- Added a consolidated continuous-G inference chapter to
  `CONTINUOUS_G_ESTIMATOR_REFERENCE.md`, organized by output object:
  `summary()`, kernel curves, scalar DML, DML-flex curves, reporting helpers,
  replicate-vs-vcov logic, and remaining limitations.
- Updated FDID Notes to align the theory discussion with the implemented
  stacked-sandwich kernel covariance, covariance-aware reporting helpers, and
  interval-average summary semantics.
- Added tests for robust kernel covariance-aware contrasts and DML-flex
  interval-average summary semantics.

## 2026-06-01 Tutorial Summary And Curve Plot Pass

Scope: local tutorial cleanup plus GitHub Pages synchronization for the public
continuous-G tutorial branch.

- Simplified Chapter 2 plotting guidance so the binary/traditional FDID chapter
  now documents only `raw`, `dynamic`, `overlap`, and method-comparison plots.
  Detailed continuous-G `plot(type = "curve")`, `eval_g`, and contrast plotting
  guidance is left to Chapters 4-6.
- Added standard `summary()` calls to the Chapter 4 kernel workflow and the
  Chapter 5 DML workflow so `summary()` is the first reporting layer before
  curve plots, reporting helpers, and custom metadata tables.
- Added a point-estimate-first DML flexible curve figure in Chapter 5.5. The
  figure plots level and derivative curves without intervals and overlays the
  known synthetic DGP truth; Chapter 5.6 remains the dedicated curve-inference
  section.
- Kept custom tables only for quantities that `summary()` does not fully cover:
  stored grid excerpts, inference metadata, density-method comparisons, and
  targeted `fdid_contrast()` / `fdid_derivative()` reporting.

Verification:

- Re-rendered the tutorial site after the source edits.
- Re-synced rendered HTML, search index, and figure assets to the GitHub staging
  checkout before push.
- Compared local and staging tutorial sources, Chapter 2/4/5/6 HTML, and all
  referenced figure assets to avoid the earlier stale-image mismatch.

## 2026-06-01 GitHub Pages Chapter 4-6 Resync

Scope: local-to-GitHub staging synchronization fix after detecting that the
published GitHub Pages branch contained a mixed tutorial state.

- Compared local `fdid-main/docs` against the GitHub staging checkout
  `/private/tmp/enhan-fdid-sync` for Chapters 4-6.
- Found that Chapter 4 HTML had been updated, but several Chapter 4 PNGs in
  staging still came from an older trimmed-grid render.
- Found that Chapter 5 staging source/HTML still reflected an older mixed
  mortality-plus-synthetic workflow, while local `fdid-main` contained the
  current synthetic-only DML chapter.
- Found that Chapter 6 source was already consistent; HTML differences were
  generated-format whitespace rather than substantive content.
- Re-synced `tutorial/04-kernel-continuous-g.Rmd`,
  `tutorial/05-dml-continuous-g.Rmd`,
  `tutorial/06-visualization.Rmd`, `docs/04-kernel-continuous-g.html`,
  `docs/05-dml-continuous-g.html`, `docs/06-visualization.html`,
  `docs/search.json`, and all Chapter 4-6 rendered PNGs from local
  `fdid-main` to the staging checkout.
- Removed obsolete staging-only Chapter 5 PNGs
  `dml-flex-plots-1.png` and `dml-blp-band-plot-1.png`, because the current
  local Chapter 5 no longer references those figures.
- Stripped trailing whitespace in the generated Chapter 5 and Chapter 6 HTML
  in both local and staging copies to keep `git diff --check` clean.

Verification:

- File equality checks passed for Chapters 4-6 tutorial sources, HTML files,
  and `docs/search.json` between local `fdid-main` and staging.
- All PNG files referenced by Chapters 4-6 HTML matched byte-for-byte between
  local `fdid-main` and staging.
- The Chapter 4-6 figure directories matched between local and staging after
  excluding local macOS `.DS_Store` files.
- `git diff --check` passed in the staging checkout.

## 2026-06-01 Kernel Tutorial Reporting And Bootstrap Effective Draws

Scope: Chapter 4 tutorial/documentation update and GitHub dev-branch refresh.

- Added a Chapter 4.4 explanation of `effective_bootstrap_draws`: finite-grid
  uniform envelopes need a whole bootstrap curve to be finite at every stored
  `eval_g` point, so internally non-finite bootstrap curves can be excluded
  from the uniform-band calculation even when pointwise percentile intervals use
  columnwise finite draws.
- Clarified that this effective-draw filter is not an additional trimming or
  winsorization step; it reflects local-linear failures at one or more grid
  points in a resampled curve, which is most likely in sparse tails and for
  derivative curves.
- Changed Chapter 4.5 reporting examples to use the analytical
  `vartype = "robust"` kernel fit (`fit_kernel_fast`) rather than the
  bootstrap object (`fit_kernel_boot`).
- Set `inference = "pointwise"` for the Chapter 4.5
  `fdid_contrast()`/`fdid_derivative()` examples and contrast coefficient plot.
  The contrast is therefore labeled `pointwise_approx`, because the robust
  kernel object currently stores endpoint SEs but not a cross-grid covariance
  matrix for `V22 + V11 - 2 * V12`.
- Preserved the existing Chapter 4.4 bootstrap table and figures to avoid
  rerunning the expensive 200-draw kernel bootstrap. Only the analytical
  reporting table and contrast plot were refreshed.
- Updated `NEWS.md`, `tutorial/04-kernel-continuous-g.Rmd`,
  `docs/04-kernel-continuous-g.html`, `docs/search.json`, and the rendered
  Chapter 4 contrast plot.

Verification:

- Recomputed the Chapter 4.5 analytical robust reporting table and contrast
  plot from `fit_kernel_fast` without rerunning the bootstrap chunk.
- Confirmed the replacement contrast plot renders with
  `Inference: pointwise_approx`.

## 2026-06-01 Kernel Bootstrap 200 And Synthetic-Only DML Chapter

Scope: local-only tutorial/documentation update. No GitHub push or commit was
performed; the user will inspect the local build first. StatsClaw run
artifacts were recorded under
`/private/tmp/statsclaw-data/workspace/Statsclaw-FDID/runs/TUTORIAL-KERNEL200-DML-SYNTHETIC-20260601-0006`.

- Updated Chapter 4 so the kernel bootstrap example now renders with
  `boot = 200` rather than a quick low-draw bootstrap. The source uses
  `kernel_boot_draws <- 200L`.
- Updated Chapter 4 text to explain that 200 bootstrap replications make
  pointwise percentile intervals and finite-grid quantile-envelope bands less
  dominated by extreme resamples, while full-support right-tail and derivative
  uncertainty can still remain wide.
- Hid manual table-construction code in Chapter 4 rendered output for support,
  scalar summaries, curve excerpts, bootstrap metadata, reporting tables, and
  option checklists.
- Rewrote Chapter 5 as a synthetic-only DML continuous-G chapter. The chapter
  no longer loads or uses the mortality data.
- Rebuilt Chapter 5 around one DGP with four periods, reference period 2,
  event periods 3 and 4, `G ~ Beta(4, 4)`, and
  `tau(g) = 0.4 + 0.6g + 0.8g^2`.
- Used the same synthetic prepared data for scalar `dml_plr`, flexible
  `dml_flex`, curve inference, `fdid_contrast()`/`fdid_derivative()`,
  `dml_incremental`, learner availability, and high-dimensional covariate
  expansion sections.
- Hid manual table-construction code in Chapter 5 rendered output while
  keeping substantive `fdid()` and related package calls visible.
- Updated `NEWS.md`, rendered Chapter 4 and Chapter 5 docs, `docs/search.json`,
  and the StatsClaw run artifacts for this local workflow.

Verification:

- Rendered `tutorial/04-kernel-continuous-g.Rmd` with `boot = 200`: pass. The
  kernel bootstrap chunk took about 25 minutes locally.
- Rendered `tutorial/05-dml-continuous-g.Rmd`: pass.
- Text audit: Chapter 5 source and rendered HTML contain no mortality-data
  workflow.
- Rendered Chapter 4 confirms `kernel_boot_draws <- 200L`; the reported
  uniform-band effective draw count is `107`, reflecting the package's usable
  finite-grid bootstrap draws after internal filtering.
- Visual checks: Chapter 4 bootstrap figures render with boxed legends below
  the plots; Chapter 5 synthetic DML level and derivative curves are now
  readable on the same synthetic support, with DGP truth overlaid.
- Hidden-code audit: manual explanatory table chunks in Chapters 4 and 5 use
  `echo = FALSE`; the rendered article keeps the main `fdid()` calls visible.
- `R CMD INSTALL .`: pass in local `fdid-main`.
- `Rscript -e 'devtools::test()'`: pass with `375` passed, `0` failed, and
  one existing warning about `target.pop` not changing `method = "did"`.
- `git diff --check`: pass for the affected tutorial/docs/news/log files.

## 2026-06-01 Kernel Full-Support Grid And DML Tutorial Consolidation

Scope: local-only tutorial/documentation update after the previous local
default-grid update. No GitHub push or commit was performed; the user will
inspect the local build first. StatsClaw run artifacts were recorded under
`/private/tmp/statsclaw-data/workspace/Statsclaw-FDID/runs/TUTORIAL-TRIM-DML-STRUCTURE-20260601-0005`.

- Revised Chapter 4 to state both the package default and the chapter-specific
  choice: `eval_g = NULL` uses a 50-point grid controlled by `trim`, and the
  package default is `trim = 0.05`, but the mortality kernel chapter now sets
  `trim = 0` so the automatic 50-point grid spans the full observed `G`
  support.
- Added support diagnostics showing the package-default 95% upper endpoint,
  the chapter's full-support upper endpoint, the number/share of observations
  above the package-default 95% endpoint, and the stored 50-point grid range.
- Simplified the Chapter 4 support plot so it no longer draws all 50
  evaluation-grid points as vertical lines. The plot now marks only the
  full-support upper endpoint, the package-default 95% endpoint, and a few
  positive-support landmarks.
- Updated both the robust and bootstrap kernel tutorial fits to keep
  `eval_g = NULL` but pass `trim = 0`.
- Consolidated Chapter 5's synthetic DML presentation. The synthetic DGP,
  flexible DML fits, curve excerpt, inference metadata, CI curve figure, and
  reporting helper now form one continuous workflow rather than separate
  repeated synthetic displays.
- Removed the redundant no-CI synthetic flexible-curve figure and redundant
  single-panel BLP level-band figure from Chapter 5.
- Restored `Learner Choice` and `High-Dimensional Expansions` as top-level
  Chapter 5 sections rather than subsections under incremental DML.
- Updated `NEWS.md`, rendered Chapter 4 and Chapter 5 docs, and the StatsClaw
  run artifacts for this local workflow.

Verification:

- Rendered `tutorial/04-kernel-continuous-g.Rmd`: pass.
- Rendered `tutorial/05-dml-continuous-g.Rmd`: pass.
- Rendered Chapter 4 text confirms `trim = 0.05` is the package default and
  the chapter uses `trim = 0`; the stored kernel grid reports 50 points with
  range `[0, 2.148173]`.
- Rendered Chapter 4 support figure no longer draws the 50 evaluation points
  as vertical lines.
- Rendered Chapter 5 text confirms `Learner Choice` is section 5.8 and
  `High-Dimensional Expansions` is section 5.9.
- Text audit: no remaining `dml-flex-plots`, `dml-blp-band-plot`, or
  `Synthetic DML CI Illustration` section in Chapter 5.
- `R CMD INSTALL .`: pass in local `fdid-main`.
- `Rscript -e 'devtools::test()'`: pass with `375` passed, `0` failed, and
  one existing warning about `target.pop` not changing `method = "did"`.
- `git diff --check`: pass for the affected tutorial/docs/news/log files.

## 2026-06-01 Default-Grid Tutorial And Outside-Legend Update

Scope: local-only update after pushed commit `d261353`, targeting
continuous-G tutorial chapters, continuous-G plot legend layout, rendered
docs, and package documentation. No GitHub push or commit was performed for
this entry; the user will inspect the local build first. StatsClaw run
artifacts were recorded under
`/private/tmp/statsclaw-data/workspace/Statsclaw-FDID/runs/TUTORIAL-EVALG-LEGEND-20260601-0004`.

- Revised Chapter 4 so the kernel tutorial no longer demonstrates a sparse
  six-point estimation grid. The mortality kernel examples now omit `eval_g`
  and therefore use the package default 50-point trimmed grid.
- Added explicit Chapter 4 text distinguishing the estimation grid from a
  plotting selection: `eval_g = NULL` evaluates
  `seq(quantile(G, trim), quantile(G, 1 - trim), length.out = 50)`, while
  user-supplied `eval_g` estimates exactly the supplied points.
- Clarified that `trim` controls default grid endpoints only. It does not drop
  observations or winsorize `G`, so a mass point at `G = 0` remains available
  to the kernel estimator and the default lower endpoint can be zero when the
  lower trimmed quantile is zero.
- Revised Chapter 5 so flexible DML curve examples also omit `eval_g` and use
  the same default 50-point grid.
- Reframed Chapter 5 around a complete synthetic DML curve workflow, including
  the data-generating process and the true level/derivative curves, because
  the mortality DML curve examples are slower and less visually stable.
- Kept mortality data for scalar DML and incremental comparisons, while using
  the synthetic design for full DML curve plotting and inference displays.
- Updated `plot.fdid()` so `legend.position = "auto"` draws boxed legends below
  the plot region for raw, overlap, and continuous-G curve plots. Explicit
  corner placements remain available for intentional inside-plot legends, and
  `"none"` still suppresses legends.
- Added a dedicated below-plot legend panel when continuous-G curves are drawn
  with `support.panel = "separate"`.
- Updated Chapter 6 to document the default 50-point continuous-G grid,
  `trim`, mass-point behavior, and the revised boxed-below-plot legend
  behavior.
- Updated `NEWS.md`, `man/plot.fdid.Rd`, rendered `docs/`, and the StatsClaw
  run artifacts for this local workflow.

Verification:

- `R CMD INSTALL .`: pass in the staging checkout.
- Rendered `tutorial/04-kernel-continuous-g.Rmd`: pass.
- Rendered `tutorial/05-dml-continuous-g.Rmd`: pass.
- Rendered `tutorial/06-visualization.Rmd`: pass.
- Re-rendered the same Chapters 4, 5, and 6 inside local `fdid-main/tutorial`
  after syncing the staging outputs back: pass.
- `R CMD INSTALL .`: pass in local `fdid-main`.
- Kernel default-grid check: `eval_g = NULL` stores 50 evaluation points with
  lower endpoint `0` in the mortality tutorial data.
- DML flexible-curve default-grid check: `eval_g = NULL` stores 50 evaluation
  points in the synthetic DML tutorial workflow.
- Visual checks: curve legends render in boxed areas below the plotting region;
  embedded support histograms and x-axis labels are no longer covered by
  legends; the kernel fixed-reference contrast plot keeps readable target
  labels.
- `Rscript -e 'devtools::test()'`: pass with `375` passed, `0` failed, and
  one existing warning about `target.pop` not changing `method = "did"`.
- `git diff --check`: pass after trimming generated Quarto HTML trailing
  whitespace.
- Synced updated `tutorial/`, `docs/`, `R/plot.R`, `man/plot.fdid.Rd`, and
  `NEWS.md` from the staging checkout back into local `fdid-main`.

## 2026-06-01 Continuous-G Plot Layout Polish

Scope: current GitHub dev branch
`continuous-g-dev-2026-06-01-tutorial-enrichment`, plotting layer,
continuous-G tutorial figures, rendered docs, and local/GitHub sync. Pushed
commit: `d261353` (`Polish continuous-G plot layouts`).

- Shortened visible curve legend labels while preserving detailed inference
  descriptions in reporting helpers and metadata.
- Added `legend.position` to `plot.fdid()` with `"auto"` placement, standard
  base-graphics corner positions, and `"none"` for suppressing legends.
- Suppressed redundant one-item curve legends in auto mode and discouraged
  automatic legend placement over embedded support displays.
- Fixed contrast plot x-axis label/tick overlap by drawing ticks and labels
  separately and expanding the bottom margin when labels are multiline.
- Updated Chapter 4 so the kernel contrast coefficient plot uses `G = 0` as
  the fixed reference and absolute level targets rather than Q1-Q5 labels.
- Updated Chapter 5 so mortality DML level and derivative plots both suppress
  wide CI layers and show point-estimate diagnostics only.
- Removed the Chapter 5 incremental point plot while retaining the incremental
  comparison table.
- Reorganized Chapter 5 so `Learner Choice` and `High-Dimensional Expansions`
  are subsections under section 5.8.
- Updated Chapter 6 to document `legend.position`.
- Updated `NEWS.md`, `man/plot.fdid.Rd`, tests, rendered `docs/`, and
  StatsClaw run artifacts under
  `/private/tmp/statsclaw-data/workspace/Statsclaw-FDID/runs/PLOT-POLISH-TUTORIAL-20260601-0003`.

Verification:

- `R CMD INSTALL .`: pass.
- Rendered `tutorial/04-kernel-continuous-g.Rmd`: pass.
- Rendered `tutorial/05-dml-continuous-g.Rmd`: pass.
- Rendered `tutorial/06-visualization.Rmd`: pass.
- Visual checks: kernel contrast labels no longer collide with tick marks;
  DML mortality level/derivative plots no longer show wide CI layers; BLP and
  synthetic CI legends are compact and no longer use long method strings.
- Text audit: no remaining `dml-incremental-point-plot` output and no old long
  legend labels in R/man/tutorial/docs/tests.
- `Rscript -e 'devtools::test()'`: pass with `375` passed, `0` failed, and
  one existing warning about `target.pop` not changing `method = "did"`.
- `git diff --check`: pass.
- GitHub push: pass to `Enhan-Liu/fdid` branch
  `continuous-g-dev-2026-06-01-tutorial-enrichment`.

## 2026-06-01 Kernel/DML Tutorial Figure Refinement

Scope: current GitHub dev branch
`continuous-g-dev-2026-06-01-tutorial-enrichment`, continuous-G tutorial
figures, rendered docs, and local/GitHub sync. Pushed commit: `f7edcb8`
(`Refine continuous-G tutorial figures`).

- Revised Chapter 4 so the kernel bootstrap, pointwise/uniform band,
  reporting, and contrast examples all keep using the mortality data rather
  than switching to a synthetic kernel data-generating process.
- Added a fixed-reference kernel contrast coefficient plot with the reference
  point shown as the algebraic zero anchor.
- Revised Chapter 5 so the noisy mortality DML derivative and incremental
  displays emphasize point estimates without confidence intervals.
- Added a separate synthetic DML CI illustration for level and derivative
  curve uncertainty, preserving the mortality application as the main running
  example.
- Replaced plot-image titles that showed raw `Xdistr = ...` argument values
  with substantive figure titles in Chapters 4, 5, and 6.
- Updated `NEWS.md`, rendered `docs/`, and the StatsClaw run artifacts under
  `/private/tmp/statsclaw-data/workspace/Statsclaw-FDID/runs/TUTORIAL-KERNEL-DML-PLOT-20260601-0002`.

Verification:

- Rendered `tutorial/04-kernel-continuous-g.Rmd`: pass.
- Rendered `tutorial/05-dml-continuous-g.Rmd`: pass.
- Rendered `tutorial/06-visualization.Rmd`: pass.
- Visual checks: kernel bootstrap, kernel contrast, DML mortality curve,
  synthetic DML CI, DML incremental point, and Chapter 6 support-display
  figures rendered with the intended titles and layout.
- Text audit: no remaining Chapter 4 kernel synthetic demo marker and no
  `main = "Xdistr ..."` plot titles in tutorial/docs output.
- `Rscript -e 'devtools::test()'`: pass with `374` passed, `0` failed, and
  one existing warning about `target.pop` not changing `method = "did"`.
- `git diff --check`: pass.
- GitHub push: pass to `Enhan-Liu/fdid` branch
  `continuous-g-dev-2026-06-01-tutorial-enrichment`.

## 2026-06-01 GRF Tutorial And Curve-Plot Refinement

Scope: current GitHub dev branch
`continuous-g-dev-2026-06-01-tutorial-enrichment`, tutorial rendering,
continuous-G curve plotting, and local/GitHub sync. Pushed commit:
`3aaa34d` (`Update GRF tutorial examples and curve plots`).

- Updated DML tutorial examples from `learner = "linear"` to
  `learner = "grf"` for binary DML, continuous PLR, flexible DML curves,
  incremental DML, and DML visualization examples.
- Added tutorial prose explaining that GRF examples use generalized random
  forests for cross-fitted nuisance learning, while the package default remains
  `learner = "linear"` for lightweight runs.
- Added a Chapter 5 note that GRF examples take noticeably longer to render
  because they fit multiple forests across folds and targets.
- Fixed a cache-consistency issue: Chapter 2 DML binary chunks and Chapter 6
  default chunks no longer reuse old cached output, preventing rendered
  summaries from disagreeing with the code.
- Changed `plot.fdid(type = "curve", curve = "both")` from vertically stacked
  panels to side-by-side level/derivative panels with embedded support
  displays.
- Updated `plot.fdid` help text and `NEWS.md` for the side-by-side
  `curve = "both"` behavior.
- Standardized Chapter 4 and Chapter 5 estimator-chapter curve plots to use
  `Xdistr = "histogram"`, `support.panel = "embedded"`, and
  `main = "Xdistr = histogram"`.
- Left the broader support-display, marker, contrast, and interval option
  demonstrations in Chapter 6.
- Removed obsolete rendered Chapter 4 support-display option images after the
  source chunk was removed.
- Synced the pushed `tutorial/`, `docs/`, `R/plot.R`, `man/plot.fdid.Rd`, and
  `NEWS.md` back into local `fdid-main`.

Verification:

- `requireNamespace("grf", quietly = TRUE)`: pass.
- `R CMD INSTALL .`: pass before rendering.
- `Rscript -e 'devtools::test()'`: pass with `374` passed, `0` failed, and
  one existing warning about `target.pop` not changing `method = "did"`.
- Full `/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto
  render tutorial`: pass; GRF rendering was slower as expected.
- Targeted re-renders for `tutorial/02-main.Rmd`,
  `tutorial/05-dml-continuous-g.Rmd`, and `tutorial/06-visualization.Rmd`:
  pass.
- Rendered-output cache audit: pass; Chapter 6 now reports `Learner: grf` in
  the displayed DML incremental summary.
- Visual checks: `curve = "both"` is side-by-side; Chapter 4 and Chapter 5
  curve figures use histogram embedded support.
- `git diff --check`: pass after trimming generated Quarto trailing
  whitespace.
- GitHub push: pass to `Enhan-Liu/fdid` branch
  `continuous-g-dev-2026-06-01-tutorial-enrichment`.

## 2026-05-31 Target-First Tutorial Cleanup

Scope: tutorial source consolidation, rendered GitHub Pages refresh, project
planning/log updates, and architecture/run artifacts. No R estimator,
inference, plotting, `man/`, `tests/`, `NAMESPACE`, or `DESCRIPTION` files
were edited.

- Consolidated the public tutorial source into one canonical `tutorial/`
  directory.
- Rebuilt the tutorial as eight target-first chapters:
  `01-start.Rmd`, `02-main.Rmd`, `03-continuous-overview.Rmd`,
  `04-kernel-continuous-g.Rmd`, `05-dml-continuous-g.Rmd`,
  `06-visualization.Rmd`, `07-reporting-inference.Rmd`, and
  `08-sensitivity.Rmd`.
- Removed obsolete tutorial source copies:
  `tutorial-full-continuous-g/` and `tutorial-next-2026-05-30/`.
- Removed stale placeholder `tutorial/summary.qmd`.
- Updated `tutorial/_quarto.yml` to render directly into `docs/` and removed
  references to deleted chapter names.
- Regenerated `docs/` from the final tutorial source and restored
  `docs/.nojekyll`.
- Centralized conservative inference language: kernel bootstrap bands are
  quantile-envelope bands, non-BLP `dml_flex` bands are practical
  signal-mapping multiplier bands, and BLP-spline DML is a signal map with
  finite-grid covariance/Gaussian max-`t` bands when available.
- Updated `TUTORIAL_CLEANUP_PLAN_2026-05-31.md`, `TODO.md`,
  `PACKAGE_PLOTTING_LAYER_UPDATE_PLAN_2026-05-31.md`, and architecture/run
  artifacts.

Verification:

- `quarto render tutorial`: blocked because `quarto` was not on `PATH`.
- `/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto render tutorial`:
  pass; output created in `docs/` with the final chapter names.
- Stale rendered chapter scan: pass; no deleted chapter pages remain in
  `docs/`.
- `git diff --cached --check`: pass after trimming generated HTML trailing
  whitespace.
- File hygiene scan for `.DS_Store`, `.Rhistory`, and `Rplots.pdf`: pass.
- Staged scope scan: pass; no `R/`, `man/`, `tests/`, `NAMESPACE`, or
  `DESCRIPTION` files changed.
- Tester dispatch was blocked by the child Codex usage limit, so the leader ran
  the same observable audit checklist as a fallback.
- R source parse and `tools::checkRd`: pass; locale startup warnings only.
- `R CMD INSTALL -l /private/tmp/fdidlib .`: pass.
- RStudio-bundled Quarto render with the installed package library: pass;
  validation regenerated whitespace-only HTML differences, which were restored
  to the committed docs after confirming the render succeeds.
- `R CMD build --no-build-vignettes --no-manual /private/tmp/enhan-fdid-sync`:
  pass; built `/private/tmp/fdid_1.0.2.tar.gz`.
- `Rscript -e 'devtools::test(".")'`: pass with `374` passed, `0` failed, and
  one existing warning about `target.pop` not changing `method = "did"`.

## 2026-05-31 Review-Driven Package-Impact Fixes

Scope: package operation and build-impact items from the StatsClaw
`REVIEW-2026-05-31-fdid-main-package` review. Process-only review blockers
such as full run-directory architecture artifacts and full tester evidence are
tracked separately from these package edits.

- Renamed the local Gaussian max-t statistic in `.dml_gaussian_max_t_crit()`
  from `T` to `max_t_stat`, satisfying the active R-package profile rule
  against using `T`/`F` as variable names.
- Added `^HANDOFF\\.md$` to `.Rbuildignore` so the local restart/provenance
  handoff note is not included in package builds.
- Updated the new tutorial copy's local-install example to use a generic
  `fdid_<version>.tar.gz` placeholder instead of stale `fdid_0.1.0.tar.gz`.
- Replaced the tutorial dependency helper's `installed.packages()` scan with
  `requireNamespace()` checks, avoiding a slow package-wide installed-package
  query in rendered tutorial code examples.

Verification:

- `Rscript` parse/checkRd over `R/` and `man/`: pass.
- `Rscript` purl/parse for `tutorial-next-2026-05-30/01-start.Rmd`: pass.
- `R CMD build --no-build-vignettes --no-manual fdid-main` from `/private/tmp`:
  pass; the built tarball excludes `HANDOFF.md`, `fdid.Rcheck`,
  tutorial-next outputs, and `tests/testthat/Rplots.pdf`.
- `R CMD INSTALL -l /private/tmp/fdidlib fdid-main`: pass.
- `R_LIBS=/private/tmp/fdidlib Rscript -e 'devtools::test("fdid-main")'`:
  374 passes, 0 failures, 1 existing target-pop warning, 0 skips.
- `R_LIBS=/private/tmp/fdidlib Rscript -e 'devtools::run_examples("fdid-main")'`:
  pass; roxygen outputs were refreshed mechanically by devtools.
- `R_LIBS=/private/tmp/fdidlib quarto render fdid-main/tutorial-next-2026-05-30`:
  pass; output created at `tutorial-next-2026-05-30/_book/index.html`.

## 2026-05-31 Contrast Plot Refinement Update

Scope: `R/plot.R`, `man/plot.fdid.Rd`, continuous-G tests, the estimator
reference, tutorial-next visualization chapter, project log, TODO tracking,
and the package plotting-layer plan. Package-wide ggplot/autoplot work remains
out of scope.

- Added `include.ref = TRUE` as the default for `plot(type = "contrast")`, so
  fixed-reference contrast plots now include the reference value as an
  algebraic zero row.
- The reference row is constructed directly with estimate, SE, and confidence
  limits equal to zero, rather than calling `fdid_contrast(ref, ref)`. This
  avoids spurious nonzero pointwise fallback intervals for a contrast that is
  exactly zero by construction.
- Added character reference rules for `ref.g`: `"median"`, `"mean"`, `"min"`,
  `"max"`, and `"first"`, while preserving numeric `ref.g`.
- Added quantile target selection through `target.type = "quantile"` and
  `quantile.source = "support"` or `"eval_g"`. For example,
  `target.g = 5`, `target.g = "5"`, or `target.g = "q5"` draws five
  quantile-selected rows and uses the first selected quantile as the default
  reference.
- Kept non-reference inference delegated to `fdid_contrast()`. The plotting
  layer still does not compute contrast SEs or CIs for estimated contrasts.
- Updated tests, help, estimator reference, tutorial copy, TODO, and the
  package plotting-layer plan to describe the new default reference row and
  quantile target mode.

Verification:

- `Rscript` parse/checkRd over `R/` and `man/`: pass.
- `R CMD INSTALL -l /private/tmp/fdidlib fdid-main`: pass.
- `R_LIBS=/private/tmp/fdidlib Rscript -e 'setwd("fdid-main"); testthat::test_file("tests/testthat/test-continuous.R")'`:
  297 passes, 0 failures, 0 warnings, 3 CRAN skips.
- `R_LIBS=/private/tmp/fdidlib Rscript -e 'devtools::test("fdid-main")'`:
  374 passes, 0 failures, 1 existing target-pop warning, 0 skips.
- `R_LIBS=/private/tmp/fdidlib quarto render fdid-main/tutorial-next-2026-05-30`:
  pass; output created at `tutorial-next-2026-05-30/_book/index.html`.
- StatsClaw sidecar review: planner pass completed with non-blocking
  implementation notes; final reviewer returned PASS WITH NOTE, and the noted
  plan-document consistency issues were addressed.

## 2026-05-31 Fixed-Reference Contrast Plot Update

Status: superseded later on 2026-05-31 by the contrast plot refinement above,
which adds the default algebraic-zero reference row, character `ref.g` rules,
and quantile target selection.

Scope: `R/plot.R`, `man/plot.fdid.Rd`, continuous-G tests, the estimator
reference, tutorial-next visualization chapter, project log, TODO tracking, and
the package plotting-layer plan. Package-wide ggplot/autoplot work remains out
of scope.

- Implemented `plot.fdid(type = "contrast")` for `method = "kernel"` and
  `method = "dml_flex"`.
- Added internal `.fdid_contrast_plot_data()`, which selects a reference
  level and target `G` values, then calls `fdid_contrast()` for every displayed
  coefficient. The plot layer does not duplicate contrast SE or CI logic.
- Added internal `.fdid_plot_contrast_base()`, a base-graphics coefficient plot
  with point estimates, optional 95 percent CI bars when available, a zero
  reference line, and an inference label. The existing `ci = FALSE` plot
  option suppresses contrast CI bars just as it suppresses interval displays
  elsewhere.
- Added public plot arguments for the new contrast path: `ref.g`, `target.g`,
  `n.contrast`, `max.contrast`, `inference`, and `interpolate`. The existing
  `contrast.values` argument is also accepted as target values for
  `type = "contrast"`, but cannot be combined with `target.g`.
- The initial default target selection used stored finite `eval_g` values and
  omitted the reference value. This was later refined so the reference row is
  shown by default and the display cap counts all displayed rows.
- Contrast plots are level-contrast displays only:
  `mu(g_j) - mu(g_ref)` for kernel and
  `theta(g_j) - theta(g_ref)` for `dml_flex`. They are not new estimators and
  do not re-estimate or densify the curve.
- Inference semantics are inherited from `fdid_contrast()`: stored covariance
  first, then stored replicate curves, then stored bands, then labeled
  pointwise fallback.
- Updated the latest tutorial visualization chapter to replace the previous
  manual contrast-curve example with the formal `plot(type = "contrast")`
  coefficient plot.

Verification:

- `Rscript -e 'invisible(parse("fdid-main/R/plot.R")); invisible(tools::checkRd("fdid-main/man/plot.fdid.Rd"))'`: pass.
- `R CMD INSTALL -l /private/tmp/fdidlib fdid-main`: pass.
- `R_LIBS=/private/tmp/fdidlib Rscript -e 'setwd("fdid-main"); testthat::test_file("tests/testthat/test-continuous.R")'`:
  297 passes, 0 failures, 0 warnings, 3 CRAN skips after the refinement.
- `R_LIBS=/private/tmp/fdidlib Rscript -e 'devtools::test("fdid-main")'`:
  349 passes, 0 failures, 1 existing target-pop warning, 0 skips.
- `R_LIBS=/private/tmp/fdidlib quarto render fdid-main/tutorial-next-2026-05-30`:
  pass; output created at `tutorial-next-2026-05-30/_book/index.html`.

## 2026-05-31 StatsClaw Package Audit And Tutorial Copy Update

Scope: package-wide continuous-G consistency audit, tutorial copy and two new
continuous-G tutorial chapters, `.Rbuildignore`, `fdid()` argument docs, and
TODO tracking. The original `tutorial/` directory was left in place; the new
work was done in `tutorial-next-2026-05-30/`.

- Ran a StatsClaw-style package audit with two read-only sidecar reviews:
  one for theory/implementation/documentation consistency across the package,
  and one for the `interflex` tutorial structure and reusable teaching
  patterns.
- No blocking theory-to-implementation conflicts were found. The audit
  reaffirmed the current inference labels:
  kernel bootstrap quantile-envelope bands are not studentized sup-\(t\)
  bands; scalar DML uses analytical score or scalar score-multiplier
  inference, not a refit bootstrap; non-BLP `dml_flex` bands remain practical
  signal-residual multiplier bands; `signal_map = "blp_spline"` stores
  finite-grid covariance matrices and uses Gaussian max-\(t\) bands when
  usable.
- Fixed a documentation inconsistency in the `boot` argument description:
  kernel bootstrap replications control both level and derivative bands and
  stored replicate curves, not only derivative bands.
- Clarified the historical DML TODO/log language so it no longer reads as if
  `dml_flex` has no covariance storage. The remaining open issue is
  theorem-level efficient covariance or influence-function storage for
  non-BLP maps and all nuisance stages; the BLP-spline path already stores
  finite-grid level and derivative covariance matrices.
- Added `tutorial-next-2026-05-30/03a-kernel-continuous-g.Rmd`. The chapter
  follows the `interflex` applied pattern: support diagnostics, estimation
  grid, basic kernel fit, level and derivative plots, pointwise versus
  bootstrap band display, reporting helpers, and an option checklist.
- Added `tutorial-next-2026-05-30/03b-dml-continuous-g.Rmd`. The chapter
  follows the `interflex` DML teaching pattern: target map, nuisance learning
  and orthogonal signal overview, scalar PLR, flexible DML curves,
  DML-specific inference labels, BLP-spline covariance bands, reporting
  helpers, incremental average derivative, and learner/design sensitivity.
- Expanded `tutorial-next-2026-05-30/04-visualization.Rmd` so Chapter 6 now
  documents the current continuous-G curve plotting controls: `curve`,
  `interval`, `show.uniform.CI`, `Xdistr`, `support.panel`, `show.eval_g`,
  marker values, support display styles, and base-graphics styling arguments.
  The chapter also explains why kernel and DML plots can look different when
  their support-display options differ.
- Replaced the earlier manual fixed-reference contrast example,
  `theta(g) - theta(g_ref)`, with the first-class
  `plot(type = "contrast")` coefficient plot, which still delegates inference
  to `fdid_contrast()`.
- Expanded the DML tutorial's conditional-density discussion. It now separates
  the three `density_method` choices (`residual_kde`, `location_scale`,
  `local_kde`) from `signal_map`, and explicitly clarifies that
  `signal_map = "blp_spline"` is a B-spline best linear projection mapping and
  covariance path rather than another density estimator.
- Added `PACKAGE_PLOTTING_LAYER_UPDATE_PLAN_2026-05-31.md` as the updated
  plotting-layer sequence. It makes the fixed-reference continuous-G contrast
  coefficient plot the next code pass, before package-wide plot-data /
  ggplot/autoplot modernization.
- Updated the new tutorial copy's `_quarto.yml`, `index.qmd`, and
  `references.bib` so the kernel and DML chapters are part of the rendered
  book and cite the CME guide, Hainmueller-Mummolo-Xu, Chernozhukov et al.,
  Semenova-Chernozhukov, Kennedy et al., and Robinson.
- Added both alternate tutorial copies to `.Rbuildignore`, matching the
  existing policy that the tutorial book is source-facing documentation and is
  not included in the R package build.

Verification:

- `Rscript -e 'invisible(lapply(list.files("fdid-main/R", pattern = "[.]R$", full.names = TRUE), parse))'`: pass.
- `Rscript -e 'rd <- list.files("fdid-main/man", pattern = "[.]Rd$", full.names = TRUE); invisible(lapply(rd, tools::checkRd))'`: pass.
- `git diff --check -- fdid-main/R/fdid.R fdid-main/man/fdid.Rd fdid-main/CONTINUOUS_G_PROJECT_LOG.md fdid-main/TODO.md fdid-main/tutorial-next-2026-05-30`: pass.
- Extracted and parsed the R code from both new tutorial chapters and the
  updated visualization chapter with `knitr::purl()`: pass.
- `R CMD INSTALL fdid-main`: pass.
- `quarto render fdid-main/tutorial-next-2026-05-30`: pass; output created at
  `tutorial-next-2026-05-30/_book/index.html`.
- `devtools::test("fdid-main", filter = "continuous")`: 256 passes, 0
  failures, 0 warnings.
- `devtools::test("fdid-main")`: 327 passes, 0 failures, 1 known warning about
  `target.pop` not changing `method = "did"`.
- `LC_ALL=C LANG=C LC_CTYPE=C R CMD check --no-manual --no-build-vignettes fdid_1.0.2.tar.gz`:
  pass with 2 warnings caused by `--no-build-vignettes` leaving `inst/doc`
  absent.

## 2026-05-30 DML Inference Update

Scope: `R/dml.R`, `R/fdid.R`, `R/reporting.R`, continuous-G tests, FDID Notes,
the estimator reference, and TODO tracking. Tutorial files, `NEWS.md`, and
`HANDOFF.md` are intentionally out of scope for this pass.

- Added DML object metadata so fitted objects record the target estimand,
  scale, population, inference method, inference scope, and orthogonal-score
  label.
- Scalar DML methods now store event-period scalar inference details in
  `scalar_event`, including influence values, split estimates, multiplier
  draws, multiplier critical values, and the scalar inference method actually
  used.
- Added public scalar-DML inference controls:
  `dml_inference = "score"` or `"score_multiplier"`, `dml_boot`, and
  `dml_multiplier = "normal"`, `"wild"`, or `"exponential"`.
- `dml_binary`, `dml_plr`, and `dml_incremental` keep analytical score
  intervals by default and switch scalar CIs to score-multiplier critical
  values when `dml_inference = "score_multiplier"` is requested and the
  multiplier process is usable.
- Added `signal_map = "blp_spline"` for `dml_flex`: a B-spline best linear
  projection of the cross-fitted signal on `G` with stored level and derivative
  grid covariance matrices.
- `dml_flex` BLP-spline fits now store `theta_vcov`, `delta_vcov`, and
  `level_derivative_vcov`; Gaussian max-t critical values from the stored
  covariance matrices are used for finite-grid level and derivative bands when
  available.
- `fdid_contrast()` and `fdid_derivative()` now use stored curve covariance
  matrices before replicate curves in `inference = "auto"`, and support
  covariance-aware linear interpolation, not only exact-grid covariance.
- Removed the previous internal cap on `dml_flex` multiplier replications:
  `boot` now controls the practical residual-multiplier draws and Gaussian
  max-t simulation draws.
- Added continuous-G tests for scalar DML multiplier inference, DML metadata,
  BLP-spline covariance storage, covariance-aware contrast SEs using
  `V11 + V22 - 2 V12`, covariance-aware derivative SEs, and new validation
  errors.
- Verification: `devtools::test("fdid-main", filter = "continuous")` produced
  233 passes, 0 warnings, and 0 failures.

## 2026-05-30 DML Plot Alignment Update

Scope: `R/plot.R`, `man/plot.fdid.Rd`, continuous-G plot tests, the estimator
reference, TODO tracking, and the DML plot plan. Package-wide ggplot/autoplot
work, tutorials, `NEWS.md`, and `HANDOFF.md` remain out of scope.

- Executed the DML plot plan's Phase A rather than the package-wide plotting
  plan: the update keeps base graphics as the active plotting path and focuses
  on aligning `dml_flex` curve displays with the current kernel curve controls.
- Extended the internal curve plot-data helper with DML-aware metadata fields:
  target estimand, inference scope, inference method, `signal_map`,
  `density_method`, pointwise label, band label, band caveat, practical-band
  method, Gaussian critical value, simultaneous critical value, and curve
  covariance availability.
- Base `dml_flex` legends now distinguish practical signal-residual multiplier
  bands from BLP-spline Gaussian max-t covariance bands.
- The plot-data caveats explicitly preserve the inference boundary: practical
  `dml_flex` bands are not final theorem-level DML uniform bands, and BLP
  Gaussian max-t bands are not DML refit bootstrap bands.
- Added focused continuous-G tests for DML curve plot controls, DML-specific
  plot-data labels, BLP-spline band labels, `show.uniform.CI` suppression, and
  rejection of scalar DML methods as curve plots.

## 2026-05-30 Kernel Plot Update

Scope: `R/plot.R`, `man/plot.fdid.Rd`, continuous-G plot tests, the estimator
reference, and TODO tracking. Tutorial files, `NEWS.md`, and `HANDOFF.md` are
intentionally out of scope for this pass.

- Added an internal continuous-G curve plotting data layer so base graphics
  rendering resolves pointwise intervals, simultaneous bands, metadata, and
  `xlim` filtering in one place.
- Added a contiguous-segment helper so pointwise ribbons and band lines do not
  connect across missing interval values.
- Added `interval = "auto"`, `"pointwise"`, `"uniform"`, `"both"`, and
  `"none"` for continuous-G curve plots.
- Added `show.uniform.CI` to mirror the interflex-style control for displaying
  stored simultaneous bands.
- Added `curve = "both"` for two-panel level and derivative displays.
- Added `Xdistr = "rug"`, `support.panel = "embedded"` or `"separate"`,
  `show.eval_g`, optional `diff.values`/`contrast.values` markers, and base
  style controls for curve lines, ribbons, bands, and support displays.
- Kernel bootstrap plots now distinguish bootstrap percentile pointwise
  intervals from finite-grid bootstrap quantile-envelope simultaneous bands.
- Kernel robust pointwise normal intervals are no longer displayed or labeled
  as uniform bands.
- Added continuous-G plot tests for interval controls, support displays,
  helper output, robust-versus-bootstrap band behavior, NA-gap segmentation,
  and graphics-state restoration after `curve = "both"`.

## 2026-05-30 Kernel Inference Update

Scope: `R/kernel.R`, continuous-G reporting helpers, tests, FDID Notes, the
continuous-G estimator reference, and TODO tracking. `NEWS.md`, `HANDOFF.md`,
and tutorial files are intentionally out of scope for this pass.

- Reworked kernel bootstrap inference toward interflex parity while preserving
  the existing local-linear point estimator.
- `vartype = "bootstrap"` now stores pointwise bootstrap standard errors and
  percentile intervals for both the level curve `mu(g)` and derivative curve
  `delta(g)`.
- Kernel bootstrap now stores separate level and derivative simultaneous bands
  using an interflex-style bootstrap quantile-envelope calibration. These bands
  are documented as quantile-envelope bands, not as fully studentized sup-t
  bands.
- Kernel bootstrap now stores grid-by-grid bootstrap covariance matrices
  `mu_boot_vcov` and `delta_boot_vcov`.
- Scalar event summaries and grid-level contrasts can use bootstrap replicate
  contrasts or the covariance identity `V22 + V11 - 2 V12` instead of treating
  endpoints as independent.
- Kernel bootstrap resampling is row-based by default and cluster-based when
  `fdid_prepare(cluster_label = ...)` supplies a cluster column.
- `plot(type = "curve")` can display kernel level bands as well as derivative
  bands, and uses stored bootstrap pointwise intervals when available.
- Remaining theory caveat: a fully studentized kernel sup-t band remains a
  separate optional extension. The implemented default follows the Practical
  Guide/interflex quantile-envelope procedure.

## 2026-05-29 Grid And Plot Clarification

Scope: `fdid-main` package docs and package-level plot behavior.

- Rechecked `R/kernel.R`, `R/dml.R`, and `R/plot.R` against the local
  `interflex-master` kernel and DML code paths.
- Confirmed that `eval_g` is an estimation grid for `kernel` and `dml_flex`,
  not a plotting-only option.
- Confirmed that FDID kernel runs one local-linear WLS fit per `eval_g` point.
  If a caller supplies five quantiles, the package estimates five curve points.
  If `eval_g = NULL`, FDID estimates the default 50-point trimmed grid.
- Confirmed that FDID `dml_flex` builds orthogonal pseudo-outcomes for all
  observations, then maps those signals to `eval_g`; the stored level and
  derivative curves are still grid-based.
- Compared this with `interflex`: kernel uses `X.eval` or a default 50-point
  `neval` grid, and DML computes pseudo-outcomes/influence values for all
  observations but predicts the CATE on a default 50-point BLP grid.
- Clarified that the important difference is default grid range and any
  user-supplied sparse grid, not whether one package estimates every observed
  moderator/treatment value by default.
- Recorded that curve support overlays should be read as descriptive displays
  of observed `G` within the plotted curve window, not as part of the
  estimation procedure.

### 2026-05-29 Plotting Updates

Today's package-level plotting changes:

- `plot.fdid(type = "curve")` now treats `xlim` as a curve-window argument:
  stored `eval_g` rows outside `xlim` are omitted before drawing.
- Curve support overlays are aligned to the displayed curve window. The
  observed `G` values used for the lower histogram or density are restricted
  to `plot_xlim`, whether that window comes from `xlim` or from the displayed
  `eval_g` range.
- The lower support display supports `Xdistr = "histogram"`, `"hist"`,
  `"density"`, or `"none"`. Histograms now use a finer 80-bin display so dense
  or clustered support is easier to read.
- The support display is descriptive only. It is not used by the estimator and
  does not define the curve grid.
- Curve plots use the fitted object's `alpha` for pointwise intervals when the
  object stores a valid `alpha`, with a 0.05 fallback for older or malformed
  objects.
- `plot(type = "curve")` remains a display method: it reads stored
  `eval_g`/`curve_event` values and does not re-estimate, densify, or
  interpolate additional curve points.

### 2026-05-29 Reporting-Layer Updates

Scope: first-class reporting for Notes-style continuous-G estimands.

- Added exported helper `fdid_contrast()` for level comparisons such as
  `theta(g1) - theta(g0)` for `dml_flex` and `mu(g1) - mu(g0)` for `kernel`.
- Added exported helper `fdid_derivative()` for fixed-`g` derivative values
  such as `theta_prime(g0)` or `delta(g0)`.
- Both helpers operate on stored `eval_g` curve rows and support explicit
  `interpolate = "linear"` or `interpolate = "none"` behavior.
- `fdid_contrast()` uses stored replicate curves when available, then
  conservative band-implied intervals when possible, and only falls back to
  endpoint-independent pointwise intervals with an approximate note.
- Kernel bootstrap fits now store `curve_event$mu_bootstrap` and
  `curve_event$delta_bootstrap`, so level contrasts and derivative values can
  use bootstrap replicate curves when `vartype = "bootstrap"`.
- `dml_flex` fits now store `curve_event$theta_multiplier_reps` and
  `curve_event$delta_multiplier_reps` from the practical multiplier routine.
- Added package tests for replicate storage and reporting helper extraction.
- Updated the continuous-G tutorial and FDID Notes to show and describe the new
  reporting layer.
- Remaining caveat: `dml_flex` multiplier intervals are still practical
  signal-mapping intervals, not final theorem-level DML uniform inference.

## Source Markdown Files Consolidated

Archived source files:

```text
architecture.md
CONTINUOUS_G_IMPLEMENTATION_DETAILS.md
CONTINUOUS_G_METHODS.md
CONTINUOUS_G_TUTORIAL_REORG_PLAN_2026-05-18.md
DML_ESTIMATOR_REVIEW_VERDICT.md
DML_FLEX_SUMMARY.md
PROJECT_CHANGE_LOG.md
TODO.md
TODO_COMPLETION_LOG_2026-05-16.md
ORACLE_SIGNAL_VALIDATION_RESULTS.md
ORACLE_SIGNAL_STRESS_RESULTS.md
NUISANCE_ORACLE_ABLATION_RESULTS.md
```

The `fdid.Rcheck/` copies were treated as generated copies, not independent
source documents.

The archive also contains generated Rcheck copies and older snapshots such as
`RCHECK_00_pkg_src_CONTINUOUS_G_METHODS.md`. Those are preserved for provenance
but the three active markdowns are the source of truth after 2026-05-19.
Some older generated snapshots predate the later `dml_flex` versus
`dml_incremental` correction; where they describe `dml_incremental` as a curve
method, that description is superseded by the current scalar
average-derivative implementation.

## Package Baseline

The original public package implemented six FDID estimators:

```text
ols1
ols2
did
ebal
ipw
aipw
```

The original variance support was:

```text
robust
bootstrap
jackknife
```

The original workflow:

1. `fdid_prepare()` reshapes long panel data to wide unit-level data.
2. `fdid()` computes FDID estimates.
3. `print()`, `summary()`, and `plot()` display results.
4. `fdid_list()` compares multiple fitted objects.

The `NEWS.md` package history records:

- `fdid 0.2.0`: initial public release with six estimators, variance support,
  parallel bootstrap/jackknife, `fdid_prepare()`, `fdid_list()`, S3 methods,
  the built-in `mortality` dataset, and the Xu, Zhao, and Ding (2026)
  reference.
- `fdid 1.0.2`: replaced global `options(warn = -1)` with
  `suppressWarnings()` inside `silent_ebalance()` for CRAN compliance.

The architecture before the continuous-G work was mostly monolithic:

- `R/fdid.R`: main estimator, original method subroutines, original variance
  wrappers, and helper closures.
- `R/fdid_prepare.R`: data reshaping.
- `R/plot.R`: plotting for raw means, dynamic effects, overlap, and later
  continuous-G curves.
- `R/summary.R`: formatted summaries.
- `R/plot.fdid_list.R`: comparison plots.
- `R/print.R`: print method.
- `R/data.r`: mortality data documentation.

Continuous-G work added major implementation code in:

```text
R/kernel.R
R/dml.R
R/fdid_prepare.R
R/plot.R
R/summary.R
R/fdid.R
```

## High-Level Continuous-G Summary

The package now adds five method values to `fdid()`:

```text
method = "kernel"
method = "dml_binary"
method = "dml_plr"
method = "dml_flex"
method = "dml_incremental"
```

Important estimand split:

| Method | Target | Output |
|---|---|---|
| `kernel` | local-linear continuous-G level and derivative curves | curve |
| `dml_plr` | partially linear scalar slope | scalar |
| `dml_flex` | fixed-G level curve `theta(g)` and derivative curve `theta'(g)` | curve |
| `dml_incremental` | observed-population average derivative `E[partial_g mu(G, X)]` | scalar |
| `dml_binary` | binary-G cross-fitted AIPW/IRM contrast | scalar |

The largest conceptual correction was separating `dml_flex` from
`dml_incremental`:

- `dml_flex` owns the flexible level curve and the derivative of that curve.
- `dml_incremental` owns the scalar observed-population average derivative and
  uses its own orthogonal score.
- `dml_incremental` is not obtained by differentiating or averaging the
  `dml_flex` curve.

This correction was made because deriving the scalar incremental target from
the flexible curve can lose the Neyman-orthogonality logic for the scalar
average derivative.

## Timeline

### Early Continuous-G Implementation Phase

Main additions:

- `method = "kernel"` for local-linear continuous-G curves.
- `method = "dml_plr"` for a scalar partially linear DML slope.
- `method = "dml_binary"` for cross-fitted binary-G AIPW/IRM.
- `method = "dml_flex"` for flexible DML level curves.
- Planner guidance and pseudocode were used, but the final API preserved the
  existing `fdid()` entry point.

Key rule preserved:

```text
New estimators enter through fdid(..., method = ...), not separate functions.
```

Existing methods were kept unchanged in behavior.

### 2026-05-04 Validation, Notes, And Visualization Phase

Main additions:

- oracle-signal validation DGP;
- stress tests varying `G | X` density shape and nonlinear/nonsmooth
  `mu(g, x)`;
- nuisance-oracle ablation tests;
- validation figures;
- FDID Notes section describing implementation, validation, caveats, and
  diagnostic conclusions;
- plotting changes distinguishing level and derivative curves.

### 2026-05-10 Planning And TODO Phase

Main additions:

- large TODO plan covering inference cleanup, diagnostics, density extensions,
  tutorial organization, validation, release hygiene, and implementation order;
- detailed Continuous G tutorial reorganization plan centered on estimator
  choice, mortality example, HD expansion, learner sensitivity, and synthetic
  sandbox.

### 2026-05-16 Completion Pass

Completed from the May 10 TODOs:

- input validation;
- build hygiene;
- cluster-aware scalar DML inference;
- DML diagnostics;
- `dml_flex` practical multiplier bands;
- plot and summary changes;
- xgboost backend fix;
- tutorial reorganization;
- mortality support diagnostics;
- HD and learner sensitivity examples;
- synthetic examples;
- FDID Notes update;
- package verification.

This pass also addressed the earlier review labels:

- Mapping Extensions: added `signal_map = "local_poly"`, `"kernel"`,
  `"spline"`, and `"gam"`;
- Density Extensions: added `density_method = "residual_kde"`,
  `"location_scale"`, and `"local_kde"`;
- Summary Change: updated `summary.fdid()` to distinguish `dml_flex` curve
  output from scalar compatibility output and to report DML diagnostics and
  curve-band metadata.

### 2026-05-18 Tutorial Reorganization Review

The tutorial organization was reviewed against:

- `fect`: https://yiqingxu.org/packages/fect/index.html
- `interflex`: https://yiqingxu.org/packages/interflex/
- `interflex` continuous treatment tutorial:
  https://yiqingxu.org/packages/interflex/articles/continuous.html
- `interflex` DML tutorial:
  https://yiqingxu.org/packages/interflex/articles/dml.html
- `fect` cheatsheet:
  https://yiqingxu.org/packages/fect/aa-cheatsheet.html

Main lesson:

The Continuous G tutorial should be an applied teaching workflow, not a
combined methods memo, validation appendix, and development notebook.

Lessons copied from the 2026-05-18 plan:

- separate the manual from the lesson;
- teach by empirical design type;
- put code and the first result before deep theory;
- use a compact cheatsheet after the reader has seen the workflow;
- start from a concrete data object and a substantive question;
- fit one model, show one table or plot, interpret it, then add diagnostics
  and options;
- do not show all learners, all density methods, all simulations, and all
  validation evidence in the first visible path.

Recommended structure:

```text
3 Continuous G
  3.1 Continuous G in one page
  3.2 Load the local package and mortality data
  3.3 Inspect treatment support
  3.4 Scalar benchmarks
  3.5 Kernel level and derivative curves
  3.6 Flexible DML level and derivative curves
  3.7 Scalar incremental DML
  3.8 Learner and HD expansion choices
  3.9 One known-truth synthetic sandbox
  3.10 Method choice cheatsheet

4 Continuous G advanced examples
  learner availability, HD expansion, learner grids, stress examples

5 Continuous G technical appendix
  estimands, signals, cross-fitting, nuisances, density, inference, validation
```

Concrete rewrite actions from the plan:

1. Start `03-continuous.Rmd` with a short continuous-G roadmap, not local
   install details.
2. Move local package reinstall/re-library instructions to a compact callout
   or appendix.
3. Keep mortality as the first real example.
4. Put support diagnostics before any curve estimator.
5. Present scalar benchmarks in one table.
6. Present kernel as the first curve method.
7. Present `dml_flex` as the flexible level-curve method, with derivative plot
   secondary.
8. Present `dml_incremental` as the scalar average derivative method.
9. Replace the large learner comparison in the main tutorial with a two-by-two
   demonstration: raw versus HD expansion, simple versus flexible learner.
10. Move full learner and HD grids to an advanced page or optional section.
11. Keep only one known-truth synthetic example in the main tutorial.
12. Move stress tests and validation plots to advanced or technical pages.
13. End with a method-choice cheatsheet.

Acceptance criteria from the plan:

- a new reader can run a first continuous-G model within 10 minutes;
- the first empirical plot appears before any long theory block;
- each estimator has one primary result object;
- the tutorial distinguishes level curves, derivative curves, and scalar
  derivative estimands;
- `dml_flex` is not presented as if its scalar compatibility summary were the
  main target;
- `dml_incremental` is described as a separate scalar average derivative target;
- optional learner failures do not break the default render;
- advanced validation material remains available but does not dominate the main
  teaching path.

### 2026-05-19 Markdown Consolidation

The project markdowns were consolidated into three active documents:

```text
TODO.md
CONTINUOUS_G_PROJECT_LOG.md
CONTINUOUS_G_ESTIMATOR_REFERENCE.md
```

Historical markdowns were moved to:

```text
historical_mds_2026-05-19/
```

`HANDOFF.md` was then added as a compact restart entry point. It intentionally
does not replace the three active source-of-truth documents.

### 2026-05-19 Documentation Audit And Code Cross-Check

A follow-up audit checked the three consolidated docs against:

```text
R/fdid.R
R/kernel.R
R/dml.R
R/fdid_prepare.R
R/plot.R
R/summary.R
R/print.R
R/fdid_list.R
R/plot.fdid_list.R
DESCRIPTION
NAMESPACE
.Rbuildignore
inst/
historical_mds_2026-05-19/
```

Important audit outcomes:

- the workspace root and `fdid-main/` are not Git repositories, so this was a
  file-level audit rather than a git-diff audit;
- validation scripts and saved outputs are still present under `inst/`;
- the active estimator reference was expanded with code-level details for
  dispatch, dependency checks, missing-data handling, learner aliases/backends,
  density floors, finite-difference rules, object fields, plotting behavior,
  and S3 helpers;
- `TODO.md` was expanded with a dated crosswalk preserving TODOs from chat and
  archived markdowns;
- this log was updated to record the audit itself and correct saved validation
  figure filenames;
- `.Rbuildignore` now ignores `historical_mds_2026-05-19/` so the markdown
  archive is not accidentally included in package builds.

Code-level caveats surfaced by the audit:

- `fdid()` checks for `estimatr` and `car` before all method dispatch,
  including DML paths;
- continuous-G dispatchers do not include cluster column `c` in their initial
  complete-case filter;
- scalar DML SEs are cluster-aware, but kernel and `dml_flex` curve inference
  are not;
- `dml_flex` curve diagnostics currently report no cluster count because the
  diagnostic helper is called with `cluster = NULL`;
- public `boot = 500` was capped internally at 100 for `dml_flex` practical
  multiplier bands at the time of this audit; this was superseded by the
  2026-05-30 DML inference update, where `boot` controls the requested draws;
- `plot(type = "curve")` draws 95 percent pointwise intervals with
  `qnorm(0.975)`, regardless of the fit's original `alpha`;
- continuous-G `raw_means` are median-split displays, not estimator targets.

## API Changes

### Method Options

Current `fdid()` methods:

```r
c(
  "ols1", "ols2", "did", "ebal", "aipw", "ipw",
  "kernel", "dml_binary", "dml_plr", "dml_flex", "dml_incremental"
)
```

### Continuous-G Arguments

| Argument | Methods | Purpose |
|---|---|---|
| `eval_g` | `kernel`, `dml_flex` | G-grid for curve estimation |
| `h0` | `kernel` | baseline bandwidth |
| `K_folds` | `kernel` | bandwidth cross-validation folds |
| `boot` | `kernel`, `dml_flex` | kernel bootstrap or DML multiplier band replications |
| `alpha` | continuous-G methods | CI level |
| `trim` | `kernel`, `dml_flex` | tail trimming for default G grid |
| `learner` | DML methods | nuisance learner |
| `K` | DML methods | cross-fitting folds |
| `S` | DML methods | repeated random splits |
| `signal_map` | `dml_flex` | signal-to-G mapping |
| `map_degree` | `dml_flex` | local polynomial or spline degree |
| `map_df` | `dml_flex` | spline or GAM basis size |
| `density_method` | `dml_flex`, `dml_incremental` | conditional-density estimator |

### Variance Support

| Method family | Supported `vartype` |
|---|---|
| `ols1`, `ols2`, `did`, `ebal`, `ipw`, `aipw` | `robust`, `bootstrap`, `jackknife` |
| `kernel` | `robust`, `bootstrap` |
| `dml_binary`, `dml_plr`, `dml_flex`, `dml_incremental` | `robust` only |

DML bootstrap and jackknife requests are rejected intentionally.

### Namespace And Public Functions

Current exported functions:

```text
fdid
fdid_prepare
fdid_expand_covariates
fdid_list
```

Current S3 registrations:

```text
plot.fdid
plot.fdid_list
print.fdid
summary.fdid
```

The continuous-G estimators remain internal implementation paths. They are not
exported as separate top-level functions.

## Implementation Changes

### Input Validation

Added early validation in `fdid()` and repeated sample-size-aware validation in
continuous-G dispatchers.

Validated:

```text
alpha > 0 and alpha < 1
trim >= 0 and trim < 0.5
eval_g finite and nonempty when supplied
h0 positive when supplied
boot >= 1
K_folds >= 2
K >= 2
S >= 1
K_folds <= n after missing-data handling
K <= n after missing-data handling
```

Previous failure modes fixed:

| Call pattern | Old behavior |
|---|---|
| `method = "dml_plr", K = 0` | returned `Estimate = 0`, `Std.Error = Inf` |
| `method = "dml_plr", S = 0` | returned all `NA` |
| `method = "dml_flex", K = 0` | returned zero curve summaries |
| `method = "kernel", K_folds = 0` | returned a kernel estimate |
| `alpha = 2` | returned invalid confidence intervals |

### Build Hygiene

Changes:

- removed local install metadata from `DESCRIPTION`;
- expanded `.Rbuildignore`;
- kept generated checks, tutorial renders, and validation result folders out of
  package builds unless intentionally shipped.
- on 2026-05-19, added `^historical_mds_2026-05-19$` to `.Rbuildignore` after
  archiving the historical markdowns.

### Kernel Estimator

Implemented in `R/kernel.R`.

Main features:

- requires continuous `G` with more than two unique values;
- builds event-period change score:

```text
DeltaY_i = mean(Y_i,t over tr_period) - Y_i,ref
```

- fits local-linear regressions over `G`;
- includes covariate adjustment through `x1`, `x2`, ... and interactions with
  local `G - g0`;
- uses adaptive bandwidths based on a pilot marginal density of `G`;
- chooses baseline bandwidth by least-squares cross-validation when `h0` is
  not supplied;
- stores level curve `mu_hat(g)` and derivative curve `delta_hat(g)`;
- stores pointwise robust SEs for both;
- supports bootstrap bands for the derivative curve when
  `vartype = "bootstrap"`;
- returns scalar `est$event` as an interval-average summary over the displayed
  grid.

### DML Learner Registry

DML methods share learner helpers in `R/dml.R`.

Supported learners:

| User value | Backend |
|---|---|
| `linear`, `lm`, `ols`, `default` | `lm()` or binomial `glm()` |
| `glmnet`, `lasso` | `glmnet::cv.glmnet(alpha = 1)` |
| `ridge` | `glmnet::cv.glmnet(alpha = 0)` |
| `elasticnet`, `enet` | `glmnet::cv.glmnet(alpha = 0.5)` |
| `ranger`, `rf`, `randomforest`, `forest` | `ranger::ranger()` |
| `grf`, `regressionforest`, `probabilityforest` | `grf::regression_forest()` |
| `gam`, `mgcv` | `mgcv::gam()` |
| `nnet`, `nn`, `neuralnetwork` | `nnet::nnet()` |
| `xgboost`, `xgb`, `boosting` | `xgboost::xgb.train()` |

Optional packages are checked lazily when the learner is requested.

The xgboost backend was updated to use `xgb.train()`, compatible with xgboost
3.x.

### `dml_binary`

Target:

```text
E[mu_1(X) - mu_0(X)]
```

Implementation:

- requires binary `G` in `{0, 1}`;
- cross-fits `mu0(X)`, `mu1(X)`, and propensity `pi(X)`;
- uses an AIPW/IRM score;
- clips propensity predictions to `[0.01, 0.99]`;
- repeats across `S` splits and uses median aggregation;
- robust SE uses orthogonal score variance plus repeated-split dispersion;
- scalar DML SEs are cluster-aware when `fdid_prepare(cluster_label = ...)`
  supplies a cluster column.

### `dml_plr`

Target:

```text
DeltaY = beta G + h(X) + U
```

Implementation:

- cross-fits `E[DeltaY | X]` and `E[G | X]`;
- residualizes outcome and treatment;
- estimates the Robinson slope:

```text
beta_hat = sum((G - m_hat(X)) * (DeltaY - l_hat(X))) /
           sum((G - m_hat(X))^2)
```

- uses score:

```text
psi_i = V_i * (U_i - beta_hat V_i)
```

- divides by Jacobian `mean(V_i^2)` for variance;
- repeats across `S` splits and uses median aggregation;
- scalar DML SEs are cluster-aware when clusters are supplied.

### `dml_flex`

Target:

```text
theta(g) = E_X[mu(g, X)]
mu(g, x) = E[DeltaY | G = g, X = x]
```

Implementation:

- requires continuous `G`;
- constructs default `eval_g` as a trimmed grid over `G` when not supplied;
- cross-fits outcome nuisance `mu(G, X)`;
- cross-fits conditional density `s(G | X)`;
- estimates marginal density `omega(G)` within each training fold;
- builds a cross-fitted orthogonal pseudo-outcome:

```text
xi_i =
  ((DeltaY_i - mu_hat(G_i, X_i)) / s_hat(G_i | X_i)) * omega_hat(G_i)
  + m_hat(G_i)
```

- `m_hat(G_i)` is estimated by averaging the fitted outcome regression at
  `G_i` over held-out covariates, with a leave-one-out adjustment inside the
  held-out fold;
- maps `xi_i` back to `G` using `signal_map`;
- stores `theta_hat(g)` and `delta_hat(g) = theta_hat'(g)`;
- stores pointwise CIs and practical multiplier bands;
- treats scalar `est$event` as a grid-average compatibility summary, not the
  main estimand.

Signal-map options:

| `signal_map` | Implementation |
|---|---|
| `local_poly` | local polynomial map, default degree 2 |
| `kernel` | local-linear map |
| `spline` | B-spline linear model |
| `gam` | `mgcv::gam()` smooth |

### `dml_incremental`

Target:

```text
E[partial_g mu(G, X)]
```

Implementation:

- requires continuous `G`;
- ignores `eval_g` and `signal_map`;
- cross-fits `mu(G, X)` and `s(G | X)`;
- computes finite-difference estimates of `partial_g mu(G, X)`;
- computes finite-difference estimates of `partial_g log s(G | X)`;
- builds the scalar average-derivative orthogonal signal:

```text
phi_i =
  partial_g mu_hat(G_i, X_i)
  - partial_g log s_hat(G_i | X_i) * (DeltaY_i - mu_hat(G_i, X_i))
```

- averages held-out `phi_i`;
- repeats across `S` splits and uses median aggregation;
- stores the signal, split estimates, density method, and diagnostics;
- scalar DML SEs are cluster-aware when clusters are supplied.

### Conditional Density Methods

Used by `dml_flex` and `dml_incremental`.

| Method | Description |
|---|---|
| `residual_kde` | fit `E[G | X]`, KDE residuals; location-shift approximation |
| `location_scale` | fit conditional mean and log residual variance, KDE standardized residuals |
| `local_kde` | product Gaussian kernel over scaled `X` and Gaussian KDE over `G` |

Validation showed that conditional density is a first-order practical nuisance.
`location_scale` behaved better than `residual_kde` in heteroskedastic and
mixture settings; `residual_kde` was competitive in heavy-tailed curve tests.

### High-Dimensional Covariate Expansion

Implemented in `fdid_prepare()` and `fdid_expand_covariates()`.

Options:

| Argument | Meaning |
|---|---|
| `basis_type = "none"` | raw covariates |
| `basis_type = "polynomial"` | powers up to `poly_degree` |
| `basis_type = "bspline"` | B-spline basis |
| `include_interactions = TRUE` | pairwise interactions among expanded columns |
| `standardize_basis = TRUE` | center and scale expanded columns |

The expansion replaces the prepared `x1`, `x2`, ... columns. It does not alter
the treatment/factor `G`.

### Plotting

`plot.fdid()` now supports:

```r
plot(result, type = "curve", curve = "level")
plot(result, type = "curve", curve = "derivative")
```

Curve availability:

| Method | Level curve | Derivative curve | `curve = "auto"` |
|---|---|---|---|
| `kernel` | `mu_hat(g)` | `delta_hat(g)` | derivative |
| `dml_flex` | `theta_hat(g)` | `theta_hat'(g)` | level |
| `dml_incremental` | none | none | none |

For continuous-G curves:

- pointwise CIs are drawn in blue;
- available method-specific bands are drawn in gray;
- kernel derivative bootstrap bands are used when available;
- `dml_flex` practical multiplier bands are used when available.
- plotted curve points are the stored `eval_g` rows; the plotting method does
  not add new estimation points;
- lower support overlays are descriptive observed-`G` histograms or densities
  restricted to the displayed curve window.

### Summary

`summary.fdid()` now:

- prints continuous-G curve blocks for `kernel` and `dml_flex`;
- prints `dml_incremental` scalar-signal information;
- prints learner, `K`, `S`, signal map, and density method when applicable;
- prints the `dml_flex` simultaneous-band method when stored;
- prints DML diagnostics when stored;
- keeps `dml_flex` scalar `est$event` explicitly secondary.

## Inference Log

### Original Methods

Original methods support:

- robust SEs;
- bootstrap;
- jackknife.

Cluster support is available where the original robust regression paths use the
cluster column.

### Kernel

Kernel inference:

- pointwise robust sandwich SEs for level and derivative curves;
- pointwise bootstrap SEs and percentile intervals for level and derivative
  curves when `vartype = "bootstrap"`;
- level and derivative curve bands through either pointwise normal intervals or
  interflex-style bootstrap quantile-envelope bands;
- bootstrap curve covariance matrices for level and derivative curves;
- scalar event summary is an interval-average over the displayed grid and uses
  bootstrap replicate contrasts or covariance-aware endpoint variance when
  available;
- bootstrap resampling is cluster-aware when a prepared cluster column exists.

Caveat:

The implemented bootstrap simultaneous band is a quantile-envelope band over
the evaluation grid. A fully studentized sup-t kernel band remains a separate
theory-aligned extension.

### Scalar DML

Scalar DML inference:

```text
se = sqrt(median(var_s + (tau_s - median(tau_s))^2))
```

where `var_s` is split-specific score variance and the second term captures
repeated-split dispersion.

Cluster-aware scalar DML:

- aggregates influence contributions by cluster;
- uses a finite-cluster correction;
- falls back to unit-level variance when clusters are unavailable or unusable.

Applies to:

```text
dml_binary
dml_plr
dml_incremental
```

### `dml_flex`

`dml_flex` stores:

```text
theta_hat
delta_hat
se_theta
se_delta
theta_ci_lower
theta_ci_upper
delta_ci_lower
delta_ci_upper
theta_band_lower
theta_band_upper
delta_band_lower
delta_band_upper
simultaneous_crit_theta
simultaneous_crit_delta
simultaneous_band_method
```

Interpretation:

- pointwise CIs are direct second-stage mapping intervals with repeated-split
  dispersion folded in;
- multiplier bands are practical signal-mapping bands;
- bands are useful diagnostics but not yet final theorem-level DML uniform
  bands.

Important open issue:

Level comparisons such as `theta(g1) - theta(g0)` are now exposed through
`fdid_contrast()`, which first uses stored curve covariance matrices when
available, then stored replicate curves, then clearly labeled fallback
intervals. The remaining open research-level issue is full theorem-level
efficient covariance or influence-function storage for non-BLP maps and all
nuisance stages; `signal_map = "blp_spline"` already stores finite-grid level
and derivative covariance matrices.

## Tutorial Updates

The tutorial was updated to:

- include local RStudio package workflow;
- add mortality continuous-G support diagnostics;
- explain that the main curve grid uses interior positive support;
- explain that selected quantiles are landmarks, not bins;
- organize mortality examples into scalar benchmarks, kernel curves,
  `dml_flex` curves, scalar `dml_incremental`, sensitivity, and a binary
  companion;
- add HD expansion and learner sensitivity for both `dml_flex` and
  `dml_incremental`;
- add known-truth synthetic examples for smooth nonlinear response, nonlinear
  confounding plus HD expansion, heteroskedastic `G | X`, kinked outcome
  functions, and learner/HD grids;
- add `FDID_TUTORIAL_HEAVY=true` as the switch for full optional learner grids.

Remaining tutorial concern:

The Continuous G tutorial still needs clearer pedagogy. The best structure is a
main applied workflow plus advanced and technical appendices.

## FDID Notes Updates

FDID Notes were updated to describe:

- all continuous-G estimators;
- cross-fitting;
- nuisance functions;
- density methods;
- validation tests;
- diagnostic conclusions;
- caveats;
- method choice guidance.

Future Notes updates should add a full discussion of:

- level comparisons `theta(g1) - theta(g0)`;
- derivative values `theta_prime(g0)`;
- which SE/CI procedures are available now;
- which are approximate or diagnostic;
- what remains to be implemented.

## Validation Results

### Early Monte Carlo Validation

Validation script:

```text
inst/simulate_mc.R
```

Older generated documentation recorded the following early Monte Carlo checks
with `n = 300` and 100 simulations. These were preliminary regression checks,
not the later oracle-signal validation battery.

DGP 1, linear moderation:

```text
DeltaY = 2G + 0.5X + error
true delta = 2
true beta = 2
theta(g) = 2g
```

| Criterion | Result | Target |
|---|---:|---:|
| Kernel RMSE on derivative curve | 0.277 | < 0.4 |
| PLR absolute bias | 0.023 | < 0.2 |
| PLR 95% coverage | 0.940 | >= 0.85 |
| Flex RMSE on level curve | 0.074 | < 0.5 |

DGP 5, null effect:

| Criterion | Result | Target |
|---|---:|---:|
| Kernel rejection rate | 0.060 | <= 0.12 |
| PLR rejection rate | 0.040 | <= 0.12 |
| DML binary rejection rate | 0.030 | <= 0.12 |

All early targets were recorded as met. Later validation work superseded these
as the primary evidence.

### Smooth Oracle-Signal Validation

Validation script:

```text
inst/oracle_signal_validation.R
```

DGP:

```text
mu(g, x) = 0.40 + 0.50 x1 - 0.25 x2 + 0.80 g + 0.25 g^2 + sin(0.70 g)
G | X ~ N(0.40 x1 - 0.30 x2, 0.70^2)
partial_g mu(g, x) = 0.80 + 0.50 g + 0.70 cos(0.70 g)
partial_g log s(g | x) = -(g - 0.40 x1 + 0.30 x2) / 0.70^2
```

Settings:

```text
Master seed: 20260503
Evaluation grid: -0.892, -0.615, -0.390, -0.190, 0.000, 0.190, 0.390, 0.615, 0.892
Estimated-nuisance learner: gam
Signal map: spline
Density method: residual_kde
```

Oracle signal sanity check:

| n | reps | level_rmse | level_coverage | derivative_rmse | derivative_coverage | scalar_bias | scalar_rmse | scalar_coverage |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 500 | 120 | 0.0536 | 0.9528 | 0.1609 | 0.9324 | -0.0009 | 0.0298 | 0.9917 |
| 1500 | 120 | 0.0342 | 0.9454 | 0.1034 | 0.9176 | 0.0024 | 0.0164 | 0.9917 |

Estimated-nuisance comparison:

| n | reps | flex_failures | incremental_failures | oracle_level_rmse | flex_level_rmse | flex_level_vs_oracle_rmse | flex_level_coverage | oracle_derivative_rmse | flex_derivative_rmse | flex_derivative_vs_oracle_rmse | flex_derivative_coverage | oracle_scalar_bias | oracle_scalar_rmse | oracle_scalar_coverage | incremental_scalar_bias | incremental_scalar_rmse | incremental_scalar_coverage | incremental_vs_oracle_bias | incremental_vs_oracle_rmse |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 400 | 40 | 0 | 0 | 0.0588 | 0.0640 | 0.0352 | 0.9722 | 0.1925 | 0.1837 | 0.1318 | 0.9583 | -0.0061 | 0.0300 | 1.0000 | -0.0061 | 0.0326 | 0.9750 | 0.0001 | 0.0146 |
| 800 | 40 | 0 | 0 | 0.0404 | 0.0401 | 0.0133 | 0.9694 | 0.1184 | 0.1131 | 0.0508 | 0.9472 | 0.0007 | 0.0239 | 0.9750 | 0.0005 | 0.0254 | 0.9500 | -0.0001 | 0.0090 |

Interpretation:

- oracle signal RMSE falls as `n` increases;
- estimated-nuisance RMSE is generally larger than oracle RMSE;
- scalar `dml_incremental` targets `E[partial_g mu(G, X)]`;
- derivative-curve comparisons belong to `dml_flex`, not
  `dml_incremental`.

### Stress Tests

Validation script:

```text
inst/oracle_signal_stress_tests.R
```

Settings:

```text
Master seed: 20260504
Oracle reps per scenario/n: 80
Estimated reps per scenario/n/density: 12
Estimated learner: gam
Signal map: spline
Map df: 7
Density methods: residual_kde, location_scale
```

Scenarios:

- `heteroskedastic_interactions`: normal `G | X` with X-dependent scale and
  smooth G-by-X interactions in `mu(g, x)`.
  Evaluation grid:
  `-0.858, -0.618, -0.414, -0.231, -0.051, 0.130, 0.333, 0.561, 0.866`.
- `mixture_kink`: two-component normal mixture `G | X` with X-dependent
  weights and kinks in `mu(g, x)`.
  Evaluation grid:
  `-0.838, -0.619, -0.430, -0.254, -0.072, 0.123, 0.349, 0.641, 1.023`.
- `heavy_tail_interactions`: scaled t(3) residuals in `G | X` with nonlinear
  trigonometric interactions in `mu(g, x)`.
  Evaluation grid:
  `-0.932, -0.619, -0.388, -0.190, -0.004, 0.186, 0.396, 0.633, 0.945`.

Oracle signal results:

| scenario | n | reps | level_rmse | level_coverage | derivative_rmse | derivative_coverage | scalar_bias | scalar_rmse | scalar_coverage |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| heavy_tail_interactions | 600 | 80 | 0.0400 | 0.9222 | 0.0972 | 0.9028 | 0.0009 | 0.0240 | 0.9875 |
| heavy_tail_interactions | 1200 | 80 | 0.0306 | 0.9056 | 0.0873 | 0.7681 | 0.0010 | 0.0171 | 0.9625 |
| heteroskedastic_interactions | 600 | 80 | 0.0484 | 0.9389 | 0.1578 | 0.9458 | 0.0020 | 0.0286 | 1.0000 |
| heteroskedastic_interactions | 1200 | 80 | 0.0372 | 0.9444 | 0.1239 | 0.9097 | -0.0017 | 0.0200 | 0.9875 |
| mixture_kink | 600 | 80 | 0.0418 | 0.9278 | 0.1865 | 0.8278 | 0.0001 | 0.0268 | 0.9875 |
| mixture_kink | 1200 | 80 | 0.0327 | 0.9319 | 0.1686 | 0.7597 | -0.0007 | 0.0170 | 1.0000 |

Estimated-nuisance results:

| scenario | density_method | n | reps | flex_failures | incremental_failures | oracle_level_rmse | flex_level_rmse | flex_level_vs_oracle_rmse | flex_level_coverage | oracle_derivative_rmse | flex_derivative_rmse | flex_derivative_vs_oracle_rmse | flex_derivative_coverage | incremental_scalar_bias | incremental_scalar_rmse | incremental_scalar_coverage | incremental_vs_oracle_rmse |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| heavy_tail_interactions | location_scale | 500 | 12 | 0 | 0 | 0.0415 | 0.0614 | 0.0361 | 0.9537 | 0.0903 | 0.1489 | 0.1133 | 0.9352 | 0.0157 | 0.0307 | 1.0000 | 0.0210 |
| heavy_tail_interactions | location_scale | 1000 | 12 | 0 | 0 | 0.0289 | 0.0350 | 0.0211 | 0.9537 | 0.0759 | 0.0855 | 0.0655 | 0.9630 | -0.0018 | 0.0320 | 0.9167 | 0.0170 |
| heavy_tail_interactions | residual_kde | 500 | 12 | 0 | 0 | 0.0420 | 0.0446 | 0.0377 | 0.9722 | 0.0864 | 0.1040 | 0.1063 | 0.9630 | 0.0099 | 0.0282 | 1.0000 | 0.0210 |
| heavy_tail_interactions | residual_kde | 1000 | 12 | 0 | 0 | 0.0262 | 0.0304 | 0.0204 | 0.9722 | 0.0827 | 0.0994 | 0.0570 | 0.8796 | -0.0033 | 0.0233 | 1.0000 | 0.0190 |
| heteroskedastic_interactions | location_scale | 500 | 12 | 0 | 0 | 0.1112 | 0.1658 | 0.1088 | 0.9444 | 0.3788 | 0.5509 | 0.4129 | 0.9537 | 0.0241 | 0.0601 | 1.0000 | 0.0494 |
| heteroskedastic_interactions | location_scale | 1000 | 12 | 0 | 0 | 0.0363 | 0.0765 | 0.0625 | 0.9074 | 0.1124 | 0.2307 | 0.2054 | 0.8241 | 0.0332 | 0.0443 | 0.9167 | 0.0343 |
| heteroskedastic_interactions | residual_kde | 500 | 12 | 0 | 0 | 0.0432 | 0.1103 | 0.1027 | 0.9537 | 0.1132 | 0.4211 | 0.4112 | 0.8981 | 0.1439 | 0.1524 | 0.0833 | 0.1528 |
| heteroskedastic_interactions | residual_kde | 1000 | 12 | 0 | 0 | 0.0659 | 0.1061 | 0.1086 | 0.8148 | 0.1784 | 0.3229 | 0.3555 | 0.6667 | 0.1461 | 0.1484 | 0.0000 | 0.1331 |
| mixture_kink | location_scale | 500 | 12 | 0 | 0 | 0.0427 | 0.0783 | 0.0690 | 1.0000 | 0.2008 | 0.3263 | 0.2638 | 0.9074 | 0.0042 | 0.0351 | 1.0000 | 0.0422 |
| mixture_kink | location_scale | 1000 | 12 | 0 | 0 | 0.0384 | 0.0681 | 0.0501 | 0.9444 | 0.1718 | 0.2825 | 0.2031 | 0.8056 | -0.0005 | 0.0334 | 1.0000 | 0.0205 |
| mixture_kink | residual_kde | 500 | 12 | 0 | 0 | 0.0548 | 0.0640 | 0.0444 | 0.8611 | 0.2081 | 0.2574 | 0.1743 | 0.7778 | -0.0715 | 0.0790 | 0.6667 | 0.0823 |
| mixture_kink | residual_kde | 1000 | 12 | 0 | 0 | 0.0352 | 0.0480 | 0.0362 | 0.8519 | 0.1599 | 0.2028 | 0.1322 | 0.7593 | -0.0689 | 0.0759 | 0.3333 | 0.0681 |

Stress-test interpretation:

- `flex_*` columns evaluate `dml_flex` level and derivative curves.
- `incremental_*` columns evaluate scalar `dml_incremental`.
- `*_vs_oracle_*` columns compare estimated-nuisance fits to oracle-signal
  fits on the same DGP.
- Higher RMSE than the smooth baseline is expected.
- Density choice matters strongly for scalar incremental targets in
  heteroskedastic and mixture settings.

### Nuisance-Oracle Ablation

Validation script:

```text
inst/nuisance_oracle_ablation.R
```

Settings:

```text
Master seed: 20260505
Replications per scenario/n: 24
Sample sizes: 600, 1000
Learner for estimated nuisances: gam
Signal map: spline
Map df: 7
Estimated density methods: residual_kde, location_scale
```

Ablation layers:

| Layer | Meaning |
|---|---|
| `oracle_all` | true outcome and true density nuisances |
| `estimated_outcome_oracle_density` | estimated outcome, true density |
| `oracle_outcome_estimated_density` | true outcome, estimated density |
| `estimated_all` | estimated outcome and estimated density |

Ablation scenario evaluation grids:

- `heteroskedastic_interactions`:
  `-0.865, -0.625, -0.424, -0.250, -0.071, 0.120, 0.330, 0.561, 0.866`.
- `mixture_kink`:
  `-0.840, -0.618, -0.426, -0.251, -0.068, 0.124, 0.353, 0.636, 1.025`.
- `heavy_tail_interactions`:
  `-0.942, -0.636, -0.404, -0.204, -0.018, 0.174, 0.372, 0.609, 0.920`.

Selected complete summary:

| scenario | n | layer | density_method | reps | failures | level_rmse | level_coverage | derivative_rmse | derivative_coverage | scalar_bias | scalar_rmse | scalar_coverage | level_vs_oracle_rmse | derivative_vs_oracle_rmse | scalar_vs_oracle_rmse |
|---|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| heavy_tail_interactions | 600 | estimated_all | location_scale | 24 | 0 | 0.0545 | 0.9352 | 0.1589 | 0.8704 | -0.0031 | 0.0318 | 1.0000 | 0.0393 | 0.1196 | 0.0249 |
| heavy_tail_interactions | 600 | estimated_all | residual_kde | 24 | 0 | 0.0550 | 0.9398 | 0.1445 | 0.8750 | -0.0014 | 0.0230 | 1.0000 | 0.0407 | 0.1164 | 0.0221 |
| heavy_tail_interactions | 600 | estimated_outcome_oracle_density | oracle | 24 | 0 | 0.0460 | 0.9259 | 0.1326 | 0.8704 | -0.0018 | 0.0248 | 0.9583 | 0.0293 | 0.0946 | 0.0146 |
| heavy_tail_interactions | 600 | oracle_all | oracle | 24 | 0 | 0.0359 | 0.9444 | 0.0897 | 0.9167 | -0.0036 | 0.0218 | 1.0000 | 0.0000 | 0.0000 | 0.0000 |
| heavy_tail_interactions | 600 | oracle_outcome_estimated_density | location_scale | 24 | 0 | 0.0426 | 0.9583 | 0.1059 | 0.9537 | -0.0031 | 0.0212 | 1.0000 | 0.0160 | 0.0468 | 0.0140 |
| heavy_tail_interactions | 600 | oracle_outcome_estimated_density | residual_kde | 24 | 0 | 0.0364 | 0.9537 | 0.0917 | 0.9444 | -0.0043 | 0.0208 | 1.0000 | 0.0091 | 0.0257 | 0.0113 |
| heavy_tail_interactions | 1000 | estimated_all | location_scale | 24 | 0 | 0.0380 | 0.9306 | 0.1102 | 0.8333 | -0.0055 | 0.0279 | 0.9583 | 0.0214 | 0.0654 | 0.0185 |
| heavy_tail_interactions | 1000 | estimated_all | residual_kde | 24 | 0 | 0.0403 | 0.9120 | 0.1181 | 0.7870 | -0.0068 | 0.0191 | 1.0000 | 0.0221 | 0.0648 | 0.0132 |
| heavy_tail_interactions | 1000 | estimated_outcome_oracle_density | oracle | 24 | 0 | 0.0423 | 0.8935 | 0.1349 | 0.7546 | -0.0067 | 0.0229 | 0.9583 | 0.0220 | 0.0746 | 0.0078 |
| heavy_tail_interactions | 1000 | oracle_all | oracle | 24 | 0 | 0.0335 | 0.8704 | 0.1009 | 0.7824 | -0.0034 | 0.0204 | 0.9583 | 0.0000 | 0.0000 | 0.0000 |
| heavy_tail_interactions | 1000 | oracle_outcome_estimated_density | location_scale | 24 | 0 | 0.0334 | 0.8704 | 0.0971 | 0.8009 | -0.0052 | 0.0202 | 1.0000 | 0.0079 | 0.0226 | 0.0071 |
| heavy_tail_interactions | 1000 | oracle_outcome_estimated_density | residual_kde | 24 | 0 | 0.0344 | 0.8843 | 0.1030 | 0.7685 | -0.0041 | 0.0209 | 1.0000 | 0.0084 | 0.0252 | 0.0079 |
| heteroskedastic_interactions | 600 | estimated_all | location_scale | 24 | 0 | 0.0946 | 0.8889 | 0.3032 | 0.8889 | 0.0350 | 0.0476 | 0.9583 | 0.0779 | 0.2632 | 0.0379 |
| heteroskedastic_interactions | 600 | estimated_all | residual_kde | 24 | 0 | 0.1247 | 0.8981 | 0.3976 | 0.8704 | 0.1308 | 0.1368 | 0.0417 | 0.1142 | 0.3731 | 0.1278 |
| heteroskedastic_interactions | 600 | estimated_outcome_oracle_density | oracle | 24 | 0 | 0.0982 | 0.9074 | 0.3483 | 0.8657 | 0.0091 | 0.0360 | 1.0000 | 0.0815 | 0.2913 | 0.0303 |
| heteroskedastic_interactions | 600 | oracle_all | oracle | 24 | 0 | 0.0429 | 0.9630 | 0.1469 | 0.9352 | 0.0065 | 0.0314 | 1.0000 | 0.0000 | 0.0000 | 0.0000 |
| heteroskedastic_interactions | 600 | oracle_outcome_estimated_density | location_scale | 24 | 0 | 0.0440 | 0.9722 | 0.1461 | 0.9676 | 0.0070 | 0.0292 | 1.0000 | 0.0214 | 0.0786 | 0.0099 |
| heteroskedastic_interactions | 600 | oracle_outcome_estimated_density | residual_kde | 24 | 0 | 0.0430 | 0.9722 | 0.1452 | 0.9769 | 0.0007 | 0.0310 | 1.0000 | 0.0264 | 0.0948 | 0.0147 |
| heteroskedastic_interactions | 1000 | estimated_all | location_scale | 24 | 0 | 0.0653 | 0.9352 | 0.2139 | 0.8935 | 0.0234 | 0.0382 | 0.9167 | 0.0548 | 0.1805 | 0.0323 |
| heteroskedastic_interactions | 1000 | estimated_all | residual_kde | 24 | 0 | 0.1214 | 0.8380 | 0.3848 | 0.7778 | 0.1328 | 0.1357 | 0.0000 | 0.1263 | 0.3770 | 0.1380 |
| heteroskedastic_interactions | 1000 | estimated_outcome_oracle_density | oracle | 24 | 0 | 0.1026 | 0.9352 | 0.3304 | 0.9028 | -0.0006 | 0.0296 | 1.0000 | 0.0787 | 0.2701 | 0.0249 |
| heteroskedastic_interactions | 1000 | oracle_all | oracle | 24 | 0 | 0.0475 | 0.9537 | 0.1450 | 0.9537 | -0.0026 | 0.0174 | 1.0000 | 0.0000 | 0.0000 | 0.0000 |
| heteroskedastic_interactions | 1000 | oracle_outcome_estimated_density | location_scale | 24 | 0 | 0.0425 | 0.9537 | 0.1337 | 0.9769 | -0.0008 | 0.0179 | 1.0000 | 0.0188 | 0.0620 | 0.0084 |
| heteroskedastic_interactions | 1000 | oracle_outcome_estimated_density | residual_kde | 24 | 0 | 0.0430 | 0.9352 | 0.1450 | 0.9398 | -0.0003 | 0.0160 | 1.0000 | 0.0415 | 0.1432 | 0.0098 |
| mixture_kink | 600 | estimated_all | location_scale | 24 | 0 | 0.0946 | 0.9213 | 0.3636 | 0.8611 | -0.0102 | 0.0316 | 1.0000 | 0.0753 | 0.2761 | 0.0215 |
| mixture_kink | 600 | estimated_all | residual_kde | 24 | 0 | 0.0612 | 0.9028 | 0.2506 | 0.7546 | -0.0745 | 0.0826 | 0.5000 | 0.0456 | 0.1606 | 0.0808 |
| mixture_kink | 600 | estimated_outcome_oracle_density | oracle | 24 | 0 | 0.0651 | 0.9583 | 0.3010 | 0.8704 | -0.0044 | 0.0371 | 1.0000 | 0.0430 | 0.1888 | 0.0266 |
| mixture_kink | 600 | oracle_all | oracle | 24 | 0 | 0.0426 | 0.9074 | 0.1879 | 0.8241 | -0.0002 | 0.0260 | 1.0000 | 0.0000 | 0.0000 | 0.0000 |
| mixture_kink | 600 | oracle_outcome_estimated_density | location_scale | 24 | 0 | 0.0602 | 0.9213 | 0.2385 | 0.8704 | 0.0006 | 0.0314 | 1.0000 | 0.0312 | 0.1079 | 0.0184 |
| mixture_kink | 600 | oracle_outcome_estimated_density | residual_kde | 24 | 0 | 0.0405 | 0.9213 | 0.1750 | 0.8241 | 0.0069 | 0.0256 | 0.9583 | 0.0120 | 0.0469 | 0.0144 |
| mixture_kink | 1000 | estimated_all | location_scale | 24 | 0 | 0.0609 | 0.8935 | 0.2399 | 0.7361 | -0.0047 | 0.0349 | 0.9583 | 0.0487 | 0.1494 | 0.0237 |
| mixture_kink | 1000 | estimated_all | residual_kde | 24 | 0 | 0.0470 | 0.8889 | 0.2065 | 0.7269 | -0.0806 | 0.0871 | 0.2500 | 0.0362 | 0.1283 | 0.0800 |
| mixture_kink | 1000 | estimated_outcome_oracle_density | oracle | 24 | 0 | 0.0466 | 0.9213 | 0.1934 | 0.8102 | -0.0076 | 0.0311 | 1.0000 | 0.0256 | 0.0944 | 0.0235 |
| mixture_kink | 1000 | oracle_all | oracle | 24 | 0 | 0.0361 | 0.9213 | 0.1637 | 0.7917 | -0.0049 | 0.0251 | 1.0000 | 0.0000 | 0.0000 | 0.0000 |
| mixture_kink | 1000 | oracle_outcome_estimated_density | location_scale | 24 | 0 | 0.0461 | 0.9491 | 0.2011 | 0.8194 | -0.0034 | 0.0228 | 1.0000 | 0.0207 | 0.0746 | 0.0106 |
| mixture_kink | 1000 | oracle_outcome_estimated_density | residual_kde | 24 | 0 | 0.0327 | 0.9120 | 0.1568 | 0.7963 | -0.0030 | 0.0225 | 1.0000 | 0.0107 | 0.0394 | 0.0103 |

Ablation interpretation:

- `*_vs_oracle_*` columns compare each layer with `oracle_all` on the same
  replication.
- Large `oracle_outcome_estimated_density` errors isolate density and
  density-derivative problems.
- Large `estimated_outcome_oracle_density` errors isolate outcome-regression
  and outcome-derivative problems.
- `dml_flex` curve diagnostics are in level and derivative columns.
- `dml_incremental` diagnostics are in scalar columns.

## Validation Figures

Visualization script:

```text
inst/plot_validation_results.R
```

Generated outputs:

```text
inst/validation_figures/oracle_rmse.png
inst/validation_figures/stress_diagnostics.png
inst/validation_figures/ablation_scalar_rmse.png
inst/validation_figures/oracle_rmse_plot_data.csv
inst/validation_figures/stress_diagnostics_plot_data.csv
inst/validation_figures/ablation_scalar_plot_data.csv
```

The figures summarize:

- oracle RMSE and coverage for level, derivative, and scalar derivative
  targets;
- stress-test sensitivity to density method;
- nuisance-oracle ablation for the scalar derivative.

## Verification

Recorded verification:

- `testthat::test_file("tests/testthat/test-continuous.R")`: 139 passed,
  3 skipped on CRAN.
- RStudio bundled Quarto render:

```text
/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto render
03-continuous.Rmd --to html
```

produced:

```text
tutorial/_book/03-continuous.html
```

- `R CMD build --no-build-vignettes fdid-main`: succeeded.
- `R CMD check --no-manual --no-build-vignettes fdid_1.0.2.tar.gz`:
  2 warnings, 0 errors, 0 notes.
- Earlier local validation recorded:
  - `devtools::test(filter = "continuous")`: 128 passed.
  - `devtools::test()`: 199 passed, with one pre-existing warning about
    `target.pop` for `did`.
  - `R CMD check --no-manual --no-build-vignettes`: no errors; two expected
    warnings from not building vignettes into `inst/doc`.
- 2026-05-19 documentation-audit packaging check:
  `R CMD build --no-build-vignettes fdid-main` succeeded and produced
  `fdid_1.0.2.tar.gz`. A tar listing confirmed that
  `historical_mds_2026-05-19/` and validation result folders were excluded,
  while the three active project markdowns remained included.

## Current Verdict

Current status:

```text
PASS WITH CAVEATS FOR DEVELOPMENT USE
```

The implemented DML signals are coherent with the intended targets:

- `dml_binary`: binary AIPW/IRM signal.
- `dml_plr`: partially linear residualized score.
- `dml_flex`: cross-fitted level-curve pseudo-outcome mapped to `G`.
- `dml_incremental`: separate scalar average-derivative orthogonal signal.

Cross-fitting is used in all DML paths.

Remaining caveats:

- `dml_flex` practical bands are not final theorem-level DML uniform bands.
- Level-comparison and fixed-`g` derivative reporting are first-class through
  `fdid_contrast()` and `fdid_derivative()`. Kernel bootstrap fits support
  replicate-based and curve-covariance-based grid contrasts; theorem-level
  efficient contrast inference remains open for other curve paths.
- Conditional density remains the most important nuisance for continuous-G DML.
- Derivative targets are more fragile than level targets.
- Scalar DML bootstrap and jackknife are unsupported.
- External benchmark validation remains future work.

## Current Important Files

Implementation:

```text
R/fdid.R
R/kernel.R
R/dml.R
R/fdid_prepare.R
R/plot.R
R/summary.R
```

Tests and validation:

```text
tests/testthat/test-continuous.R
inst/oracle_signal_validation.R
inst/oracle_signal_stress_tests.R
inst/nuisance_oracle_ablation.R
inst/plot_validation_results.R
inst/oracle_signal_validation_results/
inst/oracle_signal_stress_results/
inst/nuisance_oracle_ablation_results/
inst/validation_figures/
```

Tutorial:

```text
tutorial/03-continuous.Rmd
tutorial/04-visualization.Rmd
tutorial/05-sensitivity.Rmd
```

Current active docs:

```text
HANDOFF.md
TODO.md
CONTINUOUS_G_PROJECT_LOG.md
CONTINUOUS_G_ESTIMATOR_REFERENCE.md
```

Build/config:

```text
DESCRIPTION
NAMESPACE
.Rbuildignore
```

---

# 2026-06-01 Tutorial Enrichment And Public-Branch Cleanup

## Summary

Implemented `TUTORIAL_ENRICHMENT_PLAN_2026-06-01.md` through the StatsClaw
workflow on the public staging branch
`continuous-g-dev-2026-06-01-tutorial-enrichment` in
`/private/tmp/enhan-fdid-sync`.

Main outcomes:

- Restored the richer Chapter 1 installation/data/background tutorial material.
- Restored the richer Chapter 2 binary-G workflow and added a `dml_binary`
  section explaining it as an AIPW/IRM-style binary contrast with ML nuisance
  learners and cross-fitting.
- Removed the standalone reporting/inference chapter and kept reporting inside
  the estimator and visualization chapters.
- Restored detailed continuous-G kernel, DML, and visualization chapters.
- Renumbered sensitivity analysis to Chapter 7 and kept its content
  substantively unchanged.
- Regenerated the public `docs/` site with rendered outputs and 43 PNG figure
  assets.
- Removed local development/workflow markdown files from the public branch
  root while preserving them in this local checkout. The public branch root now
  keeps only `LICENSE.md` and `NEWS.md`.

## Public Branch Commits

```text
f27250f scriber: enrich fdid tutorial
0bd326d docs: finalize tutorial public cleanup
```

The second commit removes the generated `ARCHITECTURE.md` from the public
branch root to preserve the public/private documentation boundary.

The branch was pushed to:

```text
origin/continuous-g-dev-2026-06-01-tutorial-enrichment
```

## Validation

Validation completed in the staging repo:

```text
R CMD INSTALL -l /private/tmp/fdidlib .
env R_LIBS=/private/tmp/fdidlib Rscript -e 'devtools::test(".")'
env R_LIBS=/private/tmp/fdidlib quarto render tutorial
git diff --check
```

Results:

- Package install passed.
- `devtools::test(".")` passed with 374 passes, 0 failures, and 1 existing
  warning about `target.pop` not changing `did`/`ols1`.
- Quarto rendered all tutorial chapters successfully.
- `git diff --check` passed after mechanical cleanup of Quarto-generated
  trailing whitespace.
- A tester child-agent dispatch was attempted but blocked by a Codex usage
  limit before producing `audit.md`; leader completed the equivalent
  `test-spec.md` validation and recorded it in the StatsClaw run directory.

StatsClaw run artifacts are stored under:

```text
/private/tmp/statsclaw-data/workspace/Statsclaw-FDID/runs/TUTORIAL-ENRICH-20260601-0001
```
