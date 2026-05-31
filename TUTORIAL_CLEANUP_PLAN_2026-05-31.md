# Tutorial Cleanup Plan

Date: 2026-05-31

## Status

Implementation status: completed in the local StatsClaw scriber worktree for
branch `continuous-g-dev-2026-05-31`; reviewer and shipper follow this pass.

The public tutorial now has one canonical source directory:

```text
tutorial/
```

and one rendered GitHub Pages output directory:

```text
docs/
```

The obsolete public tutorial source copies were removed after migration:

```text
tutorial-full-continuous-g/
tutorial-next-2026-05-30/
```

## Final Tutorial Structure

```text
tutorial/01-start.Rmd
tutorial/02-main.Rmd
tutorial/03-continuous-overview.Rmd
tutorial/04-kernel-continuous-g.Rmd
tutorial/05-dml-continuous-g.Rmd
tutorial/06-visualization.Rmd
tutorial/07-reporting-inference.Rmd
tutorial/08-sensitivity.Rmd
tutorial/index.qmd
tutorial/references.bib
tutorial/references.qmd
tutorial/_quarto.yml
```

`tutorial/_quarto.yml` renders directly to `../docs` and no longer references
the deleted chapter names.

## Target-First Rewrite

The continuous-`G` tutorial now starts from the user's target:

- scalar binary/discrete FDID effect;
- scalar continuous-`G` DML PLR slope;
- smooth level curve;
- derivative curve;
- observed-population average derivative;
- fixed-reference level contrast;
- support and overlap diagnostics.

The inference language remains conservative:

- kernel analytical intervals are pointwise intervals from stored standard
  errors;
- kernel bootstrap pointwise intervals use stored bootstrap curve replicates;
- kernel bootstrap bands are finite-grid quantile-envelope bands, not
  studentized sup-`t` bands;
- scalar DML robust intervals come from orthogonal-score variance;
- scalar score-multiplier intervals use multiplier draws over stored influence
  values;
- non-BLP `dml_flex` bands are practical signal-mapping multiplier bands;
- BLP-spline DML is a signal map with finite-grid covariance/Gaussian max-`t`
  bands when available.

## Rendered Site

Rendered with:

```sh
/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto render tutorial
```

Result: pass. Quarto produced:

```text
docs/index.html
docs/01-start.html
docs/02-main.html
docs/03-continuous-overview.html
docs/04-kernel-continuous-g.html
docs/05-dml-continuous-g.html
docs/06-visualization.html
docs/07-reporting-inference.html
docs/08-sensitivity.html
docs/references.html
docs/search.json
```

`docs/.nojekyll` was restored manually after render because the render did not
create it.

## Validation Results

Scriber-run checks plus leader fallback audit after tester dispatch hit the
child Codex usage limit:

| Command | Result | Notes |
|---|---|---|
| `quarto render tutorial` | blocked | `quarto` was not on `PATH`. |
| `/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto render tutorial` | pass | Rendered final tutorial source into `docs/`; R startup emitted locale warnings only. |
| stale docs chapter scan | pass | `docs/` no longer contains pages for deleted chapter names. |
| `git diff --cached --check` | pass | No whitespace errors after trimming generated HTML trailing whitespace. |
| `find . -maxdepth 4 \( -name '.DS_Store' -o -name '.Rhistory' -o -name 'Rplots.pdf' \) -print` | pass | No matching artifacts found. |
| staged scope scan for `R/`, `man/`, `tests/`, `NAMESPACE`, `DESCRIPTION` | pass | No package implementation, help, test, namespace, or description files staged. |
| R source parse | pass | `Rscript` parse check passed; locale warning only. |
| Rd check | pass | `tools::checkRd` passed; locale warning only. |
| `R CMD INSTALL -l /private/tmp/fdidlib .` | pass | Package installed and loaded from temporary and final locations. |
| RStudio-bundled Quarto render with installed package library | pass | Render succeeded; generated whitespace-only HTML differences were restored to the committed docs. |
| `R CMD build --no-build-vignettes --no-manual /private/tmp/enhan-fdid-sync` | pass | Built `/private/tmp/fdid_1.0.2.tar.gz`. |
| `Rscript -e 'devtools::test(".")'` | pass with warning | `374` passed, `0` failed, `1` existing warning about `target.pop` not changing `method = "did"`. |

Detailed fallback audit evidence is recorded in the StatsClaw run artifact
`audit.md`.

## Changed Files

- Rewrote `tutorial/` around the final eight target-first chapters.
- Regenerated `docs/` from the final tutorial source.
- Removed `tutorial-full-continuous-g/` and `tutorial-next-2026-05-30/`.
- Removed stale `tutorial/summary.qmd`.
- Updated `CONTINUOUS_G_PROJECT_LOG.md`, `TODO.md`, and
  `PACKAGE_PLOTTING_LAYER_UPDATE_PLAN_2026-05-31.md`.
- Added/updated architecture and StatsClaw run artifacts.

## Commit And Push Plan

Local scriber commit message:

```text
scriber: rewrite target-first tutorial
```

Push target after reviewer approval:

```text
origin continuous-g-dev-2026-05-31
```

Shipper should record the final commit hash, push result, reviewer verdict,
and any GitHub Pages follow-up. No co-author trailers should be added.
