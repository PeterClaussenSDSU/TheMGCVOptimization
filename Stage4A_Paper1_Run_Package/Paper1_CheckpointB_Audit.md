# Paper 1 — Checkpoint B core-manuscript audit

Date: 2026-11-08. This checkpoint follows completion of Stage 4A and **precedes** consideration of the reviewer's preprocessing/smoothing-ladder analysis. No reviewer-sensitivity findings are incorporated.

## Audited and corrected

1. **Reference sets and contrasts:** added `tables/checkpointB_reference_sets.tex`. The 600-point Stage 3C *grid best*, the Stage 4A *best alternative excluding automatic*, and the Stage 4A *best evaluated including automatic* are now explicitly distinct. Signed GCV gaps and signed Product B coefficient differences use the same reference on both sides. The threshold `>1e-4` is explicitly numerical, not inferential.
2. **More precise coefficient interpretation:** recomputed 199-Easting medians from the canonical Stage 3C and 4A tables. Median absolute auto–grid Product B gap = `0.0021980971`; auto–best-**alternative** gap = `7.900919e-9`; auto–best-**evaluated** = `5.096700e-9`; auto–exact-**linear** = `1.1849907249`. This reveals that most tiny auto–grid differences reflect a discretized comparison and do **not** demonstrate inferior optimization. The 16 Stage 3C and Stage 4A inferior selections are unchanged.
3. **Design algebra and replication:** used 24 pass records in `source_data/Example1SmoothHarvestStrip214.csv`. Product B indicator is 12 of 24; correlation with linear Northing is `0.144403`. Intercept-residualized product design quantity is `6`; intercept + linear Northing residual design is `5.874886358670` (97.9% retained). Independently reassembled Stage 2's cubic B-spline basis and integrated curvature penalty, verifying the saved Easting-214 GCV curve to max difference `4.6925e-9`, product coefficient to `1.4558e-9`, and partialled regression identity to `6.5540e-10`. Stage 2 manual-basis `p^T M_lambda p` varies between `5.1178610659` and `5.8748859738` over log lambda `[-1,23]`. Full CSV and source code are bundled. These are design quantities, **not** standard errors nor mgcv-basis output. A separately verified `mgcv` calculation requires an R environment; a guarded R script is included but has not been executed here.
4. **Novelty and references:** checked the scope against Reiss–Ogden (2009), Li–Cao (2023), the `mgcv` model-selection documentation, Green et al. (1985), and other cited literature. The manuscript explicitly avoids claims of discovering multiple GCV minima, inventing grid search, or proving global minima. Revised a multiple-minima citation to one documenting the phenomenon directly and attributed `mgcv` documentation to Simon N. Wood. All existing cited keys are defined.
5. **Figures/tables, labels and provenance:** automated cross-reference audit found no duplicated labels or unresolved manuscript reference targets; the figure filenames exist. Stage 3C captions remain explicitly *grid comparisons*, while Stage 4A captions remain *refined-alternative comparisons*. Some figures are descriptively overlapping (e.g., F11 and F07); the original exhibits have been **preserved** pending target-journal length constraints, rather than removed. The original 4,776-row processed Yield input remains absent from this archive; Stage 4A raw refitting requires it. Manuscript explicitly limits what the completed outputs prove.

## Decisions and work NOT performed

- **Second-stage k = 8/10/12 sensitivity:** not a prerequisite for the paper's deliberately conditional fixed-specification `k=10` numerical case study; it would be a valuable supplementary robustness check if requested during review, but no new k-results have been fabricated. The existing cross-implementation example at Easting 176 is not a substitute for a same-basis k sweep.
- **Full mgcv treatment-design check:** the original fixed GAM design can be diagnosed using the included `scripts/checkpointB_mgcv_design_diagnostic.R`; it includes a strict reproduction gate. This code is unexecuted here because R/mgcv are unavailable, and no mgcv-basis design numbers are asserted in the manuscript.
- **Author details and target journal:** both remain placeholders; do not infer names, affiliation or journal. Final author confirmation, journal formatting, and a journal-specific page/figure reduction require user input.
- **Reviewer work:** separate Checkpoint C remains deferred.
- **Sections 6 and 7:** intentional placeholders until Checkpoint C is assessed.

## Source checks and citation verification

- Reiss & Ogden (2009), `https://doi.org/10.1111/j.1467-9868.2008.00695.x`.
- Li & Cao (2023), `https://doi.org/10.1007/s11222-022-10178-z`.
- Wood, `mgcv::gam.selection` documentation, `https://www.stat.ethz.ch/R-manual/R-devel/library/mgcv/html/gam.selection.html`.
- Green, Jennison & Seheult (1985), `https://doi.org/10.1111/j.2517-6161.1985.tb01358.x`.
- Claussen (2025), South Dakota State University dissertation, `https://openprairie.sdstate.edu/etd2/1720/`.
