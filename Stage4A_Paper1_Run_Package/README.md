# Stage 4A — Paper 1 GCV refinement run package

**Scope:** Original *processed* Yield, same second-stage `mgcv` GAM and GCV criterion used in Stage 3C. Reviewer preprocessing comparisons are deliberately on hold. Paper 2 is out of scope.

## Inputs

The five accepted Stage 3B/3C result files are already included in `inputs/` and were verified here. **One additional input is required:** the 4,776-row processed `Example1SmoothHarvestPasses.csv`, used in Stages 1–3C. Place it in `inputs/`, or leave it at its original location: `C:/Users/pbbs/Claussen/Example1SmoothHarvestPasses.csv`. Do **not** substitute `unsmoothed_grid.csv`.

## R dependencies and execution

Uses R, the external `mgcv` package, and Base R `stats`, `tools`, and graphics (built in). Install once if needed: `install.packages("mgcv")`.

1. Extract the ZIP, keeping all files and the `inputs/` folder together.
2. In RStudio/R console:

```r
setwd("C:/Users/pbbs/Downloads/Stage4A_Paper1_Run_Package") # change as needed
if (!requireNamespace("mgcv", quietly = TRUE)) install.packages("mgcv")
source("Stage4A_Full_Run.R", echo = FALSE)
```

Or from Windows Command Prompt opened in the extracted folder:

```bat
Rscript Stage4A_Full_Run.R
```

The script first validates all 199 saved 600-point profiles and the 218 detected dip candidates, **then reproduces fixed-sp model scores and Product B estimates at Eastings 214, 216, and 374 from the processed full dataset**. It performs bracketed Brent optimization on every detected dip. It records version and MD5 provenance, and stops if checks fail.

## Expected new outputs

Saved under `Stage4A_Results/` only after the Stage 4A script executes:

- `Stage4A_Refined_Interior_Minima.csv` — 218 refined candidate rows (if every refit succeeds), localized slopes and checks.
- `Stage4A_Full_Field_Refined_Summary.csv` — 199-Eastings score, coefficient, regime comparison.
- `Stage4A_Regime_Transitions.csv` — observed adjacent-Easting candidate classifications.
- `Stage4A_Validation_Gate.csv` — numerical checks and PASS/FAIL.
- `Stage4A_Refinement_Diagnostics.pdf` — criterion/coefficient plots.
- `Stage4A_Input_Manifest.csv`, `Stage4A_Fit_Errors.csv`, `Stage4A_Provenance.txt`.

**Preflight files in this ZIP are NOT Stage 4A refit results.** They only verify the already accepted Stage 3B/3C saved outputs and construct the 218 optimization brackets. The present execution environment has neither R/`mgcv` nor the full processed Yield CSV, so actual refitted Stage 4A minima have not been generated here.

The local refinement searches only brackets around previously detected dips. It does not prove all stationary points were located or establish a global minimizer on (0, infinity).

## Quality-control note

The Stage 4A script adds three explicit fixed-`sp` reproducibility checks relative to the original planned 04A R script. The original is included unchanged as `04A_Original_Reference_Unmodified.R`.
