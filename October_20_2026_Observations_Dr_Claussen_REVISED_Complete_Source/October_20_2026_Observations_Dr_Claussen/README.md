# October 20, 2026 Observations from Dr. Claussen — revised annotated edition

## Main deliverables
- `October_20_2026_Observations_from_Dr_Claussen.pdf`: revised 18-page report, Sections 1–11 only; no embedded appendices or code dump.
- `October_20_2026_Observations_from_Dr_Claussen.tex`: editable LaTeX source.
- `October_20_2026_Observations_from_Dr_Claussen.qmd`: Quarto reproducible results notebook.
- `generate_figures_tables.R`: Base R code invoked by Quarto to regenerate figures and tables from the saved data (not the original GAM fits).
- `figures/` and `tables/`: current charts/tables and LaTeX fragments.
- `source_data/`: original Stage 2–3C saved CSV results used to generate the exhibits.
- `verbatim/`: four source transcripts copied without changes, preserved separately rather than placed after Section 11.

## Recompiling the LaTeX PDF
Within this directory:
`xelatex -interaction=nonstopmode October_20_2026_Observations_from_Dr_Claussen.tex`
Run a second time to update the table of contents.

## Reproducing the charts from Quarto
Use RStudio/Quarto and open the QMD in the root of this directory. The first executable chunk calls:
`source("generate_figures_tables.R")`
The script expects its `source_data/` directory and creates `generated_figures_tables/`. It uses base R; Quarto output formatting uses standard knitr dependencies. The saved GAM results are inputs; model fits are not rerun.

## Revision details
- Kept exactly 11 main document sections.
- Removed appendices and inline figure-generation code from the PDF, leaving code in the separate QMD/R script.
- Expanded explanations for experimental stages, GCV definitions and transformations, Eastings 176 and 374, criterion preference and optimizer anomalies, and file provenance.
- Added 18 red margin-guide notes (11 major-section notes plus seven detailed-topic notes).
- Source CSVs are cited under the corresponding tables and figures.

## Research caveats
Values describe completed Stage 1A–3C analyses. Grid-detected minima are not proved global minimizers. No Reiss–Ogden h-function roots or later-stage simulations are claimed as completed.
