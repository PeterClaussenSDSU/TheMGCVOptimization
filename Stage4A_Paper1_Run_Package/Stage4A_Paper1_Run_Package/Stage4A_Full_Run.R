# ============================================================================
# STAGE 4A: LOCAL REFINEMENT OF GRID-DETECTED mgcv GCV MINIMA
# Reuses accepted Stage 3B and 3C outputs. Does NOT compute Reiss-Ogden h's.
# Package: mgcv (external); all other data processing is Base R.
# "optimize" requires an anonymous objective function; unlike Stage 2's
# original advisor-mandated manual script, this script uses a callback.
# ============================================================================

# === PART 1: SETTINGS AND INPUT FILES ===
if (!requireNamespace("mgcv", quietly = TRUE)) {
  stop("Install mgcv once: install.packages('mgcv')")
}
library(mgcv)

# Standalone package: run with working directory set to this folder.
# The five accepted Stage 3B/3C CSVs are shipped in ./inputs/.
# The original *processed* Example1SmoothHarvestPasses.csv must be present
# in ./inputs/ or at the established Windows project path.
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
data_dir <- file.path(base_dir, "Claussen/Stage4A_Paper1_Run_Package/inputs")
project_dir <- base_dir
out_dir <- file.path(base_dir, "Stage4A_Results")
full_file <- file.path(data_dir, "Example1SmoothHarvestPasses.csv")
if (!file.exists(full_file)) {
  windows_input <- "C:/Users/pbbs/Claussen/Stage4A_Paper1_Run_Package/inputs/Example1SmoothHarvestPasses.csv"
  if (file.exists(windows_input)) full_file <- windows_input
}
grid_file <- file.path(data_dir, "Stage3C_Full_Grid.csv")
summary_file <- file.path(data_dir, "Stage3C_Auto_vs_Grid.csv")
minima_file <- file.path(data_dir, "Stage3C_Interior_Grid_Dips.csv")
auto_file <- file.path(data_dir, "Stage3B_Automatic_Selections.csv")
linear_file <- file.path(data_dir, "Stage3B_Linear_Limit_Check.csv")
files_needed <- c(full_file, grid_file, summary_file, minima_file,
                  auto_file, linear_file)
if (!all(file.exists(files_needed))) {
  stop("Stage 4A inputs missing: ",
       paste(files_needed[!file.exists(files_needed)], collapse = "; "),
       "\nPlace original Example1SmoothHarvestPasses.csv in inputs/.
",
       "Do NOT substitute unsmoothed_grid.csv.")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# These tolerances concern numerical agreement, not statistical significance.
check_tol <- 1e-5
inferior_tol <- 1e-4
near_tie_tol <- 0.10
optimization_tol <- 1e-9   # optimize() tolerance in log(sp) units

manifest_file <- file.path(out_dir, "Stage4A_Input_Manifest.csv")
manifest <- data.frame(
  Item = c(files_needed, "R_version", "mgcv_version", "log_grid_n", "model"),
  Value = c(unname(tools::md5sum(files_needed)), R.version.string,
            as.character(packageVersion("mgcv")), "600",
            "Yield ~ Product + s(Northing, bs='bs', k=10, m=c(3,2)); gaussian; GCV.Cp")
)
if (file.exists(manifest_file)) {
  old <- read.csv(manifest_file, stringsAsFactors = FALSE)
  stopifnot(identical(as.character(old$Item), as.character(manifest$Item)),
            identical(as.character(old$Value), as.character(manifest$Value)))
} else {
  write.csv(manifest, manifest_file, row.names = FALSE)
}
cat("STAGE 4A: Input manifest checked; mgcv", as.character(packageVersion("mgcv")), "\n")

# === PART 2: READ, ORDER, AND CROSS-CHECK ACCEPTED RESULTS ===
full <- read.csv(full_file, stringsAsFactors = FALSE)
grid <- read.csv(grid_file, stringsAsFactors = FALSE)
old_summary <- read.csv(summary_file, stringsAsFactors = FALSE)
old_minima <- read.csv(minima_file, stringsAsFactors = FALSE)
auto <- read.csv(auto_file, stringsAsFactors = FALSE)
linear <- read.csv(linear_file, stringsAsFactors = FALSE)

eastings <- seq(2, 398, by = 2)
log_grid <- seq(-12, 25, length.out = 600)
stopifnot(nrow(full) == 4776, nrow(grid) == 119400,
          nrow(old_summary) == 199, nrow(auto) == 199,
          nrow(linear) == 199,
          identical(as.numeric(old_summary$Easting), as.numeric(eastings)),
          identical(as.numeric(auto$Easting), as.numeric(eastings)),
          identical(as.numeric(linear$Easting), as.numeric(eastings)),
          all(auto$Converged))

full$Product <- factor(trimws(as.character(full$Product)), levels = c("E", "B"))
stopifnot(!anyNA(full$Product))
stopifnot(nrow(old_minima) == sum(old_summary$Interior_Grid_Dip_Count))

# The grid and candidate files are converted to a consistent sort order.
grid <- grid[order(grid$Easting, grid$LogLambda), , drop = FALSE]
old_minima <- old_minima[order(old_minima$Easting, old_minima$LogLambda), , drop = FALSE]

# 199 complete 600-point profiles must be present on the prescribed grid.
for (ii in seq_along(eastings)) {
  gg <- grid[grid$Easting == eastings[ii], , drop = FALSE]
  stopifnot(nrow(gg) == 600,
            max(abs(gg$LogLambda - log_grid)) < 1e-11,
            abs(min(gg$GCV) - old_summary$Best_Grid_GCV[ii]) < check_tol,
            abs(auto$GCV[ii] - old_summary$Automatic_GCV[ii]) < check_tol,
            abs(linear$Linear_GCV[ii] - old_summary$Exact_Linear_GCV[ii]) < check_tol)
}
cat("INPUT CHECKS: 199 complete grid profiles,", nrow(old_minima),
    "original interior dip candidates.\n")

# Additional validation: recompute three saved Stage 3C fixed-sp fits
# from the original processed Yield. These are crucial provenance checks.
for (check_e in c(214, 216, 374)) {
  zz <- full[full$Easting == check_e, , drop = FALSE]
  zz <- zz[order(zz$PassNo), , drop = FALSE]
  gg <- grid[grid$Easting == check_e, , drop = FALSE]
  jj <- which.min(gg$GCV)
  check_fit <- mgcv::gam(
    Yield ~ Product + s(Northing, bs = "bs", k = 10, m = c(3, 2)),
    data = zz, family = gaussian(), method = "GCV.Cp",
    sp = exp(gg$LogLambda[jj])
  )
  check_score <- nrow(zz) * sum(residuals(check_fit)^2) /
    (nrow(zz) - sum(check_fit$edf))^2
  check_beta <- as.numeric(coef(check_fit)["ProductB"])
  stopifnot(abs(check_score - gg$GCV[jj]) < check_tol,
            abs(check_beta - gg$Product_B[jj]) < check_tol)
  cat("Fixed-sp model reproduction PASS at Easting", check_e, "\n")
}

# === PART 3: LOCALLY REFINE EVERY DETECTED INTERIOR DIP ===
# No attempt is made to search for unobserved dips outside or between grid
# brackets. Each candidate is bracketed by its immediately adjacent grid nodes.
refined_rows <- vector("list", nrow(old_minima))
fit_errors <- data.frame(Easting = numeric(0), Grid_Index = integer(0),
                         Error = character(0), stringsAsFactors = FALSE)

for (rr in seq_len(nrow(old_minima))) {
  original <- old_minima[rr, , drop = FALSE]
  e <- original$Easting
  j <- as.integer(original$Grid_Index)
  gg <- grid[grid$Easting == e, , drop = FALSE]
  stopifnot(j >= 2L, j <= 599L,
            abs(gg$LogLambda[j] - original$LogLambda) < 1e-10,
            abs(gg$GCV[j] - original$GCV) < check_tol)

  z <- full[full$Easting == e, , drop = FALSE]
  z <- z[order(z$PassNo), , drop = FALSE]
  stopifnot(nrow(z) == 24,
            identical(as.numeric(z$PassNo), as.numeric(1:24)))

  lower <- gg$LogLambda[j - 1L]
  upper <- gg$LogLambda[j + 1L]

  # optimize() minimizes in ONE dimension. Its objective is the exact same
  # fixed-sp mgcv GCV used in Stage 3C, evaluated at log(sp) = theta.
  opt <- try(stats::optimize(
    f = function(theta) {
      f <- mgcv::gam(
        Yield ~ Product + s(Northing, bs = "bs", k = 10, m = c(3, 2)),
        data = z, family = gaussian(), method = "GCV.Cp", sp = exp(theta)
      )
      n <- nrow(z)
      n * sum(residuals(f)^2) / (n - sum(f$edf))^2
    },
    interval = c(lower, upper), tol = optimization_tol
  ), silent = TRUE)

  if (inherits(opt, "try-error")) {
    fit_errors <- rbind(fit_errors, data.frame(Easting = e, Grid_Index = j,
      Error = substr(as.character(opt)[1], 1, 300)))
    next
  }

  best_log <- opt$minimum
  best_fit <- try(mgcv::gam(
    Yield ~ Product + s(Northing, bs = "bs", k = 10, m = c(3, 2)),
    data = z, family = gaussian(), method = "GCV.Cp", sp = exp(best_log)
  ), silent = TRUE)
  if (inherits(best_fit, "try-error")) {
    fit_errors <- rbind(fit_errors, data.frame(Easting = e, Grid_Index = j,
      Error = substr(as.character(best_fit)[1], 1, 300)))
    next
  }
  n <- nrow(z)
  best_gcv <- n * sum(residuals(best_fit)^2) / (n - sum(best_fit$edf))^2
  best_edf <- sum(best_fit$edf) - 2
  best_beta <- unname(as.numeric(coef(best_fit)["ProductB"]))

  # Independent local classification: probe immediately to either side.
  h <- min((upper - lower) / 10, 0.005)
  left_log <- max(lower, best_log - h)
  right_log <- min(upper, best_log + h)
  left_fit <- mgcv::gam(
    Yield ~ Product + s(Northing, bs = "bs", k = 10, m = c(3, 2)),
    data = z, family = gaussian(), method = "GCV.Cp", sp = exp(left_log)
  )
  right_fit <- mgcv::gam(
    Yield ~ Product + s(Northing, bs = "bs", k = 10, m = c(3, 2)),
    data = z, family = gaussian(), method = "GCV.Cp", sp = exp(right_log)
  )
  left_gcv <- n * sum(residuals(left_fit)^2) / (n - sum(left_fit$edf))^2
  right_gcv <- n * sum(residuals(right_fit)^2) / (n - sum(right_fit$edf))^2
  slope_left <- (best_gcv - left_gcv) / (best_log - left_log)
  slope_right <- (right_gcv - best_gcv) / (right_log - best_log)
  local_check <- is.finite(best_gcv) && best_gcv <= left_gcv + 1e-5 &&
    best_gcv <= right_gcv + 1e-5 && best_gcv <= gg$GCV[j] + 1e-5 &&
    abs(best_gcv - opt$objective) < check_tol && best_edf > 1.01 - 1e-4

  refined_rows[[rr]] <- data.frame(
    Easting = e, Original_Grid_Index = j,
    Grid_LogSP = original$LogLambda, Grid_GCV = original$GCV,
    Bracket_Lower = lower, Bracket_Upper = upper,
    Refined_LogSP = best_log, Refined_SP = exp(best_log),
    Refined_GCV = best_gcv, Refined_Smooth_EDF = best_edf,
    Refined_Product_B = best_beta,
    Grid_Minus_Refined_GCV = original$GCV - best_gcv,
    Local_Left_Slope = slope_left, Local_Right_Slope = slope_right,
    Local_Minimum_Check = local_check, stringsAsFactors = FALSE
  )
  if (rr %% 30L == 0L || rr == nrow(old_minima)) {
    cat("Refinement:", rr, "/", nrow(old_minima), "candidates checked.\n")
  }
}
write.csv(fit_errors, file.path(out_dir, "Stage4A_Fit_Errors.csv"), row.names = FALSE)
if (length(fit_errors$Easting) > 0L || any(vapply(refined_rows, is.null, logical(1)))) {
  stop("Local refinement failed; inspect Stage4A_Fit_Errors.csv.")
}
refined <- do.call(rbind, refined_rows)
row.names(refined) <- NULL
refined <- refined[order(refined$Easting, refined$Refined_LogSP), , drop = FALSE]
cat("Refinement complete:", nrow(refined), "candidates; local checks pass:",
    sum(refined$Local_Minimum_Check), "/", nrow(refined), "\n")

# === PART 4: COMPARE REFINED INTERIOR MINIMA WITH THE EXACT LINEAR LIMIT ===
# Distinguish: (1) structural preference; (2) automatic optimization gap;
# (3) coefficient sensitivity. The latter two are NOT equivalent.
summary_rows <- vector("list", length(eastings))
for (ii in seq_along(eastings)) {
  e <- eastings[ii]
  rr <- refined[refined$Easting == e, , drop = FALSE]
  old <- old_summary[ii, , drop = FALSE]
  a <- auto[ii, , drop = FALSE]
  l <- linear[ii, , drop = FALSE]

  interior_gcv <- NA_real_
  interior_beta <- NA_real_
  interior_logsp <- NA_real_
  if (nrow(rr) > 0L) {
    jj <- which.min(rr$Refined_GCV)
    interior_gcv <- rr$Refined_GCV[jj]
    interior_beta <- rr$Refined_Product_B[jj]
    interior_logsp <- rr$Refined_LogSP[jj]
  }
  structural_gap <- interior_gcv - l$Linear_GCV
  preferred_class <- if (is.na(interior_gcv)) "linear_no_dip" else
    if (interior_gcv < l$Linear_GCV) "interior" else "linear_with_dip"

  # Best alternative to automatic choice includes exact linear, grid minimum,
  # and refined interior minima. The best evaluated score also includes auto.
  options <- data.frame(
    Type = c("grid", "exact_linear", if (nrow(rr) > 0L) "refined_interior"),
    GCV = c(old$Best_Grid_GCV, l$Linear_GCV,
            if (nrow(rr) > 0L) interior_gcv),
    Product_B = c(old$Best_Grid_Product_B, l$Linear_Product_B,
                  if (nrow(rr) > 0L) interior_beta),
    stringsAsFactors = FALSE
  )
  best_alt <- options[which.min(options$GCV), , drop = FALSE]
  best_eval_gcv <- min(best_alt$GCV, a$GCV)
  best_eval_beta <- if (a$GCV < best_alt$GCV) a$Product_B else best_alt$Product_B

  summary_rows[[ii]] <- data.frame(
    Easting = e, Interior_Dip_Count = nrow(rr),
    Best_Interior_Refined_LogSP = interior_logsp,
    Best_Interior_Refined_GCV = interior_gcv,
    Best_Interior_Refined_Beta = interior_beta,
    Exact_Linear_GCV = l$Linear_GCV,
    Exact_Linear_Beta = l$Linear_Product_B,
    Structural_Gap_Refined_Interior_Minus_Linear = structural_gap,
    Preferred_Regime = preferred_class,
    Near_Score_Tie = is.finite(structural_gap) && abs(structural_gap) < near_tie_tol,
    Automatic_GCV = a$GCV, Automatic_Beta = a$Product_B,
    Best_Grid_GCV = old$Best_Grid_GCV, Best_Grid_Beta = old$Best_Grid_Product_B,
    Best_Alternative_Source = best_alt$Type,
    Best_Alternative_GCV = best_alt$GCV,
    Best_Alternative_Beta = best_alt$Product_B,
    Auto_Minus_Best_Alternative_GCV = a$GCV - best_alt$GCV,
    Auto_Excess_GCV = max(0, a$GCV - best_alt$GCV),
    Auto_Inferior = (a$GCV - best_alt$GCV) > inferior_tol,
    Auto_Minus_Best_Alternative_Beta = a$Product_B - best_alt$Product_B,
    Best_Evaluated_GCV = best_eval_gcv,
    Best_Evaluated_Beta = best_eval_beta,
    stringsAsFactors = FALSE
  )
}
field <- do.call(rbind, summary_rows)
row.names(field) <- NULL
cat("Refined GCV: inferior automatic selections:", sum(field$Auto_Inferior), "\n")
cat("Regime counts:\n")
print(table(field$Preferred_Regime))
cat("Near-ties (GCV difference <", near_tie_tol, "):",
    sum(field$Near_Score_Tie), "\n")

# === PART 5: RECORD REGIME TRANSITIONS WITHOUT CLAIMING A CONTINUOUS CROSSING ===
transitions <- data.frame()
for (ii in 2:nrow(field)) {
  a <- field[ii - 1L, , drop = FALSE]
  b <- field[ii, , drop = FALSE]
  class_a <- if (a$Preferred_Regime == "interior") "interior" else "linear"
  class_b <- if (b$Preferred_Regime == "interior") "interior" else "linear"
  if (class_a != class_b) {
    transitions <- rbind(transitions, data.frame(
      Easting_From = a$Easting, Easting_To = b$Easting,
      Class_From = class_a, Class_To = class_b,
      Structural_Gap_From = a$Structural_Gap_Refined_Interior_Minus_Linear,
      Structural_Gap_To = b$Structural_Gap_Refined_Interior_Minus_Linear,
      Auto_Beta_Change = b$Automatic_Beta - a$Automatic_Beta,
      Best_Evaluated_Beta_Change = b$Best_Evaluated_Beta - a$Best_Evaluated_Beta,
      stringsAsFactors = FALSE
    ))
  }
}
cat("Observed adjacent preferred-regime changes:", nrow(transitions), "\n")
print(transitions, row.names = FALSE, digits = 8)

# === PART 6: VALIDATION GATE AND SAVE ACCEPTED TABLES ===
checks <- data.frame(
  Check = c("199 full-field input profiles", "One refinement per original dip",
            "All candidate evaluations finite", "Refined GCV not above grid dip",
            "All local neighborhood checks", "199 summary rows",
            "Best-evaluated <= automatic", "No fit errors", "Input manifest recorded"),
  PASS = c(nrow(old_summary) == 199,
           nrow(refined) == nrow(old_minima),
           all(is.finite(as.matrix(refined[, c("Refined_GCV", "Refined_Product_B",
                                              "Refined_LogSP")]))),
           all(refined$Grid_Minus_Refined_GCV > -check_tol),
           all(refined$Local_Minimum_Check),
           nrow(field) == 199,
           all(field$Best_Evaluated_GCV <= field$Automatic_GCV + check_tol),
           nrow(fit_errors) == 0,
           file.exists(manifest_file))
)
print(checks, row.names = FALSE)
write.csv(checks, file.path(out_dir, "Stage4A_Validation_Gate.csv"), row.names = FALSE)
if (!all(checks$PASS)) {
  write.csv(refined, file.path(out_dir, "FAILED_Stage4A_Refined_Candidates.csv"),
            row.names = FALSE)
  stop("Stage 4A validation FAILED. Inspect candidate checks before acceptance.")
}
write.csv(refined, file.path(out_dir, "Stage4A_Refined_Interior_Minima.csv"),
          row.names = FALSE)
write.csv(field, file.path(out_dir, "Stage4A_Full_Field_Refined_Summary.csv"),
          row.names = FALSE)
write.csv(transitions, file.path(out_dir, "Stage4A_Regime_Transitions.csv"),
          row.names = FALSE)
cat("STAGE 4A VALIDATION: PASS. Tables saved.\n")

# === PART 7: DIAGNOSTIC FIGURES (BASE GRAPHICS) ===
pdf(file.path(out_dir, "Stage4A_Refinement_Diagnostics.pdf"),
    width = 11, height = 8.5, onefile = TRUE)
par(mfrow = c(1, 1), mar = c(5, 5, 3, 2))
plot(field$Easting, field$Auto_Excess_GCV, type = "h", lwd = 1.7,
     xlab = "Easting", ylab = "GCV(auto) - GCV(best evaluated alternative), clipped at zero",
     main = "Stage 4A: criterion gap after local refinement")
abline(h = inferior_tol, lty = 2, col = "gray40")
plot(field$Easting, field$Auto_Minus_Best_Alternative_Beta,
     type = "h", xlab = "Easting", ylab = "Product B: auto minus best alternative",
     main = "Stage 4A: coefficient sensitivity to selection")
abline(h = 0, lty = 2)
plot(field$Easting, field$Structural_Gap_Refined_Interior_Minus_Linear,
     type = "p", pch = 16, cex = 0.6,
     xlab = "Easting", ylab = "Best refined interior GCV - exact linear GCV",
     main = "Structural preference diagnostic; NA = no detected dip")
abline(h = 0, lty = 2)
plot(field$Easting, field$Best_Evaluated_Beta, type = "l", lwd = 1.8,
     xlab = "Easting", ylab = "Product B coefficient",
     main = "Automatic, best evaluated, and exact linear coefficient")
lines(field$Easting, field$Automatic_Beta, lty = 2, lwd = 1.6, col = "firebrick")
lines(field$Easting, field$Exact_Linear_Beta, lty = 3, col = "steelblue")
legend("topright", c("Best evaluated", "Automatic", "Exact linear"),
       lty = c(1, 2, 3), col = c("black", "firebrick", "steelblue"), bty = "n")

# Zoom on representative local minima. All curves below are actual Stage 3C
# grid values, and symbols indicate newly reevaluated Stage 4A candidates.
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (ee in c(48, 174, 216, 374)) {
  gg <- grid[grid$Easting == ee, , drop = FALSE]
  rr <- refined[refined$Easting == ee, , drop = FALSE]
  plot(gg$LogLambda, gg$GCV, type = "l", lwd = 1.4,
       xlab = "log(mgcv sp)", ylab = "GCV", main = paste("Easting", ee))
  if (nrow(rr) > 0L) {
    points(rr$Refined_LogSP, rr$Refined_GCV, pch = 19, col = "darkgreen")
  }
  aa <- auto[auto$Easting == ee, , drop = FALSE]
  if (aa$LogLambda >= min(log_grid) && aa$LogLambda <= max(log_grid)) {
    points(aa$LogLambda, aa$GCV, pch = 17, col = "firebrick")
  }
  abline(h = linear$Linear_GCV[match(ee, linear$Easting)], lty = 3)
}
par(mfrow = c(1, 1))
dev.off()

# === PART 8: PROVENANCE AND HANDOFF ===
writeLines(c(
  paste("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("R:", R.version.string),
  paste("mgcv:", as.character(packageVersion("mgcv"))),
  "Stage 4A refines ONLY previously grid-detected interior dips.",
  "Optimization is Brent interval minimization in log(sp) over adjacent grid nodes.",
  "It does NOT establish every mathematical stationary point or a global minimum.",
  "It does NOT compute Reiss-Ogden h_1/h_2/h_3 or test a statistical null.",
  "A regime change across two observed Eastings is not a proven continuous crossing.",
  "Basis/model match Stage 3C: bs='bs', k=10, m=c(3,2); gaussian; GCV.Cp.",
  "Reviewer sensitivity analyses deliberately deferred until after Paper 1 Stage 4A.",
  "Stage 4A validation: PASS"
), file.path(out_dir, "Stage4A_Provenance.txt"))
cat("STAGE 4A COMPLETE. Output:", out_dir, "\n")
cat("NEXT: Review Stage 4A outputs, then revise Paper 1; no Paper 2 work.\n")
