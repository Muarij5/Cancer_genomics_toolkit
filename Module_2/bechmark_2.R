# ============================================================
# BENCHMARK SCRIPT — Module 2 (Bliss & Loewe Synergy Calculator)
# Validates your tool against the `synergyfinder` package's own
# built-in reference dataset (mathews_screening_data).
# No external download needed — the reference data ships inside
# the package itself.
# ============================================================

# ------------------------------------------------------------
# 0. CONFIG — edit this one path for your setup
# ------------------------------------------------------------
TOOL_PATH     <- "D:/Tool_genomic/module_1/finalized_glm_version.R"   # path to your Shiny app script
DOSE_UNIT     <- "uM"                        # mathews_screening_data doses are in nM
RESPONSE_TYPE <- "inhibition"                # mathews_screening_data is %viability

# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("synergyfinder", quietly = TRUE)) BiocManager::install("synergyfinder")

required_pkgs <- c("shiny", "drc", "synergyfinder")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  suppressMessages(library(p, character.only = TRUE))
}

# ------------------------------------------------------------
# 2. Load ONLY the calculation functions from your tool
#    (stops before the Shiny UI/server code so this runs headless,
#     with no browser/app launch)
# ------------------------------------------------------------
tool_lines <- readLines(TOOL_PATH, warn = FALSE)
ui_marker  <- grep("^ui <- fluidPage", tool_lines)[1]
if (is.na(ui_marker)) stop("Could not find the UI section marker in TOOL_PATH — check the file.")

core_lines <- tool_lines[1:(ui_marker - 1)]

# Strip out any setwd(...) calls in the source file — a common cause of
# "cannot change working directory" errors when sourced from elsewhere.
setwd_hits <- grep("^\\s*setwd\\(", core_lines)
if (length(setwd_hits) > 0) {
  message("Skipping ", length(setwd_hits), " setwd() call(s) found in TOOL_PATH (not needed here).")
  core_lines <- core_lines[-setwd_hits]
}

core_file <- tempfile(fileext = ".R")
writeLines(core_lines, core_file)
source(core_file)

message("Loaded calculation functions from: ", TOOL_PATH)
# Shiny's validate()/need() only work inside a live Shiny app. Since we're
# calling analyze_single_pair() directly from a plain script, replace them
# with plain-R equivalents that behave the same way (stop with message if
# the check fails) but don't require a Shiny reactive context.
validate <- function(...) {
  results <- list(...)
  for (r in results) if (!is.null(r)) stop(r, call. = FALSE)
}
need <- function(expr, message = "") {
  if (isTRUE(expr)) return(NULL)
  message
}
# ------------------------------------------------------------
# 3. Load the built-in reference dataset and compute ground-truth
#    Bliss/Loewe scores using the reference package itself
# ------------------------------------------------------------
data("ONEIL_screening_data", package = "synergyfinder")
sf_raw <- ONEIL_screening_data

pick_col <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) {
    stop("None of these columns found: ", paste(candidates, collapse = ", "),
         "\nActual columns are: ", paste(colnames(df), collapse = ", "))
  }
  hit[1]
}

col_block <- pick_col(sf_raw, c("block_id", "BlockId", "PairIndex"))
col_d1    <- pick_col(sf_raw, c("drug_row", "DrugRow", "Drug1", "drug1"))
col_d2    <- pick_col(sf_raw, c("drug_col", "DrugCol", "Drug2", "drug2"))
col_c1    <- pick_col(sf_raw, c("conc_r", "ConcRow", "Conc1", "conc1"))
col_c2    <- pick_col(sf_raw, c("conc_c", "ConcCol", "Conc2", "conc2"))
col_resp  <- pick_col(sf_raw, c("response", "Response", "inhibition", "Inhibition"))
col_cell_candidates <- intersect(c("cell_line_name", "CellLine", "cell_line"), colnames(sf_raw))
col_cell  <- if (length(col_cell_candidates) > 0) col_cell_candidates[1] else NA

std <- data.frame(
  Pair_ID   = as.character(sf_raw[[col_block]]),  # named Pair_ID so YOUR tool's
  # detect_pairs() preserves this exact
  # grouping instead of rebuilding its own
  # Drug_A/Drug_B/Cell_Line-based key
  Drug_A    = sf_raw[[col_d1]],
  Drug_B    = sf_raw[[col_d2]],
  Cell_Line = if (!is.na(col_cell)) sf_raw[[col_cell]] else "TMD8",
  Dose_A    = sf_raw[[col_c1]],
  Dose_B    = sf_raw[[col_c2]],
  Response  = sf_raw[[col_resp]],
  stringsAsFactors = FALSE
)

message("Loaded ", nrow(std), " rows from mathews_screening_data across ",
        length(unique(std$Pair_ID)), " blocks.")

# Compute reference Bliss & Loewe scores with the reference package
reshaped      <- ReshapeData(sf_raw, data_type = RESPONSE_TYPE)
ref_bliss_obj <- CalculateSynergy(reshaped, method = "Bliss")
ref_loewe_obj <- CalculateSynergy(reshaped, method = "Loewe")

# Extract per-dose-pair reference scores into a long data frame.
# synergyfinder's internal object layout has changed across versions,
# so we try a few known layouts and fall back to a diagnostic dump
# (rather than a bare crash) if none match.
extract_ref_scores <- function(obj, method_label) {
  cand_tbls <- list(obj$synergy_scores, obj$scores, obj$drug_pair)
  for (tbl in cand_tbls) {
    if (is.data.frame(tbl)) {
      score_col <- grep(paste0("^", method_label, "_synergy$|^", method_label, "$"),
                        colnames(tbl), value = TRUE, ignore.case = TRUE)
      conc1_col <- grep("^conc1$|^conc_r$", colnames(tbl), value = TRUE, ignore.case = TRUE)
      conc2_col <- grep("^conc2$|^conc_c$", colnames(tbl), value = TRUE, ignore.case = TRUE)
      block_col <- grep("^block_id$", colnames(tbl), value = TRUE, ignore.case = TRUE)
      if (length(score_col) && length(conc1_col) && length(conc2_col) && length(block_col)) {
        return(data.frame(
          Pair_Key = as.character(tbl[[block_col[1]]]),
          Dose_A   = tbl[[conc1_col[1]]],
          Dose_B   = tbl[[conc2_col[1]]],
          RefScore = tbl[[score_col[1]]]
        ))
      }
    }
  }
  message("Could not auto-locate ", method_label, " scores. Structure of the object:")
  str(obj, max.level = 2)
  stop("extract_ref_scores() needs a small tweak to match the structure printed above — ",
       "paste that printed structure back and I'll fix the field names.")
}

ref_bliss_long <- extract_ref_scores(ref_bliss_obj, "Bliss")
ref_loewe_long <- extract_ref_scores(ref_loewe_obj, "Loewe")
names(ref_bliss_long)[names(ref_bliss_long) == "RefScore"] <- "Ref_Bliss"
names(ref_loewe_long)[names(ref_loewe_long) == "RefScore"] <- "Ref_Loewe"

# IMPORTANT: your tool's analyze_single_pair() converts doses to uM internally
# via convert_to_uM(df$Dose_A, unit) before storing them in its output --
# so its returned Dose_A/Dose_B are always in uM regardless of DOSE_UNIT.
# The reference table above is still in mathews_screening_data's native units
# (nM), so convert it to uM here too, to match what your tool will report.
uM_conversion_factor <- switch(DOSE_UNIT,
                               "nM"    = 1/1000,
                               "uM"    = 1,
                               "mM"    = 1000,
                               "mg/mL" = 1000,
                               1
)
ref_bliss_long$Dose_A <- ref_bliss_long$Dose_A * uM_conversion_factor
ref_bliss_long$Dose_B <- ref_bliss_long$Dose_B * uM_conversion_factor
ref_loewe_long$Dose_A <- ref_loewe_long$Dose_A * uM_conversion_factor
ref_loewe_long$Dose_B <- ref_loewe_long$Dose_B * uM_conversion_factor

# ------------------------------------------------------------
# 4. Run YOUR tool on every drug-pair block
# ------------------------------------------------------------
std <- detect_pairs(std)
pair_keys <- unique(std$Pair_Key)

message("Running your tool on ", length(pair_keys), " blocks...")

all_details <- list()
n_failed <- 0

for (pk in pair_keys) {
  sub_df <- std[std$Pair_Key == pk, ]
  res <- tryCatch(
    analyze_single_pair(sub_df, unit = DOSE_UNIT, response_type = RESPONSE_TYPE),
    error = function(e) {
      message("  Skipped block '", pk, "': ", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(res)) {
    res$detail$Pair_Key <- pk
    # Add curve-fitting diagnostics so we can tell if a Bliss/Loewe mismatch
    # comes from single-agent curve fitting vs. the synergy formula itself.
    # (inh_A_alone / inh_B_alone are already inside res$detail from analyze_single_pair)
    res$detail$Curve_A_Used_LL4 <- isTRUE(res$curve_A_fit)
    res$detail$Curve_B_Used_LL4 <- isTRUE(res$curve_B_fit)
    res$detail$EC50_A <- res$ec50_A
    res$detail$EC50_B <- res$ec50_B
    all_details[[as.character(pk)]] <- res$detail
  } else {
    n_failed <- n_failed + 1
  }
}

if (length(all_details) == 0) stop("Your tool produced no usable results on this dataset — check analyze_single_pair().")

merged <- do.call(rbind, all_details)
message(n_failed, " of ", length(pair_keys), " blocks could not be analyzed (see messages above).")

# ------------------------------------------------------------
# 5. Attach reference scores by matching on Pair_Key + doses
# ------------------------------------------------------------
# Round doses on both sides before merging — exact floating-point equality
# often fails even for "the same" dose value once it's passed through
# reshaping/unit-conversion steps on either side.
round_doses <- function(df, digits = 6) {
  df$Dose_A <- signif(as.numeric(df$Dose_A), digits)
  df$Dose_B <- signif(as.numeric(df$Dose_B), digits)
  df
}
merged         <- round_doses(merged)
ref_bliss_long <- round_doses(ref_bliss_long)
ref_loewe_long <- round_doses(ref_loewe_long)

merged <- merge(merged, ref_bliss_long, by = c("Pair_Key", "Dose_A", "Dose_B"), all.x = TRUE)
merged <- merge(merged, ref_loewe_long, by = c("Pair_Key", "Dose_A", "Dose_B"), all.x = TRUE)

n_matched <- sum(!is.na(merged$Ref_Bliss))
message(n_matched, " of ", nrow(merged), " dose pairs matched to a reference score.")

if (n_matched == 0) {
  # Don't just die here — show exactly what doesn't line up so the actual
  # cause (rounding vs. a real unit scale mismatch, e.g. nM vs uM) is visible.
  cat("\n--- DIAGNOSTIC: dose values from YOUR tool (merged) ---\n")
  print(head(sort(unique(merged$Dose_A)), 10))
  print(head(sort(unique(merged$Dose_B)), 10))
  cat("\n--- DIAGNOSTIC: dose values from the REFERENCE table ---\n")
  print(head(sort(unique(ref_bliss_long$Dose_A)), 10))
  print(head(sort(unique(ref_bliss_long$Dose_B)), 10))
  cat("\n--- DIAGNOSTIC: Pair_Key values on each side ---\n")
  cat("Your tool's Pair_Key values: ", paste(head(unique(merged$Pair_Key), 10), collapse=", "), "\n")
  cat("Reference Pair_Key values:   ", paste(head(unique(ref_bliss_long$Pair_Key), 10), collapse=", "), "\n")
  stop("No dose pairs matched. Compare the diagnostic values printed above:\n",
       "- If doses differ by a constant factor (e.g. 1000x), your tool is likely ",
       "converting units (nM<->uM) internally — check DOSE_UNIT / any conversion code in your tool.\n",
       "- If Pair_Key values look completely different in format, detect_pairs() built a ",
       "different key than expected — paste both lists back and I will fix the join.")
}

# Derive reference classifications using YOUR tool's own threshold,
# so the comparison is apples-to-apples (same cutoff on both sides)
merged$Ref_Bliss_Class <- classify_score(merged$Ref_Bliss)
merged$Ref_Loewe_Class <- classify_score(merged$Ref_Loewe)

# ------------------------------------------------------------
# 6. Continuous-score agreement: correlation & error metrics
# ------------------------------------------------------------
score_agreement <- function(predicted, reference, label) {
  ok <- is.finite(predicted) & is.finite(reference)
  n_ok <- sum(ok)
  cat("\n---", label, "score agreement (n =", n_ok, ") ---\n")
  if (n_ok < 3) {
    cat("Not enough overlapping points to compute correlation.\n")
    return(invisible(NULL))
  }
  pearson  <- cor(predicted[ok], reference[ok], method = "pearson")
  spearman <- cor(predicted[ok], reference[ok], method = "spearman")
  mae  <- mean(abs(predicted[ok] - reference[ok]))
  rmse <- sqrt(mean((predicted[ok] - reference[ok])^2))
  cat(sprintf("Pearson r   : %.3f\n", pearson))
  cat(sprintf("Spearman rho: %.3f\n", spearman))
  cat(sprintf("MAE         : %.3f\n", mae))
  cat(sprintf("RMSE        : %.3f\n", rmse))
  invisible(list(pearson = pearson, spearman = spearman, mae = mae, rmse = rmse, n = n_ok))
}

bliss_agreement <- score_agreement(merged$Bliss_Score, merged$Ref_Bliss, "BLISS")
loewe_agreement <- score_agreement(merged$Loewe_Score, merged$Ref_Loewe, "LOEWE")

# ------------------------------------------------------------
# 7. Classification agreement: precision, recall, F1, accuracy
# ------------------------------------------------------------
compute_classification_metrics <- function(predicted, actual,
                                           classes = c("Synergy", "Antagonism", "Additive"),
                                           label = "") {
  ok <- !is.na(predicted) & !is.na(actual)
  predicted <- factor(predicted[ok], levels = classes)
  actual    <- factor(actual[ok],    levels = classes)
  
  cm <- table(Predicted = predicted, Actual = actual)
  
  metrics <- data.frame(Class = classes, Precision = NA_real_, Recall = NA_real_,
                        F1 = NA_real_, Support = NA_integer_)
  for (i in seq_along(classes)) {
    cl <- classes[i]
    tp <- cm[cl, cl]
    fp <- sum(cm[cl, ]) - tp
    fn <- sum(cm[, cl]) - tp
    precision <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
    recall    <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
    f1 <- if (is.na(precision) || is.na(recall) || (precision + recall) == 0) {
      NA_real_
    } else {
      2 * precision * recall / (precision + recall)
    }
    metrics[i, c("Precision", "Recall", "F1", "Support")] <-
      c(precision, recall, f1, sum(actual == cl, na.rm = TRUE))
  }
  
  accuracy <- sum(diag(cm)) / sum(cm)
  macro_f1 <- mean(metrics$F1, na.rm = TRUE)
  
  cat("\n---", label, "classification metrics ---\n")
  cat("Confusion matrix (rows = predicted, cols = actual):\n")
  print(cm)
  cat("\nPer-class metrics:\n")
  print(metrics, row.names = FALSE)
  cat(sprintf("\nOverall accuracy: %.3f\n", accuracy))
  cat(sprintf("Macro-average F1: %.3f\n", macro_f1))
  
  invisible(list(confusion_matrix = cm, per_class = metrics, accuracy = accuracy, macro_f1 = macro_f1))
}

bliss_class_metrics <- compute_classification_metrics(merged$Bliss_Class, merged$Ref_Bliss_Class, label = "BLISS")
loewe_class_metrics <- compute_classification_metrics(merged$Loewe_Class, merged$Ref_Loewe_Class, label = "LOEWE")

# ------------------------------------------------------------
# 8. Save the full row-by-row comparison for manual review
# ------------------------------------------------------------
out_cols <- c("Pair_Key", "Drug_A", "Drug_B", "Cell_Line", "Dose_A", "Dose_B",
              "inh_A_alone", "inh_B_alone", "Curve_A_Used_LL4", "Curve_B_Used_LL4",
              "EC50_A", "EC50_B",
              "Bliss_Score", "Ref_Bliss", "Bliss_Class", "Ref_Bliss_Class",
              "Loewe_Score", "Ref_Loewe", "Loewe_Class", "Ref_Loewe_Class")
out_cols <- intersect(out_cols, colnames(merged))
write.csv(merged[, out_cols], "benchmark_comparison_results2.csv", row.names = FALSE)

cat("\nFull comparison table written to: benchmark_comparison_results.csv\n")

# ------------------------------------------------------------
# 9. Diagnostic: raw single-agent inhibition points (pre-curve-fit)
#    This shows the actual data going INTO the curve fit, independent
#    of any fitting logic -- so we can tell whether a flat/odd curve
#    is a data characteristic or a fitting artifact.
# ------------------------------------------------------------
raw_diag <- list()
for (pk in pair_keys) {
  sub_df <- std[std$Pair_Key == pk, ]
  sub_df$Dose_A_uM <- convert_to_uM(sub_df$Dose_A, DOSE_UNIT)
  sub_df$Dose_B_uM <- convert_to_uM(sub_df$Dose_B, DOSE_UNIT)
  control_rows <- sub_df[sub_df$Dose_A == 0 & sub_df$Dose_B == 0, ]
  control_val <- if (nrow(control_rows) > 0) mean(control_rows$Response, na.rm = TRUE) else 100
  sub_df$inhibition <- 1 - (sub_df$Response / control_val)
  
  a_alone <- sub_df[sub_df$Dose_B_uM == 0 & sub_df$Dose_A_uM > 0, c("Dose_A_uM","inhibition")]
  if (nrow(a_alone) > 0) {
    a_alone$Pair_Key <- pk; a_alone$Which_Drug <- "A"
    names(a_alone)[1] <- "Dose_uM"
    raw_diag[[paste0(pk, "_A")]] <- a_alone
  }
  b_alone <- sub_df[sub_df$Dose_A_uM == 0 & sub_df$Dose_B_uM > 0, c("Dose_B_uM","inhibition")]
  if (nrow(b_alone) > 0) {
    b_alone$Pair_Key <- pk; b_alone$Which_Drug <- "B"
    names(b_alone)[1] <- "Dose_uM"
    raw_diag[[paste0(pk, "_B")]] <- b_alone
  }
}
raw_diag_df <- do.call(rbind, raw_diag)
raw_diag_df <- raw_diag_df[order(raw_diag_df$Pair_Key, raw_diag_df$Which_Drug, raw_diag_df$Dose_uM), ]
write.csv(raw_diag_df, "D:/Tool_genomic/module_1/raw_single_agent_points2.csv", row.names = FALSE)
cat("Raw single-agent data (pre-fit) written to: raw_single_agent_points.csv\n")

# ------------------------------------------------------------
# 10. Direct curve comparison: fit synergyfinder's OWN single-agent
#     curve (via its exported FitDoseResponse/PredictResponse functions)
#     on the exact same raw points, and predict at the exact same doses
#     your tool used. This isolates curve-fitting differences from
#     everything else (Bliss/Loewe formula, thresholds, etc.).
# ------------------------------------------------------------
fit_and_predict_ref_curve <- function(dose_vals, response_vals, predict_at) {
  df <- data.frame(dose = dose_vals, response = response_vals)
  model <- tryCatch(
    FitDoseResponse(df, Emin = NA, Emax = NA),
    error = function(e) { message("FitDoseResponse() failed: ", conditionMessage(e)); NULL }
  )
  if (is.null(model)) return(rep(NA_real_, length(predict_at)))
  
  # Confirmed signature: PredictResponse(df, dose) -- it takes the RAW
  # data.frame (not the fitted model object) and re-fits internally.
  pred <- tryCatch(
    sapply(predict_at, function(d) PredictResponse(df, dose = d)),
    error = function(e) {
      message("PredictResponse(df, dose) still failed. Signature is:")
      print(args(PredictResponse))
      message("Original error: ", conditionMessage(e))
      rep(NA_real_, length(predict_at))
    }
  )
  as.numeric(pred)
}

ref_curve_compare <- list()
for (pk in pair_keys) {
  a_data <- raw_diag_df[raw_diag_df$Pair_Key == pk & raw_diag_df$Which_Drug == "A", ]
  b_data <- raw_diag_df[raw_diag_df$Pair_Key == pk & raw_diag_df$Which_Drug == "B", ]
  block_rows <- merged[merged$Pair_Key == pk, ]
  if (nrow(a_data) < 3 || nrow(b_data) < 3 || nrow(block_rows) == 0) next
  
  ref_inh_A <- fit_and_predict_ref_curve(a_data$Dose_uM, a_data$inhibition, block_rows$Dose_A)
  ref_inh_B <- fit_and_predict_ref_curve(b_data$Dose_uM, b_data$inhibition, block_rows$Dose_B)
  
  ref_curve_compare[[pk]] <- data.frame(
    Pair_Key = pk,
    Dose_A = block_rows$Dose_A, Dose_B = block_rows$Dose_B,
    inh_A_alone = block_rows$inh_A_alone, ref_inh_A_alone = ref_inh_A,
    inh_B_alone = block_rows$inh_B_alone, ref_inh_B_alone = ref_inh_B
  )
}

if (length(ref_curve_compare) > 0) {
  curve_compare_df <- do.call(rbind, ref_curve_compare)
  write.csv(curve_compare_df, "D:/Tool_genomic/module_1/curve_fit_comparison.csv", row.names = FALSE)
  cat("Single-agent curve comparison (yours vs. synergyfinder\'s own fit) written to: curve_fit_comparison2.csv\n")
  
  ok_a <- is.finite(curve_compare_df$inh_A_alone) & is.finite(curve_compare_df$ref_inh_A_alone)
  ok_b <- is.finite(curve_compare_df$inh_B_alone) & is.finite(curve_compare_df$ref_inh_B_alone)
  if (sum(ok_a) >= 3) {
    cat(sprintf("\nDrug A curve agreement -- MAE: %.4f, Pearson r: %.3f (n=%d)\n",
                mean(abs(curve_compare_df$inh_A_alone[ok_a] - curve_compare_df$ref_inh_A_alone[ok_a])),
                cor(curve_compare_df$inh_A_alone[ok_a], curve_compare_df$ref_inh_A_alone[ok_a]),
                sum(ok_a)))
  }
  if (sum(ok_b) >= 3) {
    cat(sprintf("Drug B curve agreement -- MAE: %.4f, Pearson r: %.3f (n=%d)\n",
                mean(abs(curve_compare_df$inh_B_alone[ok_b] - curve_compare_df$ref_inh_B_alone[ok_b])),
                cor(curve_compare_df$inh_B_alone[ok_b], curve_compare_df$ref_inh_B_alone[ok_b]),
                sum(ok_b)))
  }
} else {
  message("Not enough single-agent points in any block to run the curve comparison.")
}

cat("Done.\n")