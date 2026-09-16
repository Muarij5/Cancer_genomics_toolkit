# ============================================================
# CANCER COMBINATION ANALYSIS TOOL
# Module 2: Drug Response — Bliss & Loewe Synergy Calculator
# PIPELINE VERSION — Shiny app structure removed.
# Same calculation logic as finalized_module_2.R, callable as a
# plain R script/function instead of a web app.
# ============================================================

library(ggplot2)
library(drc)

# ============================================================
# CONSTANTS
# ============================================================

SYN_THRESHOLD <- 10
MIN_POINTS_FOR_CURVE_FIT <- 3

# ============================================================
# CALCULATION FUNCTIONS
# ============================================================

convert_to_uM <- function(value, unit) {
  switch(unit,
         "nM"    = value / 1000,
         "uM"    = value,
         "mM"    = value * 1000,
         "mg/mL" = value * 1000,
         value
  )
}

# ------------------------------------------------------------
# Response type conversion
# ------------------------------------------------------------
# Converts raw Response values into 0–1 inhibition space.
#   "viability"  : control ≈ 100, drops with dose → inhibition = 1 - (resp/100)
#   "inhibition" : control ≈ 0, rises with dose   → inhibition = resp/100
#
# If a true control row (Dose_A=0, Dose_B=0) exists, we normalize to
# that actual value rather than assuming exactly 100, which is more
# robust to solvent effects and assay drift.

# ------------------------------------------------------------
# Replicate handling
# ------------------------------------------------------------
average_replicates <- function(df) {
  agg <- aggregate(Response ~ Dose_A + Dose_B, data = df, FUN = mean)
  agg
}

# ------------------------------------------------------------
# Single-agent curve fitting
# ------------------------------------------------------------
fit_single_agent_curve <- function(conc_vals, inhibition_vals) {
  
  raw_predict <- function(dose) {
    sapply(dose, function(d) {
      if (any(abs(d - conc_vals) < 1e-9)) return(inhibition_vals[which.min(abs(d - conc_vals))])
      approx(conc_vals, inhibition_vals, xout = d, rule = 2)$y
    })
  }
  
  if (length(unique(conc_vals)) < MIN_POINTS_FOR_CURVE_FIT) {
    idx <- which.min(abs(inhibition_vals - 0.5))
    fallback_ec50 <- if (length(idx) > 0) conc_vals[idx] else mean(conc_vals, na.rm = TRUE)
    return(list(ec50 = fallback_ec50, predict_fn = raw_predict, used_fit = FALSE, inverse_fn = NULL))  }
  
  # Fit LL.4 on VIABILITY (1 - inhibition), which decreases with dose --
  # this matches the direction LL.4 is built for. Fitting inhibition directly
  # (which increases with dose) fights the model's assumed shape and causes
  # bad/flat/step-like curves, especially with sparse single-agent data.
  viability_vals <- 1 - inhibition_vals
  
  fit <- tryCatch(
    drc::drm(viability_vals ~ conc_vals,
             fct = drc::LL.4(names = c("Slope", "Lower", "Upper", "EC50"))),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    idx <- which.min(abs(inhibition_vals - 0.5))
    fallback_ec50 <- if (length(idx) > 0) conc_vals[idx] else mean(conc_vals, na.rm = TRUE)
    return(list(ec50 = fallback_ec50, predict_fn = raw_predict, used_fit = FALSE))
  }
  
  ec50 <- abs(coef(fit)[4])
  if (is.na(ec50) || ec50 <= 0) {
    idx <- which.min(abs(inhibition_vals - 0.5))
    fallback_ec50 <- if (length(idx) > 0) conc_vals[idx] else mean(conc_vals, na.rm = TRUE)
    return(list(ec50 = fallback_ec50, predict_fn = raw_predict, used_fit = FALSE))
  }
  
  fit_predict <- function(dose) {
    pred_viability <- tryCatch(
      as.numeric(predict(fit, newdata = data.frame(conc_vals = dose))),
      error = function(e) NA
    )
    if (anyNA(pred_viability)) return(raw_predict(dose))
    pred_inhibition <- 1 - pred_viability
    pmin(pred_inhibition, 1)
  }
  # Analytic inverse of the fitted LL.4 curve: given target inhibition,
  # solve for dose directly from the model equation instead of interpolating
  # over a numeric grid. This is what makes Loewe's dose-inversion stable
  # even when the curve is flat/saturated across the tested dose range.
  b_par <- coef(fit)[1]; c_par <- coef(fit)[2]; d_par <- coef(fit)[3]
  fit_inverse <- function(target_inhibition) {
    y <- 1 - target_inhibition  # convert target inhibition back to viability
    sapply(y, function(yy) {
      if (is.na(yy)) return(NA_real_)
      base <- (d_par - yy) / (yy - c_par)
      if (is.na(base) || base <= 0 || !is.finite(base)) return(NA_real_)
      dose <- ec50 * base^(1 / b_par)
      if (!is.finite(dose)) return(NA_real_)
      dose
    })
  }
  
  list(ec50 = ec50, predict_fn = fit_predict, used_fit = TRUE, inverse_fn = fit_inverse,
       b_par = b_par, c_par = c_par, d_par = d_par)
}

# ------------------------------------------------------------
# Approximate Loewe — SynergyFinder-style fallback
# ------------------------------------------------------------
# This does NOT solve the Loewe equation exactly. It scans candidate effect
# levels between what the two single agents alone can reach, and picks
# whichever candidate the actual combo dose sits closest to. This mirrors
# what the `synergyfinder` R package does internally (.SolveLoewe()) — it
# always returns a number, even when no exact answer exists.
# Used ONLY as a labeled fallback when the exact solve below returns NA —
# never presented as the real answer.
loewe_expected_inhibition_approx <- function(dose_A, dose_B, curve_A, curve_B,
                                             inv_dose_A_fn, inv_dose_B_fn,
                                             nsteps = 100) {
  if (!isTRUE(curve_A$used_fit) || !isTRUE(curve_B$used_fit)) return(NA_real_)
  
  bounds_A <- range(1 - curve_A$c_par, 1 - curve_A$d_par)
  bounds_B <- range(1 - curve_B$c_par, 1 - curve_B$d_par)
  min_y <- min(bounds_A[1], bounds_B[1])
  max_y <- max(bounds_A[2], bounds_B[2])
  if (!is.finite(min_y) || !is.finite(max_y) || min_y >= max_y) return(NA_real_)
  
  y_test <- seq(min_y, max_y, length.out = nsteps)
  best_dist <- Inf
  best_y <- NA_real_
  
  for (y in y_test) {
    xa <- inv_dose_A_fn(y)
    xb <- inv_dose_B_fn(y)
    if (is.na(xa) || is.na(xb) || xa <= 0 || xb <= 0) next
    w1 <- 1 / xa; w2 <- 1 / xb
    dist <- abs(dose_A * w1 + dose_B * w2 - 1) / sqrt(w1^2 + w2^2)
    if (is.finite(dist) && dist < best_dist) {
      best_dist <- dist
      best_y <- y
    }
  }
  best_y
}

# ------------------------------------------------------------
# Classification
# ------------------------------------------------------------
classify_score <- function(s, threshold = SYN_THRESHOLD) {
  ifelse(s > threshold, "Synergy", ifelse(s < -threshold, "Antagonism", "Additive"))
}

consensus_verdict <- function(bliss_class, loewe_class) {
  ifelse(bliss_class == loewe_class, bliss_class, "Conflicting")
}

detect_pairs <- function(df) {
  if ("Pair_ID" %in% colnames(df)) {
    df$Pair_Key <- as.character(df$Pair_ID)
  } else {
    df$Pair_Key <- paste(df$Drug_A, df$Drug_B, df$Cell_Line, sep = " || ")
  }
  df
}

# ------------------------------------------------------------
# Core analysis engine
# ------------------------------------------------------------
# NOTE: the original app version used Shiny's validate(need(...)) here so
# a bad-file message would show up nicely in the UI. That only works inside
# a running Shiny session, so in this pipeline version those checks are
# plain stop() calls instead — same guardrails, just a normal R error.
analyze_single_pair <- function(df, unit, response_type = "viability") {
  
  drug_a <- unique(df$Drug_A)[1]
  drug_b <- unique(df$Drug_B)[1]
  cell_line <- unique(df$Cell_Line)[1]
  
  df$Dose_A_uM <- convert_to_uM(df$Dose_A, unit)
  df$Dose_B_uM <- convert_to_uM(df$Dose_B, unit)
  
  # Step 1: Find control value from ORIGINAL data (before averaging)
  # so we normalize to the real untreated control, not an assumed 100
  control_rows <- df[df$Dose_A == 0 & df$Dose_B == 0, ]
  control_val <- if (nrow(control_rows) > 0) {
    mean(control_rows$Response, na.rm = TRUE)
  } else {
    100
  }
  
  # Guard: for "viability" input, the untreated control should sit roughly on
  # a 0-100% scale. If it's near zero, negative, or otherwise implausible,
  # the Response column is NOT plain 0-100% viability for this dataset, and
  # dividing by it silently produces nonsense (thousands/billions) instead
  # of a real inhibition value. Stop with a clear error instead of that.
  if (identical(response_type, "viability")) {
    if (!(is.finite(control_val) && control_val > 1)) {
      stop(paste0("Control value (Dose_A=0, Dose_B=0) is ", round(control_val, 4),
                  " -- too small/invalid for a 0-100% viability scale. ",
                  "This dataset's Response column is not on the scale this tool expects. ",
                  "Check the data source/units before trusting any result."))
    }
  }
  
  # Step 2: Collapse replicate rows to mean response
  agg <- average_replicates(df[, c("Dose_A_uM", "Dose_B_uM", "Response")] |>
                              setNames(c("Dose_A", "Dose_B", "Response")))
  
  # Step 3: NOW convert the averaged Response → 0–1 inhibition
  # (doing it after averaging is correct — we want inhibition of the mean,
  # not the mean of individual inhibitions)
  if (identical(response_type, "inhibition")) {
    agg$inhibition <- agg$Response / 100
  } else {
    agg$inhibition <- 1 - (agg$Response / 100)
  }
  
  rows_A_alone <- agg[agg$Dose_B == 0 & agg$Dose_A > 0, ]
  rows_B_alone <- agg[agg$Dose_A == 0 & agg$Dose_B > 0, ]
  rows_combo   <- agg[agg$Dose_A > 0  & agg$Dose_B > 0, ]
  
  if (nrow(rows_A_alone) < 2) stop("Need at least 2 single-agent dose points for Drug A (Dose_B = 0)")
  if (nrow(rows_B_alone) < 2) stop("Need at least 2 single-agent dose points for Drug B (Dose_A = 0)")
  if (nrow(rows_combo)   < 1) stop("Need at least 1 combination row (both doses > 0)")
  
  curve_A <- fit_single_agent_curve(rows_A_alone$Dose_A, rows_A_alone$inhibition)
  curve_B <- fit_single_agent_curve(rows_B_alone$Dose_B, rows_B_alone$inhibition)
  
  # Bliss uses the RAW observed single-agent value at the matching dose, not a
  # fitted curve. Confirmed from synergyfinder's own source (Bliss() in
  # calculate_synergy_score.R): it pulls straight from ExtractSingleDrug(),
  # never fits a model. Only Loewe (below) needs a fitted curve.
  raw_lookup_fn <- function(rows_alone, dose_col) {
    function(doses) {
      sapply(doses, function(d) {
        idx <- which(abs(rows_alone[[dose_col]] - d) < 1e-9)
        if (length(idx) > 0) return(rows_alone$inhibition[idx[1]])
        nearest <- which.min(abs(rows_alone[[dose_col]] - d))
        rows_alone$inhibition[nearest]
      })
    }
  }
  inh_A_raw_fn <- raw_lookup_fn(rows_A_alone, "Dose_A")
  inh_B_raw_fn <- raw_lookup_fn(rows_B_alone, "Dose_B")
  
  rows_combo$inh_A_alone <- inh_A_raw_fn(rows_combo$Dose_A)
  rows_combo$inh_B_alone <- inh_B_raw_fn(rows_combo$Dose_B)
  
  bliss_expected <- rows_combo$inh_A_alone + rows_combo$inh_B_alone -
    (rows_combo$inh_A_alone * rows_combo$inh_B_alone)
  rows_combo$Bliss_Score <- round((rows_combo$inhibition - bliss_expected) * 100, 3)
  
  build_inverse_fn <- function(curve, conc_vals) {
    if (isTRUE(curve$used_fit) && !is.null(curve$inverse_fn)) {
      # Analytic LL.4 inverse — stable, extrapolates beyond the tested doses
      curve$inverse_fn
    } else {
      # Fallback for sparse data (no LL.4 fit): numeric lookup over the
      # actual observed points only.
      inh_vals <- curve$predict_fn(conc_vals)
      function(x) {
        sapply(x, function(xx) {
          if (is.na(xx) || xx <= 0) return(Inf)
          if (xx >= max(inh_vals)) return(max(conc_vals))
          approx(inh_vals, conc_vals, xout = xx, rule = 2)$y
        })
      }
    }
  }
  
  inv_dose_A_fn <- build_inverse_fn(curve_A, rows_A_alone$Dose_A)
  inv_dose_B_fn <- build_inverse_fn(curve_B, rows_B_alone$Dose_B)
  
  loewe_expected_inhibition <- function(dose_A, dose_B, inv_dose_A_fn, inv_dose_B_fn) {
    inv_dose_A <- function(x) { if (x <= 0) return(Inf); inv_dose_A_fn(x) }
    inv_dose_B <- function(x) { if (x <= 0) return(Inf); inv_dose_B_fn(x) }
    f <- function(x) {
      da <- inv_dose_A(x); db <- inv_dose_B(x)
      if (is.na(da) || is.na(db)) return(NA_real_)
      (dose_A / da) + (dose_B / db) - 1
    }
    
    tryCatch({
      probe_x <- c(1e-4, 0.01, 0.05, seq(0.1, 0.95, by = 0.05), 0.999)
      f_vals  <- sapply(probe_x, f)
      
      root <- NA_real_
      for (i in seq_len(length(probe_x) - 1)) {
        v1 <- f_vals[i]; v2 <- f_vals[i + 1]
        if (is.na(v1) || is.na(v2)) next
        if (sign(v1) != sign(v2)) {
          root <- uniroot(f, interval = c(probe_x[i], probe_x[i + 1]), tol = 1e-6)$root
          break
        }
      }
      if (is.na(root)) stop("no valid root in range")
      max(0, min(1, root))
    }, error = function(e) {
      NA
    })
  }
  
  rows_combo$loewe_expected <- mapply(loewe_expected_inhibition,
                                      rows_combo$Dose_A, rows_combo$Dose_B,
                                      MoreArgs = list(inv_dose_A_fn = inv_dose_A_fn, inv_dose_B_fn = inv_dose_B_fn))
  rows_combo$Loewe_Score <- round((rows_combo$inhibition - rows_combo$loewe_expected) * 100, 3)
  
  # Approximate fallback — ONLY computed/shown where the exact solve above
  # returned NA (i.e. no dose exists for either single agent to reach the
  # required effect). Exact result is always preferred when available.
  rows_combo$loewe_expected_approx <- mapply(
    loewe_expected_inhibition_approx,
    rows_combo$Dose_A, rows_combo$Dose_B,
    MoreArgs = list(curve_A = curve_A, curve_B = curve_B,
                    inv_dose_A_fn = inv_dose_A_fn, inv_dose_B_fn = inv_dose_B_fn)
  )
  rows_combo$Loewe_Approx_Score <- round((rows_combo$inhibition - rows_combo$loewe_expected_approx) * 100, 3)
  rows_combo$Loewe_Approx_Score[!is.na(rows_combo$Loewe_Score)] <- NA
  rows_combo$Loewe_Approx_Class <- ifelse(is.na(rows_combo$Loewe_Approx_Score), NA,
                                          classify_score(rows_combo$Loewe_Approx_Score))
  
  rows_combo$Loewe_Note <- ifelse(
    is.na(rows_combo$Loewe_Score) & !is.na(rows_combo$Loewe_Approx_Score),
    "No exact Loewe value -- neither drug alone can reach this effect level. Approximate value shown (SynergyFinder-style estimate).",
    ifelse(is.na(rows_combo$Loewe_Score), "No exact or approximate Loewe value available.", "")
  )
  
  rows_combo$Bliss_Class <- classify_score(rows_combo$Bliss_Score)
  rows_combo$Loewe_Class <- classify_score(rows_combo$Loewe_Score)
  rows_combo$Consensus   <- mapply(consensus_verdict, rows_combo$Bliss_Class, rows_combo$Loewe_Class)
  
  rows_combo$Drug_A <- drug_a
  rows_combo$Drug_B <- drug_b
  rows_combo$Cell_Line <- cell_line
  
  list(
    detail       = rows_combo,
    ec50_A       = curve_A$ec50,
    ec50_B       = curve_B$ec50,
    curve_A_fit  = curve_A$used_fit,
    curve_B_fit  = curve_B$used_fit,
    n_combo      = nrow(rows_combo),
    control_used = control_val
  )
}

# Multi-pair dispatcher
analyze_combination <- function(df, unit, response_type = "viability") {
  df <- detect_pairs(df)
  pair_keys <- unique(df$Pair_Key)
  if (length(pair_keys) < 1) stop("No valid drug-pair groups found in the uploaded file.")
  
  results_by_pair <- list()
  for (pk in pair_keys) {
    sub_df <- df[df$Pair_Key == pk, ]
    res <- tryCatch(analyze_single_pair(sub_df, unit, response_type), error = function(e) NULL)
    if (!is.null(res)) {
      res$drug_a    <- unique(sub_df$Drug_A)[1]
      res$drug_b    <- unique(sub_df$Drug_B)[1]
      res$cell_line <- unique(sub_df$Cell_Line)[1]
      res$pair_key  <- pk
      results_by_pair[[pk]] <- res
    }
  }
  
  if (length(results_by_pair) < 1) {
    stop("None of the drug-pair groups had enough data to analyze (need at least 2 single-agent points per drug and 1 combo point).")
  }
  
  results_by_pair
}

# ------------------------------------------------------------
# Per-pair summary
# ------------------------------------------------------------
# Extracted from the app's server-side `pair_summaries` reactive. Takes one
# element of the analyze_combination() output list and collapses it to a
# single verdict for that drug pair (average Bliss/Loewe, consensus, etc).
summarize_pair <- function(res, unit, response_type) {
  detail <- res$detail
  avg_bliss <- round(mean(detail$Bliss_Score, na.rm = TRUE), 2)
  avg_loewe <- round(mean(detail$Loewe_Score, na.rm = TRUE), 2)
  bliss_class <- classify_score(avg_bliss)
  loewe_class <- classify_score(avg_loewe)
  n_conflicting_rows <- sum(detail$Consensus == "Conflicting", na.rm = TRUE)
  list(
    pair_key      = res$pair_key,
    detail        = detail,
    avg_bliss     = avg_bliss,
    avg_loewe     = avg_loewe,
    bliss_class   = bliss_class,
    loewe_class   = loewe_class,
    consensus     = consensus_verdict(bliss_class, loewe_class),
    n_conflicting_rows = n_conflicting_rows,
    drug_a        = res$drug_a,
    drug_b        = res$drug_b,
    cell_line     = res$cell_line,
    unit          = unit,
    response_type = response_type,
    control_used  = res$control_used,
    n             = res$n_combo,
    curve_A_fit   = res$curve_A_fit,
    curve_B_fit   = res$curve_B_fit
  )
}

# ------------------------------------------------------------
# Bliss heatmap plot
# ------------------------------------------------------------
# Extracted from the app's renderPlot("heatmap"). Takes one summarize_pair()
# result and returns a ggplot object (call print() or save it, e.g. with
# ggsave()) instead of rendering to a Shiny plotOutput.
plot_bliss_heatmap <- function(pair_summary) {
  df <- pair_summary$detail
  
  ggplot(df, aes(x = factor(Dose_A), y = factor(Dose_B), fill = Bliss_Score)) +
    geom_tile(color = "#ffffff", linewidth = 1) +
    geom_text(aes(label = round(Bliss_Score, 1),
                  color = abs(Bliss_Score) > 15),
              size = 4.5, fontface = "bold", show.legend = FALSE) +
    scale_fill_gradient2(low = "#dc2626", mid = "#f1f5f9", high = "#16a34a",
                         midpoint = 0, name = "Bliss") +
    scale_color_manual(values = c("TRUE" = "white", "FALSE" = "#1a1f2e")) +
    labs(x = paste0(pair_summary$drug_a, " dose (", pair_summary$unit, ")"),
         y = paste0(pair_summary$drug_b, " dose (", pair_summary$unit, ")")) +
    theme_minimal() +
    theme(panel.grid = element_blank())
}

# ------------------------------------------------------------
# Template CSV — useful as a schema reference / for making test data
# ------------------------------------------------------------
template_csv <- "Drug_A,Drug_B,Cell_Line,Dose_A,Dose_B,Response
Everolimus,Dactolisib,BT-549,0,0,90
Everolimus,Dactolisib,BT-549,0,1,95
Everolimus,Dactolisib,BT-549,0,10,81
Everolimus,Dactolisib,BT-549,0,30,59
Everolimus,Dactolisib,BT-549,0,100,39
Everolimus,Dactolisib,BT-549,1,0,96
Everolimus,Dactolisib,BT-549,10,0,80
Everolimus,Dactolisib,BT-549,30,0,65
Everolimus,Dactolisib,BT-549,100,0,45
Everolimus,Dactolisib,BT-549,1,1,88
Everolimus,Dactolisib,BT-549,10,10,50
Everolimus,Dactolisib,BT-549,30,30,21
Everolimus,Dactolisib,BT-549,100,100,8"

# ============================================================
# PIPELINE ENTRY POINT
# ============================================================
# run_pipeline(): the one function you call from the outside.
#
#   csv_path      — path to a CSV with columns:
#                   Drug_A, Drug_B, Cell_Line, Dose_A, Dose_B, Response
#                   (optionally Pair_ID to group multiple drug pairs)
#   unit          — dose unit in the file: "nM", "uM", "mM", or "mg/mL"
#   response_type — "viability" (control ~100%, drops with dose) or
#                    "inhibition" (control ~0%, rises with dose)
#
# Returns a list of per-pair summaries (see summarize_pair()). Each has
# $avg_bliss, $avg_loewe, $consensus ("Synergy"/"Antagonism"/"Additive"/
# "Conflicting"), plus $detail (the full per-dose-pair table).
run_pipeline <- function(csv_path, unit = "uM", response_type = "viability") {
  df <- read.csv(csv_path, stringsAsFactors = FALSE)
  
  required_cols <- c("Drug_A", "Drug_B", "Cell_Line", "Dose_A", "Dose_B", "Response")
  missing <- setdiff(required_cols, colnames(df))
  if (length(missing) > 0) {
    stop(paste("Missing columns:", paste(missing, collapse = ", ")))
  }
  
  results_by_pair <- analyze_combination(df, unit, response_type)
  
  lapply(results_by_pair, summarize_pair, unit = unit, response_type = response_type)
}

# ============================================================
# EXAMPLE USAGE (uncomment to run)
# ============================================================
# writeLines(template_csv, "template.csv")
# summaries <- run_pipeline("template.csv", unit = "nM", response_type = "viability")
# for (s in summaries) {
#   cat(sprintf("%s + %s (%s): Bliss=%.2f Loewe=%.2f -> %s\n",
#               s$drug_a, s$drug_b, s$cell_line, s$avg_bliss, s$avg_loewe, s$consensus))
# }
# print(plot_bliss_heatmap(summaries[[1]]))
# write.csv(summaries[[1]]$detail, "results.csv", row.names = FALSE)
