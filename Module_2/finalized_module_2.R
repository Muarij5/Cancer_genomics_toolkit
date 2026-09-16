# ============================================================
# CANCER COMBINATION ANALYSIS TOOL
# Module 2: Drug Response — Bliss & Loewe Synergy Calculator
# v3: added Viability / Inhibition input toggle with clear UX
# ============================================================

library(shiny)
library(ggplot2)
library(plotly)
library(DT)
library(drc)
#setwd("D:/Tool_genomic/finalized_glm_version.R")
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
  # NEW:
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
# ------------------------------------------------------------
# Core analysis engine
# ------------------------------------------------------------
analyze_single_pair <- function(df, unit, response_type = "viability") {
  
  drug_a <- unique(df$Drug_A)[1]
  drug_b <- unique(df$Drug_B)[1]
  cell_line <- unique(df$Cell_Line)[1]
  
  df$Dose_A_uM <- convert_to_uM(df$Dose_A, unit)
  df$Dose_B_uM <- convert_to_uM(df$Dose_B, unit)
  
  # Step 1: Find control value from ORIGINAL data (before averaging)
  # so we normalize to the real untreated control, not an assumed 100
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
    validate(need(
      is.finite(control_val) && control_val > 1,
      paste0("Control value (Dose_A=0, Dose_B=0) is ", round(control_val, 4),
             " \u2014 too small/invalid for a 0\u2013100% viability scale. ",
             "This dataset's Response column is not on the scale this tool expects. ",
             "Check the data source/units before trusting any result.")
    ))
  }
  # Step 2: Collapse replicate rows to mean response
  # average_replicates only keeps Dose_A, Dose_B, Response — that's fine,
  # we'll add inhibition right after this
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
  
  validate(need(nrow(rows_A_alone) >= 2, "Need at least 2 single-agent dose points for Drug A (Dose_B = 0)"))
  validate(need(nrow(rows_B_alone) >= 2, "Need at least 2 single-agent dose points for Drug B (Dose_A = 0)"))
  validate(need(nrow(rows_combo)   >= 1, "Need at least 1 combination row (both doses > 0)"))
  
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
  
  # NEW:
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
    "No exact Loewe value \u2014 neither drug alone can reach this effect level. Approximate value shown (SynergyFinder-style estimate).",
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
  validate(need(length(pair_keys) >= 1, "No valid drug-pair groups found in the uploaded file."))
  
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
  
  validate(need(length(results_by_pair) >= 1,
                "None of the drug-pair groups had enough data to analyze (need at least 2 single-agent points per drug and 1 combo point)."))
  
  results_by_pair
}

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
# UI
# ============================================================
ui <- fluidPage(
  tags$head(
    tags$title("CancerComb Analyzer"),
    tags$style(HTML("
      @import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@300;400;500;600;700&family=IBM+Plex+Mono:wght@400;500&display=swap');
      *,*::before,*::after{box-sizing:border-box;margin:0;padding:0;}
      body{font-family:'IBM Plex Sans',sans-serif;background:#f4f6f9;color:#1a1f2e;min-height:100vh;}

      .topbar{background:#1a1f2e;height:58px;display:flex;align-items:center;padding:0 32px;position:sticky;top:0;z-index:100;box-shadow:0 2px 8px rgba(0,0,0,.18);}
      .tool-logo{font-size:17px;font-weight:700;color:#fff;margin-right:40px;}
      .tool-logo span{color:#4f8ef7;}
      .nav-right{margin-left:auto;font-size:11px;color:#8892a4;}

      .page-content{max-width:1280px;margin:0 auto;padding:32px 28px;}
      .module-header{margin-bottom:28px;}
      .module-tag{display:inline-flex;align-items:center;gap:6px;background:#e8f0fe;color:#4f8ef7;font-size:11px;font-weight:600;padding:4px 12px;border-radius:20px;margin-bottom:10px;letter-spacing:.4px;text-transform:uppercase;}
      .module-title{font-size:26px;font-weight:700;letter-spacing:-.5px;margin-bottom:6px;}
      .module-desc{font-size:14px;color:#5a6378;}

      .two-col{display:grid;grid-template-columns:320px 1fr;gap:24px;align-items:start;}
      .input-panel{background:#fff;border:1px solid #e2e8f0;border-radius:12px;overflow:hidden;}
      .input-panel-header{background:#1a1f2e;padding:14px 20px;font-size:12px;font-weight:600;color:#8892a4;text-transform:uppercase;letter-spacing:1px;}
      .input-panel-body{padding:20px;display:flex;flex-direction:column;gap:20px;}
      .field-label{font-size:11px;font-weight:600;color:#8892a4;text-transform:uppercase;letter-spacing:.8px;margin-bottom:8px;}

      .upload-area{border:2px dashed #cbd5e1;border-radius:10px;padding:20px 14px;text-align:center;background:#f8fafc;}
      .upload-area:hover{border-color:#4f8ef7;background:#eff6ff;}
      .upload-icon-wrap{width:44px;height:44px;background:#e8f0fe;border-radius:10px;display:flex;align-items:center;justify-content:center;margin:0 auto 10px;font-size:20px;}
      .upload-hint{font-size:11px;color:#94a3b8;margin-top:6px;}
      .upload-area .form-group{margin:0;}
      .upload-area input[type=file]{font-size:12px;color:#5a6378;}
      .progress{display:none!important;}

      .form-control{border:1.5px solid #e2e8f0!important;border-radius:8px!important;font-size:13px!important;color:#1a1f2e!important;padding:9px 12px!important;background:#f8fafc!important;}
      .form-control:focus{border-color:#4f8ef7!important;box-shadow:0 0 0 3px #4f8ef722!important;}

      .btn-run{width:100%;background:#1a1f2e;border:none;color:#fff;font-size:13px;font-weight:600;padding:12px;border-radius:8px;cursor:pointer;display:flex;align-items:center;justify-content:center;gap:8px;}
      .btn-run:hover{background:#2d3550;}
      .btn-template{width:100%;background:transparent;border:1.5px solid #e2e8f0;color:#5a6378;font-size:12px;font-weight:500;padding:9px;border-radius:8px;cursor:pointer;}
      .btn-template:hover{border-color:#4f8ef7;color:#4f8ef7;background:#eff6ff;}

      .legend-box{background:#f8fafc;border:1px solid #e2e8f0;border-radius:8px;padding:14px;}
      .legend-row{display:flex;align-items:center;gap:8px;font-size:12px;color:#5a6378;margin-bottom:6px;line-height:1.4;}
      .legend-row:last-child{margin-bottom:0;}
      .legend-dot{width:8px;height:8px;border-radius:50%;flex-shrink:0;}

      /* Response type toggle - prominent styling */
      .response-type-panel{background:#f8fafc;border:1.5px solid #e2e8f0;border-radius:10px;padding:16px;}
      .response-type-panel .shiny-options-group{display:flex;flex-direction:column;gap:10px;}
      .response-type-panel .radio-inline{margin:0;padding:0;display:flex;align-items:flex-start;gap:10px;}
      .response-type-panel .radio-inline input[type=radio]{width:16px;height:16px;margin:2px 0 0 0;cursor:pointer;accent-color:#4f8ef7;flex-shrink:0;}
      .response-type-panel .radio-inline label{font-size:12px;font-weight:500;color:#1a1f2e;cursor:pointer;margin:0;line-height:1.4;}
      .response-type-panel .radio-inline label b{display:block;font-size:13px;font-weight:600;margin-bottom:2px;}
      .response-type-example{font-size:11px;color:#64748b;font-style:italic;margin-top:1px;}
      .response-type-note{font-size:11px;color:#475569;margin-top:12px;padding:10px;background:#eff6ff;border-radius:6px;line-height:1.5;border-left:3px solid #4f8ef7;}

      .results-panel{display:flex;flex-direction:column;gap:20px;min-height:400px;}
      .empty-state{background:#fff;border:1px solid #e2e8f0;border-radius:12px;padding:64px 32px;text-align:center;}
      .empty-icon{width:64px;height:64px;background:#f1f5f9;border-radius:16px;display:flex;align-items:center;justify-content:center;margin:0 auto 16px;font-size:28px;}
      .empty-title{font-size:16px;font-weight:600;margin-bottom:6px;}
      .empty-sub{font-size:13px;color:#94a3b8;max-width:360px;margin:0 auto;line-height:1.6;}

      .summary-card{border-radius:12px;padding:24px 28px;border:1px solid;position:relative;overflow:hidden;}
      .summary-card::before{content:'';position:absolute;top:0;left:0;right:0;height:4px;}
      .card-synergy{background:#f0fdf4;border-color:#bbf7d0;} .card-synergy::before{background:linear-gradient(90deg,#22c55e,#16a34a);}
      .card-antagonism{background:#fff1f2;border-color:#fecdd3;} .card-antagonism::before{background:linear-gradient(90deg,#ef4444,#dc2626);}
      .card-additive{background:#fffbeb;border-color:#fde68a;} .card-additive::before{background:linear-gradient(90deg,#f59e0b,#d97706);}
      .card-conflicting{background:#faf5ff;border-color:#e9d5ff;} .card-conflicting::before{background:linear-gradient(90deg,#a855f7,#9333ea);}

      .card-top-row{display:flex;align-items:flex-start;justify-content:space-between;margin-bottom:20px;}
      .combo-label{font-size:11px;font-weight:600;color:#94a3b8;text-transform:uppercase;letter-spacing:.8px;margin-bottom:4px;}
      .combo-name{font-size:22px;font-weight:700;font-family:'IBM Plex Mono',monospace;letter-spacing:-.5px;}
      .verdict-badge{font-size:13px;font-weight:700;padding:8px 20px;border-radius:8px;letter-spacing:.3px;text-transform:uppercase;white-space:nowrap;}
      .badge-synergy{background:#22c55e;color:#fff;} .badge-antagonism{background:#ef4444;color:#fff;}
      .badge-additive{background:#f59e0b;color:#fff;} .badge-conflicting{background:#a855f7;color:#fff;}

      .meta-row{font-size:12px;color:#94a3b8;margin-bottom:20px;display:flex;gap:16px;flex-wrap:wrap;}
      .meta-item{display:flex;align-items:center;gap:4px;}

      .scores-grid{display:grid;grid-template-columns:1fr 1fr;gap:14px;margin-bottom:20px;}
      .score-tile{background:#fff;border:1px solid #e2e8f0;border-radius:10px;padding:16px 18px;}
      .score-tile-label{font-size:11px;font-weight:600;color:#94a3b8;text-transform:uppercase;letter-spacing:.8px;margin-bottom:6px;}
      .score-tile-value{font-size:32px;font-weight:700;font-family:'IBM Plex Mono',monospace;letter-spacing:-1px;line-height:1;margin-bottom:4px;}
      .score-tile-class{font-size:12px;font-weight:600;}
      .col-synergy{color:#16a34a;} .col-antagonism{color:#dc2626;} .col-additive{color:#d97706;} .col-conflicting{color:#9333ea;}

      .interp-box{background:#fff;border:1px solid #e2e8f0;border-radius:8px;padding:14px 16px;font-size:13px;color:#5a6378;line-height:1.7;}
      .conflict-note{background:#faf5ff;border:1px solid #e9d5ff;border-radius:8px;padding:12px 16px;font-size:12px;color:#7e22ce;line-height:1.6;margin-top:12px;}

      .section-card{background:#fff;border:1px solid #e2e8f0;border-radius:12px;overflow:hidden;}
      .section-card-header{padding:14px 20px;border-bottom:1px solid #f1f5f9;display:flex;align-items:center;justify-content:space-between;}
      .section-card-title{font-size:13px;font-weight:600;}
      .section-card-sub{font-size:11px;color:#94a3b8;margin-top:1px;}
      .section-card-body{padding:20px;}

      .btn-download{background:#eff6ff;border:1.5px solid #bfdbfe;color:#2563eb;font-size:12px;font-weight:600;padding:7px 16px;border-radius:6px;cursor:pointer;}
      .btn-download:hover{background:#dbeafe;}

      .dataTables_wrapper{color:#1a1f2e!important;font-size:13px;}
      table.dataTable thead th{background:#f8fafc!important;color:#64748b!important;font-size:11px!important;font-weight:600!important;text-transform:uppercase!important;letter-spacing:.5px!important;border-bottom:1px solid #e2e8f0!important;padding:10px 14px!important;}
      table.dataTable tbody td{padding:10px 14px!important;border-bottom:1px solid #f1f5f9!important;color:#1a1f2e!important;}
      table.dataTable tbody tr:hover td{background:#f8fafc!important;}
      .dataTables_info,.dataTables_length,.dataTables_filter,.dataTables_paginate{font-size:12px!important;color:#94a3b8!important;}
      .dataTables_filter input{border:1px solid #e2e8f0!important;border-radius:6px!important;padding:5px 10px!important;font-size:12px!important;color:#1a1f2e!important;}
      .paginate_button{color:#64748b!important;border-radius:6px!important;}
      .paginate_button.current{background:#1a1f2e!important;color:#fff!important;border-color:#1a1f2e!important;}
    "))
  ),
  
  div(class = "topbar",
      div(class = "tool-logo", "Cancer", tags$span("Comb"), " Analyzer"),
      div(class = "nav-right", "Module 2 · Drug Combination Response")
  ),
  
  div(class = "page-content",
      div(class = "module-header",
          div(class = "module-tag", "Module 02 · Drug Response"),
          div(class = "module-title", "Drug Combination Response Analyzer"),
          div(class = "module-desc",
              "Enter two drugs, their doses, and the measured response — the tool calculates Bliss & Loewe synergy scores")
      ),
      div(class = "two-col",
          
          div(class = "input-panel",
              div(class = "input-panel-header", "Input Configuration"),
              div(class = "input-panel-body",
                  
                  div(
                    div(class = "field-label", "Step 1 — Get Template"),
                    downloadButton("downloadTemplate", "⬇  Download CSV Template", class = "btn-template")
                  ),
                  
                  div(
                    div(class = "field-label", "Step 2 — Upload Your Data"),
                    div(class = "upload-area",
                        div(class = "upload-icon-wrap", "📂"),
                        fileInput("fileInput", NULL, accept = ".csv",
                                  placeholder = "Choose CSV file", buttonLabel = "Browse"),
                        div(class = "upload-hint", "Drug A, Drug B, Cell line, Dose A, Dose B, Response")
                    )
                  ),
                  
                  # ============================================
                  # RESPONSE TYPE TOGGLE — CLEAR & PROMINENT
                  # ============================================
                  div(
                    div(class = "field-label", "Step 3 — What Does Your Response Column Contain?"),
                    div(class = "response-type-panel",
                        radioButtons("responseType", NULL,
                                     choices = c(
                                       "viability" = "viability",
                                       "inhibition" = "inhibition"
                                     ),
                                     selected = "viability"),
                        # Custom labels with examples — much clearer than default
                        tags$div(style = "margin-top:-4px;",
                                 tags$div(class = "response-type-example", 
                                          style = "padding-left:26px;color:#16a34a;",
                                          "● Untreated control (0,0) ≈ 100% — values DROP with dose"),
                                 tags$div(class = "response-type-example", 
                                          style = "padding-left:26px;color:#2563eb;margin-top:6px;",
                                          "● Untreated control (0,0) ≈ 0% — values RISE with dose")
                        ),
                        div(class = "response-type-note",
                            "ℹ️ <b>Both are converted to 0–1 inhibition internally</b> before Bliss & Loewe are calculated. If your file has a (0,0) control row, we normalize to that value; otherwise we assume 100%.")
                    )
                  ),
                  
                  div(
                    div(class = "field-label", "Step 4 — Dose Unit"),
                    selectInput("unit", NULL, choices = c("nM","uM","mM","mg/mL"), selected = "nM")
                  ),
                  
                  div(
                    div(class = "field-label", "Step 5 — Run"),
                    actionButton("analyze", "▶  Run Analysis", class = "btn-run")
                  ),
                  
                  div(
                    div(class = "field-label", "Verdict Logic"),
                    div(class = "legend-box",
                        div(class = "legend-row", div(class="legend-dot", style="background:#22c55e;"),
                            HTML(paste0("<b style='color:#16a34a'>Synergy</b> — Bliss AND Loewe both &gt; +", SYN_THRESHOLD))),
                        div(class = "legend-row", div(class="legend-dot", style="background:#f59e0b;"),
                            HTML(paste0("<b style='color:#d97706'>Additive</b> — Bliss AND Loewe both within ±", SYN_THRESHOLD))),
                        div(class = "legend-row", div(class="legend-dot", style="background:#ef4444;"),
                            HTML(paste0("<b style='color:#dc2626'>Antagonism</b> — Bliss AND Loewe both &lt; -", SYN_THRESHOLD))),
                        div(class = "legend-row", div(class="legend-dot", style="background:#a855f7;"),
                            HTML("<b style='color:#9333ea'>Conflicting</b> — Bliss and Loewe disagree")),
                        div(style="font-size:11px;color:#94a3b8;margin-top:8px;line-height:1.5;",
                            paste0("±", SYN_THRESHOLD, " threshold · Models classified independently · Never averaged"))
                    )
                  )
              )
          ),
          
          div(class = "results-panel",
              uiOutput("pairSelectorUI"),
              uiOutput("resultsUI"))
      )
  )
)

# ============================================================
# SERVER
# ============================================================
server <- function(input, output, session) {
  
  output$downloadTemplate <- downloadHandler(
    filename = "drug_response_template.csv",
    content  = function(file) writeLines(template_csv, file)
  )
  
  raw_data <- reactive({
    req(input$fileInput)
    tryCatch(read.csv(input$fileInput$datapath, stringsAsFactors = FALSE),
             error = function(e) NULL)
  })
  
  all_results <- eventReactive(input$analyze, {
    df <- raw_data()
    req(!is.null(df))
    
    required_cols <- c("Drug_A","Drug_B","Cell_Line","Dose_A","Dose_B","Response")
    missing <- setdiff(required_cols, colnames(df))
    validate(need(length(missing) == 0,
                  paste("Missing columns:", paste(missing, collapse = ", "))))
    
    # Pass the response_type toggle value through to the analysis
    analyze_combination(df, input$unit, input$responseType)
  })
  
  pair_summaries <- reactive({
    res_list <- all_results()
    lapply(names(res_list), function(pk) {
      r <- res_list[[pk]]
      detail <- r$detail
      avg_bliss <- round(mean(detail$Bliss_Score, na.rm = TRUE), 2)
      avg_loewe <- round(mean(detail$Loewe_Score, na.rm = TRUE), 2)
      bliss_class <- classify_score(avg_bliss)
      loewe_class <- classify_score(avg_loewe)
      n_conflicting_rows <- sum(detail$Consensus == "Conflicting", na.rm = TRUE)
      list(
        pair_key      = pk,
        detail        = detail,
        avg_bliss     = avg_bliss,
        avg_loewe     = avg_loewe,
        bliss_class   = bliss_class,
        loewe_class   = loewe_class,
        consensus     = consensus_verdict(bliss_class, loewe_class),
        n_conflicting_rows = n_conflicting_rows,
        drug_a        = r$drug_a,
        drug_b        = r$drug_b,
        cell_line     = r$cell_line,
        unit          = input$unit,
        response_type = input$responseType,
        control_used  = r$control_used,
        n             = r$n_combo,
        curve_A_fit   = r$curve_A_fit,
        curve_B_fit   = r$curve_B_fit
      )
    })
  })
  
  output$pairSelectorUI <- renderUI({
    summaries <- pair_summaries()
    if (length(summaries) <= 1) return(NULL)
    
    choices <- setNames(
      sapply(summaries, function(s) s$pair_key),
      sapply(summaries, function(s) paste0(s$drug_a, " + ", s$drug_b, " (", s$cell_line, ") — ", s$consensus))
    )
    div(
      div(class = "field-label", paste0("Drug Pairs Found (", length(summaries), ")")),
      selectInput("selectedPair", NULL, choices = choices)
    )
  })
  
  results <- reactive({
    summaries <- pair_summaries()
    if (length(summaries) == 1) return(summaries[[1]])
    sel <- input$selectedPair
    if (is.null(sel)) return(summaries[[1]])
    match_idx <- which(sapply(summaries, function(s) s$pair_key) == sel)
    if (length(match_idx) == 0) return(summaries[[1]])
    summaries[[match_idx[1]]]
  })
  
  output$resultsUI <- renderUI({
    if (is.null(input$fileInput) || input$analyze == 0) {
      return(div(class = "empty-state",
                 div(class = "empty-icon", "🔬"),
                 div(class = "empty-title", "Ready for Analysis"),
                 div(class = "empty-sub",
                     "Download the template, fill in Drug A, Drug B, Cell line, doses and response, upload, then click Run Analysis")
      ))
    }
    
    res <- tryCatch(results(), error = function(e) NULL)
    if (is.null(res)) {
      return(div(class = "empty-state",
                 div(class = "empty-icon", "⚠️"),
                 div(class = "empty-title", "Check your file format"),
                 div(class = "empty-sub", "Make sure all required columns are present and include single-agent rows (one dose = 0).")
      ))
    }
    
    card_suffix <- switch(res$consensus,
                          "Synergy"="synergy","Antagonism"="antagonism",
                          "Additive"="additive","conflicting")
    card_class  <- paste("summary-card", paste0("card-", card_suffix))
    badge_class <- paste("verdict-badge", paste0("badge-", card_suffix))
    
    col_class <- function(cls) switch(cls,"Synergy"="col-synergy","Antagonism"="col-antagonism",
                                      "Additive"="col-additive","col-conflicting")
    bliss_col <- col_class(res$bliss_class)
    loewe_col <- col_class(res$loewe_class)
    
    interp <- switch(res$consensus,
                     "Synergy" = paste0(res$drug_a," + ",res$drug_b," shows synergistic effect on ",res$cell_line,
                                        " across ",res$n," dose pairs. Both Bliss (",res$avg_bliss,
                                        ") and Loewe (",res$avg_loewe,") agree the combination kills more cells than expected."),
                     "Antagonism" = paste0(res$drug_a," + ",res$drug_b," shows antagonistic effect on ",res$cell_line,
                                           " across ",res$n," dose pairs. Both Bliss (",res$avg_bliss,
                                           ") and Loewe (",res$avg_loewe,") agree the drugs interfere with each other."),
                     "Additive" = paste0(res$drug_a," + ",res$drug_b," shows an additive effect on ",res$cell_line,
                                         ". Both Bliss (",res$avg_bliss,") and Loewe (",res$avg_loewe,
                                         ") fall within the no-interaction zone (±",SYN_THRESHOLD,")."),
                     paste0(res$drug_a," + ",res$drug_b," gives a CONFLICTING signal on ",res$cell_line,
                            ": Bliss (",res$avg_bliss," — ",res$bliss_class,") and Loewe (",res$avg_loewe," — ",res$loewe_class,
                            ") disagree. Review the heatmap and consider additional reference models before drawing conclusions.")
    )
    
    tagList(
      div(class = card_class,
          div(class = "card-top-row",
              div(div(class="combo-label","Drug Combination"),
                  div(class="combo-name", paste(res$drug_a,"+",res$drug_b))),
              div(class = badge_class, res$consensus)
          ),
          div(class = "meta-row",
              div(class="meta-item","🧫 Cell line:", tags$b(res$cell_line)),
              div(class="meta-item","🔢 Dose pairs:", tags$b(res$n)),
              div(class="meta-item","📐 Unit:", tags$b(res$unit)),
              div(class="meta-item","🧪 Input type:", tags$b(ifelse(res$response_type == "inhibition", "Inhibition", "Viability"))),
              div(class="meta-item","📊 Control:", tags$b(round(res$control_used, 1), "%")),
              div(class="meta-item","📈 Curve:", tags$b(paste0("A=",ifelse(res$curve_A_fit,"LL.4","interp"), " / B=",ifelse(res$curve_B_fit,"LL.4","interp"))))
          ),
          div(class = "scores-grid",
              div(class="score-tile",
                  div(class="score-tile-label","Bliss Score"),
                  div(class=paste("score-tile-value",bliss_col), res$avg_bliss),
                  div(class=paste("score-tile-class",bliss_col), res$bliss_class)),
              div(class="score-tile",
                  div(class="score-tile-label","Loewe Score"),
                  div(class=paste("score-tile-value",loewe_col), res$avg_loewe),
                  div(class=paste("score-tile-class",loewe_col), res$loewe_class))
          ),
          div(class = "interp-box", interp),
          if (res$consensus != "Conflicting" && res$n_conflicting_rows > 0)
            div(class = "conflict-note",
                paste0("⚠ ", res$n_conflicting_rows, " of ", res$n,
                       " individual dose pairs showed conflicting signals even though the overall average agrees."))
      ),
      
      div(class = "section-card",
          div(class = "section-card-header",
              div(div(class="section-card-title","Bliss Score Heatmap"),
                  div(class="section-card-sub","Synergy landscape across dose pairs"))),
          div(class = "section-card-body",
              plotOutput("heatmap", height = "300px"))
      ),
      
      div(class = "section-card",
          div(class = "section-card-header",
              div(div(class="section-card-title","Detailed Results"),
                  div(class="section-card-sub","Per dose pair breakdown")),
              downloadButton("downloadResults", "⬇ Export CSV", class = "btn-download")),
          div(class = "section-card-body", DTOutput("resultsTable"))
      )
    )
  })
  
  output$heatmap <- renderPlot({
    res <- results()
    df <- res$detail
    
    ggplot(df, aes(x = factor(Dose_A), y = factor(Dose_B), fill = Bliss_Score)) +
      geom_tile(color = "#ffffff", linewidth = 1) +
      # Dynamic text color: dark on light backgrounds, white on dark
      geom_text(aes(label = round(Bliss_Score, 1), 
                    color = abs(Bliss_Score) > 15), 
                size = 4.5, fontface = "bold", show.legend = FALSE) +
      scale_fill_gradient2(low = "#dc2626", mid = "#f1f5f9", high = "#16a34a", 
                           midpoint = 0, name = "Bliss") +
      scale_color_manual(values = c("TRUE" = "white", "FALSE" = "#1a1f2e")) +
      labs(x = paste0(res$drug_a, " dose (", res$unit, ")"),
           y = paste0(res$drug_b, " dose (", res$unit, ")")) +
      theme_minimal(base_family = "IBM Plex Sans") +
      theme(plot.background = element_rect(fill = "#fff", color = NA),
            panel.background = element_rect(fill = "#fff", color = NA),
            panel.grid = element_blank(),
            axis.text = element_text(color = "#64748b", size = 11),
            axis.title = element_text(color = "#64748b", size = 11),
            legend.text = element_text(color = "#64748b"),
            legend.title = element_text(color = "#64748b"))
  }, height = 300, width = 600)
  
  output$resultsTable <- renderDT({
    res <- results()
    df  <- res$detail
    
    # Add a note column explaining WHY Loewe is blank when it's NA
    display <- df[, c("Drug_A","Drug_B","Cell_Line","Dose_A","Dose_B","Response",
                      "Bliss_Score","Bliss_Class",
                      "Loewe_Score","Loewe_Class","Loewe_Approx_Score","Loewe_Approx_Class",
                      "Loewe_Note","Consensus")]
    colnames(display) <- c("Drug A","Drug B","Cell Line","Dose A","Dose B","Response (%)",
                           "Bliss Score","Bliss Class",
                           "Loewe Score (exact)","Loewe Class (exact)",
                           "Loewe Score (approx.)","Loewe Class (approx.)",
                           "Note","Verdict")
    
    datatable(display,
              options = list(pageLength = 8, scrollX = TRUE, dom = "ftp", autoWidth = FALSE),
              rownames = FALSE) %>%
      formatStyle("Verdict",
                  backgroundColor = styleEqual(c("Synergy","Antagonism","Additive","Conflicting"),
                                               c("#f0fdf4","#fff1f2","#fffbeb","#faf5ff")),
                  color = styleEqual(c("Synergy","Antagonism","Additive","Conflicting"),
                                     c("#16a34a","#dc2626","#d97706","#9333ea")),
                  fontWeight = "bold") %>%
      formatStyle("Note",
                  color = "#9333ea",
                  fontStyle = "italic") %>%
      formatStyle("Bliss Class",
                  color = styleEqual(c("Synergy","Antagonism","Additive"), c("#16a34a","#dc2626","#d97706")),
                  fontWeight = "bold") %>%
      formatStyle("Loewe Class (exact)",
                  color = styleEqual(c("Synergy","Antagonism","Additive"), c("#16a34a","#dc2626","#d97706")),
                  fontWeight = "bold")
  })
  
  output$downloadResults <- downloadHandler(
    filename = function() paste0(results()$drug_a,"_",results()$drug_b,"_results.csv"),
    content  = function(file) write.csv(results()$detail, file, row.names = FALSE)
  )
}

shinyApp(ui = ui, server = server)