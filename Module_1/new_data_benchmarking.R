# ============================================================
# run_external_benchmark.R  (v2 — matched to Parrish et al. 2021,
# Cell Rep. 36(9):109597, Table S4 — pgPEN dual-KO screen,
# PC9 + HeLa cells. Independent of DepMap.)
#
# Run AFTER Module2_finalized.R is sourced (APP_DATA loaded).
# Needs external_pairs.csv (built from Table S4) in DATA_DIR.
#
# Ground truth used here is the paper's OWN consolidated call
# (GI_flag column) -- no cutoffs were guessed by me:
#   buffering        -> fr_label = FR_positive
#   synthetic_lethal  -> sl_label = True
#   neither           -> FR_negative / False
# ============================================================

DATA_DIR <- "D:/Tool_genomic/depmap_data"
EXTERNAL_PAIRS_PATH <- file.path(DATA_DIR, "external_pairs.csv")

ext <- read.csv(EXTERNAL_PAIRS_PATH, stringsAsFactors = FALSE)
stopifnot(all(c("gene_a", "gene_b", "fr_label", "sl_label") %in% colnames(ext)))

n_pairs <- nrow(ext)
cat(sprintf("\n=== Starting external benchmark: %d pairs to process ===\n", n_pairs))
cat(sprintf("Estimated time at ~13 sec/pair: ~%.1f minutes\n\n", n_pairs * 13 / 60))

results <- vector("list", n_pairs)
t_start <- Sys.time()
for (i in seq_len(n_pairs)) {
  ga <- ext$gene_a[i]; gb <- ext$gene_b[i]
  
  elapsed <- as.numeric(Sys.time() - t_start, units = "mins")
  cat(sprintf("[%d/%d] %s + %s  (elapsed: %.1f min)\n", i, n_pairs, ga, gb, elapsed))
  
  res <- tryCatch(
    analyze_gene_pair_v2(ga, gb, APP_DATA$dependency_df, APP_DATA$mutation_df, APP_DATA$expression_df,
                         APP_DATA$paralog_df, compounds_df = APP_DATA$compounds_df,
                         model_df = APP_DATA$model_df, cn_df = APP_DATA$cn_df, cgc_df = APP_DATA$cgc_df),
    error = function(e) list(error = conditionMessage(e))
  )
  if (!is.null(res$error)) {
    results[[i]] <- data.frame(gene_a = ga, gene_b = gb, fr_call = NA, sl_call = NA, error = res$error)
  } else {
    results[[i]] <- data.frame(gene_a = ga, gene_b = gb, fr_call = res$FR$call, sl_call = res$SL$call, error = NA)
  }
}
cat(sprintf("\n=== Done: %d/%d pairs processed in %.1f minutes ===\n\n",
            n_pairs, n_pairs, as.numeric(Sys.time() - t_start, units = "mins")))

out_ext <- do.call(rbind, results)
out_ext <- merge(ext, out_ext, by = c("gene_a", "gene_b"))

n_errored <- sum(!is.na(out_ext$error))
cat(sprintf("\n%d / %d external pairs failed to run (missing genes / coverage gaps) -- excluded from scoring.\n",
            n_errored, nrow(out_ext)))
if (n_errored > 0) {
  cat("Top reasons (first 10 unique error messages):\n")
  print(head(unique(out_ext$error[!is.na(out_ext$error)]), 10))
}
scored_ext <- out_ext[is.na(out_ext$error), ]

score_ext <- function(pred_pos, truth_pos, label) {
  tp <- sum(pred_pos & truth_pos); fp <- sum(pred_pos & !truth_pos)
  fn <- sum(!pred_pos & truth_pos); tn <- sum(!pred_pos & !truth_pos)
  cat(sprintf("\n-- %s --\n", label))
  cat(sprintf("TP=%d FP=%d FN=%d TN=%d\n", tp, fp, fn, tn))
  cat(sprintf("Precision=%.3f  Recall=%.3f  Specificity=%.3f\n",
              tp / (tp + fp), tp / (tp + fn), tn / (tn + fp)))
}

# ---- FR external scoring ----
fr_ext <- scored_ext
fr_ext$fr_true_positive <- fr_ext$fr_label == "FR_positive"

cat("\n=== EXTERNAL VALIDATION (Parrish 2021) -- FR contingency table ===\n")
print(table(Predicted = fr_ext$fr_call, Truth = fr_ext$fr_label))

score_ext(fr_ext$fr_call %in% c("Functional Redundancy", "Functional Redundancy (moderate)"),
          fr_ext$fr_true_positive, "EXTERNAL FR STRICT")
score_ext(fr_ext$fr_call %in% c("Functional Redundancy", "Functional Redundancy (moderate)", "Borderline FR"),
          fr_ext$fr_true_positive, "EXTERNAL FR LENIENT")

# ---- SL external scoring ----
sl_ext <- scored_ext
sl_ext$sl_true_positive <- toupper(trimws(sl_ext$sl_label)) == "TRUE"

cat("\n=== EXTERNAL VALIDATION (Parrish 2021) -- SL contingency table ===\n")
print(table(Predicted = sl_ext$sl_call, Truth = sl_ext$sl_label))

score_ext(sl_ext$sl_call %in% c("Synthetic Lethal", "Synthetic Lethal (moderate)"),
          sl_ext$sl_true_positive, "EXTERNAL SL STRICT")

write.csv(out_ext, "external_benchmark_results_0.5.csv", row.names = FALSE)
cat("\nSaved: external_benchmark_results.csv\n")
cat("This is TRUE independent validation -- Parrish 2021 is a different assay\n",
    "(pgPEN dual-KO screen) in different cell lines (PC9/HeLa), not DepMap-derived.\n",
    "Compare these numbers to Section 5 (original) and the paralogy-identity-fix\n",
    "result (Precision=0.308, Recall=0.250) -- consistency across an independent\n",
    "dataset is real evidence the fix generalizes, not just overfitting.\n")












