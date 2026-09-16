# ============================================================
# Benchmark harness for run_fr_pipeline() / run_sl_pipeline()
# Ground truth: benchmark_pairs.csv
#   - fr_label:         FR_positive / ambiguous / FR_negative   (Dede & Hart 2020 zdLFC, wet-lab dual KO)
#   - sl_label_dekegel: True / False / not_in_mmc6              (De Kegel 2021, DepMap-derived)
#
# Assumes APP_DATA is already loaded in this session (the same list
# object your app builds via load_all_real_data() / make_mock_data()),
# with elements APP_DATA$dependency_df, APP_DATA$mutation_df,
# APP_DATA$expression_df, APP_DATA$paralog_df, APP_DATA$model_df,
# APP_DATA$cn_df, APP_DATA$cgc_df, APP_DATA$compounds_df — exactly as
# used by analyze_gene_pair_v2().
# ============================================================
DATA_DIR <- "D:/Tool_genomic/depmap_data"
bench <- read.csv(file.path(DATA_DIR, "benchmark_pairs.csv"), stringsAsFactors = FALSE)

results <- vector("list", nrow(bench))

for (i in seq_len(nrow(bench))) {
  ga <- bench$gene_a[i]
  gb <- bench$gene_b[i]
  
  res <- tryCatch(
    analyze_gene_pair_v2(ga, gb, APP_DATA$dependency_df, APP_DATA$mutation_df, APP_DATA$expression_df,
                         APP_DATA$paralog_df, compounds_df = APP_DATA$compounds_df,
                         model_df = APP_DATA$model_df, cn_df = APP_DATA$cn_df, cgc_df = APP_DATA$cgc_df),
    error = function(e) list(error = conditionMessage(e))
  )
  
  if (!is.null(res$error)) {
    results[[i]] <- data.frame(gene_a = ga, gene_b = gb,
                               fr_call = NA, sl_call = NA, error = res$error)
  } else {
    results[[i]] <- data.frame(gene_a = ga, gene_b = gb,
                               fr_call = res$FR$call, sl_call = res$SL$call,
                               error = NA)
  }
}

out <- do.call(rbind, results)
out <- merge(bench, out, by = c("gene_a", "gene_b"))

n_errored <- sum(!is.na(out$error))
if (n_errored > 0) {
  message(sprintf("%d / %d pairs failed to run (missing genes / coverage gaps) — excluded from scoring.",
                  n_errored, nrow(out)))
}
scored <- out[is.na(out$error), ]

# ------------------------------------------------------------
# FR scoring — against Dede/Hart wet-lab labels
# ------------------------------------------------------------
fr_eval <- scored[scored$fr_label != "ambiguous", ]  # drop 1-cell-line-only hits, too weak to score against
fr_eval$fr_true_positive <- fr_eval$fr_label == "FR_positive"

cat("\n=== FR: full contingency table (predicted call x ground truth) ===\n")
print(table(Predicted = fr_eval$fr_call, Truth = fr_eval$fr_label))

score_fr <- function(pred_positive_calls, label) {
  pred_pos <- fr_eval$fr_call %in% pred_positive_calls
  truth_pos <- fr_eval$fr_true_positive
  tp <- sum(pred_pos & truth_pos); fp <- sum(pred_pos & !truth_pos)
  fn <- sum(!pred_pos & truth_pos); tn <- sum(!pred_pos & !truth_pos)
  cat(sprintf("\n-- %s (positive = %s) --\n", label, paste(pred_positive_calls, collapse = " | ")))
  cat(sprintf("TP=%d FP=%d FN=%d TN=%d\n", tp, fp, fn, tn))
  cat(sprintf("Precision=%.3f  Recall(Sensitivity)=%.3f  Specificity=%.3f\n",
              tp / (tp + fp), tp / (tp + fn), tn / (tn + fp)))
}

# Strict: only real FR calls count as positive (Borderline FR treated as negative)
score_fr(c("Functional Redundancy", "Functional Redundancy (moderate)"), "STRICT")
# Lenient: Borderline FR also counted as positive (shows cost of the loose threshold)
score_fr(c("Functional Redundancy", "Functional Redundancy (moderate)", "Borderline FR"), "LENIENT (incl. Borderline FR)")

# ------------------------------------------------------------
# SL scoring — against De Kegel DepMap-derived labels (secondary/less independent, see caveat below)
# ------------------------------------------------------------
sl_eval <- scored[toupper(trimws(scored$sl_label_dekegel)) %in% c("TRUE", "FALSE"), ]
sl_eval$sl_true_positive <- toupper(trimws(sl_eval$sl_label_dekegel)) == "TRUE"

cat("\n\n=== SL: full contingency table (predicted call x ground truth) ===\n")
print(table(Predicted = sl_eval$sl_call, Truth = sl_eval$sl_label_dekegel))

score_sl <- function(pred_positive_calls, label) {
  pred_pos <- sl_eval$sl_call %in% pred_positive_calls
  truth_pos <- sl_eval$sl_true_positive
  tp <- sum(pred_pos & truth_pos); fp <- sum(pred_pos & !truth_pos)
  fn <- sum(!pred_pos & truth_pos); tn <- sum(!pred_pos & !truth_pos)
  cat(sprintf("\n-- %s (positive = %s) --\n", label, paste(pred_positive_calls, collapse = " | ")))
  cat(sprintf("TP=%d FP=%d FN=%d TN=%d\n", tp, fp, fn, tn))
  cat(sprintf("Precision=%.3f  Recall(Sensitivity)=%.3f  Specificity=%.3f\n",
              tp / (tp + fp), tp / (tp + fn), tn / (tn + fp)))
}

score_sl(c("Synthetic Lethal", "Synthetic Lethal (moderate)"), "STRICT")
score_sl(c("Synthetic Lethal", "Synthetic Lethal (moderate)", "Borderline SL"), "LENIENT (incl. Borderline SL)")

write.csv(out, "benchmark_results_full.csv", row.names = FALSE)
cat("\nFull per-pair results written to benchmark_results_full.csv\n")
cat("NOTE: sl_eval draws on De Kegel's DepMap-derived SL labels — same underlying data type as your\n",
    "pipeline's own tests, so treat SL precision/recall here as a secondary check, not independent\n",
    "validation. The FR scoring above (Dede/Hart) is the wet-lab-grounded number that matters most.\n")