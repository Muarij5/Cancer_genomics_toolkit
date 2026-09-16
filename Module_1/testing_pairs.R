# ============================================================
# post_fix_check.R
# Run AFTER Module2_finalized.R (with your edited line) is sourced
# and APP_DATA is loaded. Uses current/original thresholds.
# ============================================================

COHEN_D_THRESHOLD <<- 0.35
PARTIAL_COHEN_D_THRESHOLD <<- 0.25
CORRELATION_THRESHOLD <<- 0.25

full <- read.csv("D:/Tool_genomic/depmap_data/benchmark_results_full.csv", stringsAsFactors = FALSE)
full <- full[full$fr_label != "ambiguous" & is.na(full$error), ]

# ------------------------------------------------------------
# PART 1 -- did the fix break any previously-correct TP pair?
# ------------------------------------------------------------
tp_rows <- full[full$fr_label == "FR_positive" &
                  full$fr_call %in% c("Functional Redundancy", "Functional Redundancy (moderate)"), ]

cat("=== PART 1: Re-checking", nrow(tp_rows), "previously-correct TP pairs ===\n\n")
tp_still_ok <- 0
for (i in seq_len(nrow(tp_rows))) {
  ga <- tp_rows$gene_a[i]; gb <- tp_rows$gene_b[i]
  r <- tryCatch(
    analyze_gene_pair_v2(ga, gb, APP_DATA$dependency_df, APP_DATA$mutation_df, APP_DATA$expression_df,
                         APP_DATA$paralog_df, model_df = APP_DATA$model_df, cn_df = APP_DATA$cn_df),
    error = function(e) NULL
  )
  new_call <- if (is.null(r)) "ERROR" else r$FR$call
  still_positive <- new_call %in% c("Functional Redundancy", "Functional Redundancy (moderate)")
  if (still_positive) tp_still_ok <- tp_still_ok + 1
  cat(sprintf("%-10s/%-10s | old_call=%-35s | NEW_call=%-35s | still_positive=%s\n",
              ga, gb, tp_rows$fr_call[i], new_call, still_positive))
}
cat(sprintf("\nTPs retained: %d / %d\n\n", tp_still_ok, nrow(tp_rows)))

# ------------------------------------------------------------
# PART 2 -- full FR benchmark on the ~80-pair subset (fast check)
# ------------------------------------------------------------
wrong_or_positive <- full[full$fr_label == "FR_positive" |
                            full$fr_call %in% c("Functional Redundancy",
                                                "Functional Redundancy (moderate)",
                                                "Borderline FR"), ]
set.seed(42)
tn_pool <- full[full$fr_label == "FR_negative" &
                  full$fr_call == "Not Functionally Redundant", ]
tn_sample <- tn_pool[sample(nrow(tn_pool), min(60, nrow(tn_pool))), ]
bench <- unique(rbind(wrong_or_positive, tn_sample)[, c("gene_a", "gene_b", "fr_label")])

cat("=== PART 2: Full re-score on", nrow(bench), "pairs (subset) ===\n")
calls <- character(nrow(bench))
for (i in seq_len(nrow(bench))) {
  res <- tryCatch(
    analyze_gene_pair_v2(bench$gene_a[i], bench$gene_b[i],
                         APP_DATA$dependency_df, APP_DATA$mutation_df, APP_DATA$expression_df,
                         APP_DATA$paralog_df, compounds_df = APP_DATA$compounds_df,
                         model_df = APP_DATA$model_df, cn_df = APP_DATA$cn_df, cgc_df = APP_DATA$cgc_df),
    error = function(e) NULL
  )
  calls[i] <- if (is.null(res)) NA else res$FR$call
}

pred_pos  <- calls %in% c("Functional Redundancy", "Functional Redundancy (moderate)")
truth_pos <- bench$fr_label == "FR_positive"
ok <- !is.na(calls)
tp <- sum(pred_pos[ok] & truth_pos[ok]); fp <- sum(pred_pos[ok] & !truth_pos[ok]); fn <- sum(!pred_pos[ok] & truth_pos[ok])

cat(sprintf("\nAFTER FIX -- TP=%d FP=%d FN=%d\n", tp, fp, fn))
cat(sprintf("Precision=%.3f  Recall=%.3f\n", tp/(tp+fp), tp/(tp+fn)))
cat("\n(Compare to BEFORE: TP=5 FP=33 FN=11, Precision=0.132, Recall=0.312)\n")











full <- read.csv("D:/Tool_genomic/depmap_data/benchmark_results_full.csv", stringsAsFactors = FALSE)
full <- full[full$fr_label != "ambiguous" & is.na(full$error), ]
tp_rows <- full[full$fr_label == "FR_positive" &
                  full$fr_call %in% c("Functional Redundancy", "Functional Redundancy (moderate)"), ]

for (i in seq_len(nrow(tp_rows))) {
  ga <- tp_rows$gene_a[i]; gb <- tp_rows$gene_b[i]
  r <- tryCatch(
    analyze_gene_pair_v2(ga, gb, APP_DATA$dependency_df, APP_DATA$mutation_df, APP_DATA$expression_df,
                         APP_DATA$paralog_df, model_df = APP_DATA$model_df, cn_df = APP_DATA$cn_df),
    error = function(e) NULL
  )
  new_call <- if (is.null(r)) "ERROR" else r$FR$call
  cat(sprintf("%-10s/%-10s | NEW_call=%s\n", ga, gb, new_call))
}