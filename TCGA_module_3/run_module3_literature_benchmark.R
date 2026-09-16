# =====================================================================
# Runs module3_literature_benchmark_pairs.csv against your Module 3 tool
# and reports, per row, whether the tool's actual verdict matches the
# expected literature-based direction.
#
# HOW TO USE:
#   1. Make sure Module3 (Finalized_module3_CORRECTED.R) is already
#      sourced in this R session (or set MODULE3_PATH below and this
#      script will source it for you).
#   2. Put module3_literature_benchmark_pairs.csv in your working
#      directory (or set BENCHMARK_CSV_PATH below).
#   3. Run: source("run_module3_literature_benchmark.R")
# =====================================================================

MODULE3_PATH      <- "Finalized_module3_CORRECTED.R"   # set to your actual path, or leave as-is if already sourced
BENCHMARK_CSV_PATH <- Sys.getenv("TCGA_BENCHMARK_CSV", unset = "module3_literature_benchmark_pairs.csv")

if (!exists("analyze_gene_pair")) {
  if (file.exists(MODULE3_PATH)) {
    message("Sourcing Module 3 from '", MODULE3_PATH, "'...")
    source(MODULE3_PATH)
  } else {
    stop("analyze_gene_pair() not found and MODULE3_PATH does not exist. ",
         "Source your Module 3 file first, or fix MODULE3_PATH above.")
  }
}

benchmark <- read.csv(BENCHMARK_CSV_PATH, stringsAsFactors = FALSE)

# Maps each CSV "term" code to the exact term label summarize_results() uses
TERM_LABEL_MAP <- c(
  Term1_SubtypeSpecificity = "Tumor subtype specificity",
  Term2_CombinationSignal  = "Biomarker combination signal",
  Term3_MutualExclusivity  = "Mutually exclusive (same pathway)",
  Term4_DoubleHit          = "Tumor suppressor double hit (per-gene, not pair)",
  Term5_CellCycleE2F       = "Cell-cycle deregulation",
  Term6_Angiogenesis       = "Angiogenesis promotion"
)

run_one <- function(row) {
  cat(sprintf("\n--- %s: %s + %s in %s ---\n", row$term, row$gene_a, row$gene_b, row$study_id))

  result <- tryCatch(
    analyze_gene_pair(row$gene_a, row$gene_b, row$study_id),
    error = function(e) e
  )

  if (inherits(result, "error")) {
    return(data.frame(
      term = row$term, gene_a = row$gene_a, gene_b = row$gene_b, study_id = row$study_id,
      expected_result = row$expected_result, actual_verdict = NA_character_,
      actual_p_value = NA_real_, actual_adj_p_value = NA_real_,
      run_status = paste("ERROR:", conditionMessage(result)),
      stringsAsFactors = FALSE
    ))
  }

  tab <- summarize_results(result$main)
  target_label <- TERM_LABEL_MAP[[row$term]]
  tab_row <- tab[tab$term == target_label, ]

  if (nrow(tab_row) == 0) {
    return(data.frame(
      term = row$term, gene_a = row$gene_a, gene_b = row$gene_b, study_id = row$study_id,
      expected_result = row$expected_result, actual_verdict = NA_character_,
      actual_p_value = NA_real_, actual_adj_p_value = NA_real_,
      run_status = "Term label not found in summarize_results() output -- check TERM_LABEL_MAP",
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    term = row$term, gene_a = row$gene_a, gene_b = row$gene_b, study_id = row$study_id,
    expected_result = row$expected_result,
    actual_verdict = tab_row$verdict[1],
    actual_p_value = tab_row$p_value[1],
    actual_adj_p_value = tab_row$adj_p_value[1],
    run_status = "OK",
    stringsAsFactors = FALSE
  )
}

results_list <- lapply(seq_len(nrow(benchmark)), function(i) run_one(benchmark[i, ]))
results_df <- do.call(rbind, results_list)

cat("\n\n===================== BENCHMARK SUMMARY =====================\n")
print(results_df, row.names = FALSE)

out_path <- "module3_literature_benchmark_RESULTS.csv"
write.csv(results_df, out_path, row.names = FALSE)
cat(sprintf("\nSaved full results to '%s'.\n", out_path))
cat("\nRead 'actual_verdict' against 'expected_result' row by row -- the tool doesn't\n")
cat("auto-score match/mismatch because 'expected_result' is free text (e.g. 'Significant,\n")
cat("mutually exclusive (OR<1)'), not a fixed label. Compare by eye, same as you did for\n")
cat("the VHL/PBRM1 and RB1/CDK4 spot-checks already in your validation summary.\n")


# Run this AFTER run_module3_literature_benchmark.R has already populated
# results_df in your session (or read the saved CSV back in).
#
# It just prints the FULL detail/interpretation text for the rows that
# looked like a miss, so you can tell bug vs. known-limitation apart.

if (!exists("results_df")) {
  results_df <- read.csv("module3_literature_benchmark_RESULTS.csv", stringsAsFactors = FALSE)
}

# Rows to inspect closely (edit this list based on your own eyeballing of
# actual_verdict vs expected_result -- these are the ones that looked like
# misses in the run you already did):
mismatch_rows <- c(
  "GATA3-GAPDH", "FGFR3-GAPDH",       # Term 1
  "KRAS-EGFR", "TP53-MDM2",           # Term 3
  "CCNE1-GAPDH", "RB1-GAPDH-gbm",     # Term 5 (RB1 alone, gbm_tcga_pub)
  "VHL-GAPDH-kirc-angio"              # Term 6
)

# Re-run just these specific pairs and print the FULL detail text
pairs_to_check <- data.frame(
  gene_a   = c("GATA3", "FGFR3", "KRAS", "TP53", "CCNE1", "RB1", "VHL"),
  gene_b   = c("GAPDH", "GAPDH", "EGFR", "MDM2", "GAPDH", "GAPDH", "GAPDH"),
  study_id = c("brca_tcga_pub", "blca_tcga_pan_can_atlas_2018",
               "luad_tcga_pan_can_atlas_2018", "sarc_tcga_pan_can_atlas_2018",
               "ov_tcga_pan_can_atlas_2018", "gbm_tcga_pub",
               "kirc_tcga_pan_can_atlas_2018"),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(pairs_to_check))) {
  ga <- pairs_to_check$gene_a[i]; gb <- pairs_to_check$gene_b[i]; sid <- pairs_to_check$study_id[i]
  cat(sprintf("\n\n========== %s + %s (%s) ==========\n", ga, gb, sid))
  r <- tryCatch(analyze_gene_pair(ga, gb, sid), error = function(e) e)
  if (inherits(r, "error")) { cat("ERROR:", conditionMessage(r), "\n"); next }
  tab <- summarize_results(r$main)
  print(tab[, c("term", "verdict", "p_value", "detail")], row.names = FALSE, width = Inf)
}