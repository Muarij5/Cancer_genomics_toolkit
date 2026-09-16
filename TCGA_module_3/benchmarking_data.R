# ============================================================
# Benchmark harness for analyze_gene_pair() (gene_pair_terminal_app_v5.R)
# Ground truth: Suppl_Table.xlsx
#   Sheets: "GBM (HRN1)", "OVCA (HRN1)", "GBM (HRN2)", "OVCA (HRN2)"
#   Each row = a published mutual-exclusivity MODULE (3-5 genes) that was
#   found to be significantly co-altered in a mutually-exclusive pattern
#   in TCGA GBM / OVCA cohorts (classic MEMo-style modules, e.g.
#   RB1/CDK4/CDKN2A in GBM cell-cycle pathway).
#
# ------------------------------------------------------------
# WHAT THIS BENCHMARK CAN AND CANNOT DO -- READ BEFORE TRUSTING NUMBERS
# ------------------------------------------------------------
# 1) GROUND TRUTH IS POSITIVE-ONLY.
#    The sheet lists modules that WERE found significant. It does not list
#    gene pairs that were tested and found NOT mutually exclusive. So this
#    script can only measure SENSITIVITY/RECALL for Term 3
#    (mutual_exclusivity): "of pairs known from the literature to be
#    mutually exclusive, how many does the tool also call MUTUALLY
#    EXCLUSIVE?" It CANNOT compute precision, specificity, or a full
#    contingency table -- that needs a negative set, which this file
#    does not provide. Do not report specificity/precision from this
#    script; they are not in this data.
#
# 2) MODULES -> PAIRS.
#    analyze_gene_pair() takes exactly 2 genes. Modules here have 3-5
#    genes. Each module is expanded into all pairwise combinations
#    (choose 2) and every pair inherits the module's ground-truth label.
#    This means: if a module of 3 genes is a real 3-way exclusivity
#    pattern but any ONE pairwise sub-relationship is weak (e.g. two of
#    the three genes rarely co-occur with each other specifically),
#    testing that pair alone may legitimately fail even though the module
#    is real. That is an expected/known limitation of pairwise
#    decomposition, not a benchmark bug -- flag it in your writeup,
#    don't hide it.
#
# 3) STUDY ID MAPPING IS AN ASSUMPTION, NOT FACT.
#    You told me you don't know the exact cBioPortal study these modules
#    came from. I defaulted to the current TCGA PanCancer Atlas studies
#    (gbm_tcga_pan_can_atlas_2018, ov_tcga_pan_can_atlas_2018) because
#    they are the most complete/current GBM and OVCA studies on
#    cBioPortal. This is almost certainly NOT the exact cohort/version
#    the original module-finding paper used (module significance was
#    computed on an older TCGA freeze). Practical effect: percent-altered
#    numbers in your Excel sheet may not exactly match what the pipeline
#    recomputes now, and some genes' mutation/CNA rates will differ from
#    the original paper. The mutual-exclusivity CALL (direction) is
#    usually more robust to this than exact p-values, but you should
#    treat this as an approximate re-validation, not a bit-for-bit
#    reproduction. CHANGE STUDY_ID_MAP BELOW if you know the real IDs.
#
# 4) HRN1 vs HRN2 sheets are two different module-finding runs on the
#    same two cancer types. They are pooled together per cancer type
#    below; duplicate pairs (same 2 genes appearing in both) are
#    deduplicated, keeping the first occurrence.
#
# 5) OTHER 5 TERMS (subtype specificity, combination signal, double hit,
#    cell-cycle, angiogenesis) have NO ground truth in this file. This
#    script logs their verdicts per pair for your manual inspection but
#    does NOT score them right/wrong. Don't compute accuracy on them.
#
# ------------------------------------------------------------
# FIX LOG (this version vs. the one that crashed on distinct()):
#   The original parse_sheet() assumed FIXED column positions:
#       col1 = blank, col2 = Module ID, col3 = Genes, col4 = % altered,
#       col5 = p-value, col6 = q/pstar
#   If your sheet's real layout doesn't match that exact position/offset
#   (e.g. no leading blank column, or columns in a different order), every
#   row silently gets the WRONG data in the "genes" field. Every module
#   then fails to split into >=2 genes in expand_module_to_pairs(), which
#   returns an empty 0-column tibble() for every row. bind_rows() of all
#   those empties produces a table with literally no gene_a/gene_b/
#   cancer_type columns -- which is why distinct(gene_a, gene_b,
#   cancer_type) blew up with "Must use existing variables."
#
#   FIX: parse_sheet() now finds each needed column by matching the
#   HEADER TEXT itself ("Module ID", "Genes", etc.) instead of assuming a
#   fixed position. It also prints, per sheet, exactly which column index
#   got matched to which field -- so if your file still doesn't parse,
#   the printed mapping tells you immediately why, instead of an opaque
#   dplyr error three steps downstream. A hard stop() with a clear
#   message now fires immediately if Module ID / Genes columns can't be
#   found, and again if zero pairs could be expanded, rather than letting
#   the pipeline limp forward into a cryptic crash.
# ------------------------------------------------------------
#
# Usage:
#   1. source() your gene_pair_terminal_app_v5.R functions WITHOUT letting
#      run_terminal_app() auto-execute -- comment out the last line
#      `run_terminal_app()` in that file, or wrap this benchmark's source
#      call so the interactive prompt loop doesn't block. Simplest fix:
#      duplicate the file, delete/comment the final `run_terminal_app()`
#      call, save as gene_pair_terminal_app_v5_funcs.R, and source that.
#   2. Update PIPELINE_FILE and SUPPL_TABLE_PATH below.
#   3. Rscript run_benchmark_gene_pair_v5_FIXED.R
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
})

# ---- 0. Config -- CHANGE THESE ----
PIPELINE_FILE     <- Sys.getenv("TCGA_PIPELINE_FILE", unset = "Finalized_module3.R")
SUPPL_TABLE_PATH  <- Sys.getenv("TCGA_SUPPL_TABLE", unset = "Suppl.Table.xlsx")
OUT_CSV           <- Sys.getenv("TCGA_OUT_CSV", unset = "benchmark_results_mutual_exclusivity.csv")

# ASSUMPTION (see note 3 above) -- edit if you know the real study IDs
STUDY_ID_MAP <- c(
  "GBM"  = "gbm_tcga_pub",
  "OVCA" = "ov_tcga_pub"
)

# ---- 0b. Source the pipeline file WITHOUT letting it auto-launch the
#          interactive terminal app. The original v5 file ends with a bare
#          `run_terminal_app()` call, which -- if sourced as-is -- starts
#          reading readline()/stdin prompts and silently swallows every
#          line below this point in THIS script as if it were typed into
#          those prompts (that's what happened last run). To avoid needing
#          you to hand-edit the pipeline file, we read it as text, drop
#          only the trailing auto-run call, and source the rest from a
#          temp file.
source(PIPELINE_FILE)

pipeline_lines[is_autorun_call] <- "# [benchmark harness removed auto-run call here]"
tmp_pipeline <- tempfile(fileext = ".R")
writeLines(pipeline_lines, tmp_pipeline)
source(tmp_pipeline)
cat(sprintf("Sourced %d functions/definitions from %s (auto-run line stripped).\n",
            sum(grepl("<- function", pipeline_lines)), PIPELINE_FILE))

# ---- 1. Parse Suppl_Table.xlsx into a pair-level ground-truth table ----
sheet_names <- excel_sheets(SUPPL_TABLE_PATH)
# sheet name format: "GBM (HRN1)", "OVCA (HRN2)" etc.

parse_sheet <- function(sheet) {
  raw <- suppressMessages(read_excel(SUPPL_TABLE_PATH, sheet = sheet, col_names = FALSE))
  raw <- as.data.frame(raw, stringsAsFactors = FALSE)
  
  # ---- locate header row (any row containing "Module ID" in any cell) ----
  header_row <- which(apply(raw, 1, function(r) any(grepl("module\\s*id", r, ignore.case = TRUE))))[1]
  if (is.na(header_row)) {
    warning(sprintf("Sheet '%s': no 'Module ID' header found anywhere -- skipping this sheet.", sheet))
    return(tibble())
  }
  header_vals <- as.character(raw[header_row, ])
  
  # ---- locate each needed column BY HEADER TEXT, not by fixed position.
  #      This is the actual fix: the old code assumed col2/col3/col4/... in
  #      a specific order and broke if the sheet's real layout differed. ----
  find_col <- function(pattern) {
    idx <- which(sapply(header_vals, function(x) !is.na(x) && grepl(pattern, x, ignore.case = TRUE)))[1]
    if (length(idx) == 0) NA_integer_ else idx
  }
  col_module <- find_col("module\\s*id")
  col_genes  <- find_col("gene")
  col_pct    <- find_col("%|alter")
  col_pval   <- find_col("^p[\\s._-]*val|p-?value")
  col_q      <- find_col("^q[\\s._-]*val|fdr|pstar|q-?value")
  
  missing <- c(module_id = is.na(col_module), genes = is.na(col_genes))
  if (any(missing)) {
    stop(sprintf(
      paste0("Sheet '%s': could not locate required column(s) [%s] by header text.\n",
             "  Header row %d actually contains: [%s]\n",
             "  Fix: either rename the header cells to include the words 'Module ID' / 'Genes',\n",
             "  or edit find_col() patterns above to match your actual header wording."),
      sheet, paste(names(missing)[missing], collapse = ", "),
      header_row, paste(header_vals, collapse = " | ")))
  }
  
  df <- raw[(header_row + 1):nrow(raw), , drop = FALSE]
  out <- tibble(
    module_id   = as.character(df[[col_module]]),
    genes       = as.character(df[[col_genes]]),
    pct_altered = if (!is.na(col_pct))  as.character(df[[col_pct]])  else NA_character_,
    p_value     = if (!is.na(col_pval)) as.character(df[[col_pval]]) else NA_character_,
    q_or_pstar  = if (!is.na(col_q))    as.character(df[[col_q]])    else NA_character_
  ) %>%
    filter(!is.na(module_id), module_id != "", !is.na(genes), genes != "")
  
  out$cancer_type     <- str_extract(sheet, "^[A-Za-z]+")
  out$network_version <- str_extract(sheet, "HRN\\d")
  out <- out %>%
    mutate(sig_value = suppressWarnings(as.numeric(gsub("<", "", q_or_pstar)))) %>%
    filter(!is.na(sig_value) & sig_value < 0.05) %>%
    select(-sig_value)
  message(sprintf(
    "Sheet '%s': header row %d -> module_id=col%d, genes=col%d, pct_altered=col%s, p_value=col%s, q=col%s | %d data rows kept",
    sheet, header_row, col_module, col_genes,
    ifelse(is.na(col_pct), "NA", col_pct),
    ifelse(is.na(col_pval), "NA", col_pval),
    ifelse(is.na(col_q), "NA", col_q),
    nrow(out)))
  
  out
}

modules <- map_dfr(sheet_names, parse_sheet)

if (nrow(modules) == 0) {
  stop("No modules were parsed from ANY sheet in ", SUPPL_TABLE_PATH,
       " -- check the per-sheet messages/warnings printed above for why.")
}

# ---- 2. Expand each module into all pairwise gene combinations ----
expand_module_to_pairs <- function(row) {
  genes <- str_split(row$genes, ",\\s*")[[1]] %>% str_trim()
  genes <- genes[nzchar(genes)] %>% unique()
  if (length(genes) < 2) {
    warning(sprintf(
      "Module '%s' (%s, %s): only %d gene(s) parsed from genes field '%s' -- skipped (need >=2).",
      row$module_id, row$cancer_type, row$network_version, length(genes), row$genes))
    return(tibble())
  }
  combos <- combn(genes, 2, simplify = FALSE)
  tibble(
    gene_a = map_chr(combos, 1),
    gene_b = map_chr(combos, 2),
    module_id = row$module_id,
    cancer_type = row$cancer_type,
    network_version = row$network_version,
    module_pct_altered = row$pct_altered,
    module_p_value = as.character(row$p_value)
  )
}

# NOTE: map over row indices directly instead of split(seq_len(nrow(.))) --
# functionally equivalent, but avoids relying on split() preserving a
# data-frame class per chunk, and makes the empty-result guard below
# possible to check cleanly.
pairs_raw <- map_dfr(seq_len(nrow(modules)), ~ expand_module_to_pairs(modules[.x, ]))

if (nrow(pairs_raw) == 0 || !all(c("gene_a", "gene_b", "cancer_type") %in% names(pairs_raw))) {
  stop("No valid gene pairs could be expanded from ANY module (every module had < 2 parsed genes).\n",
       "  Check the warnings above -- they show the raw 'genes' field value the parser saw.\n",
       "  Most likely cause: the 'genes' column detected by parse_sheet() isn't actually the\n",
       "  comma-separated gene-symbol column in your sheet. Check the 'genes=colN' mapping\n",
       "  printed per sheet above against your spreadsheet.")
}

pairs <- pairs_raw %>%
  distinct(gene_a, gene_b, cancer_type, .keep_all = TRUE)  # dedup HRN1/HRN2 overlap per cancer type

pairs$study_id <- unname(STUDY_ID_MAP[pairs$cancer_type])
unmapped <- pairs %>% filter(is.na(study_id)) %>% distinct(cancer_type)
if (nrow(unmapped) > 0) {
  stop("No study_id mapped for cancer_type(s): ", paste(unmapped$cancer_type, collapse = ", "),
       ". Add them to STUDY_ID_MAP.")
}

cat(sprintf("Parsed %d modules -> %d unique (gene pair, cancer type) combinations to test.\n",
            nrow(modules), nrow(pairs)))

# ---- 3. Run the pipeline on every pair ----
results <- vector("list", nrow(pairs))

for (i in seq_len(nrow(pairs))) {
  ga <- pairs$gene_a[i]; gb <- pairs$gene_b[i]; sid <- pairs$study_id[i]
  cat(sprintf("[%d/%d] %s + %s (%s)\n", i, nrow(pairs), ga, gb, sid))
  
  res <- tryCatch(
    analyze_gene_pair(ga, gb, sid),
    error = function(e) list(error = conditionMessage(e))
  )
  
  if (!is.null(res$error)) {
    results[[i]] <- tibble(
      gene_a = ga, gene_b = gb, study_id = sid, error = res$error,
      mutual_exclusivity_verdict = NA, mutual_exclusivity_p = NA,
      subtype_verdict = NA, combination_verdict = NA, double_hit_verdict = NA,
      cell_cycle_verdict = NA, angiogenesis_verdict = NA
    )
  } else {
    tab <- summarize_results(res$main)
    get_v <- function(term_label) tab$verdict[tab$term == term_label][1]
    get_p <- function(term_label) tab$p_value[tab$term == term_label][1]
    results[[i]] <- tibble(
      gene_a = ga, gene_b = gb, study_id = sid, error = NA,
      mutual_exclusivity_verdict = get_v("Mutually exclusive (same pathway)"),
      mutual_exclusivity_p       = get_p("Mutually exclusive (same pathway)"),
      subtype_verdict     = get_v("Tumor subtype specificity"),
      combination_verdict = get_v("Biomarker combination signal"),
      double_hit_verdict  = get_v("Tumor suppressor double hit (per-gene, not pair)"),
      cell_cycle_verdict  = get_v("Cell-cycle deregulation"),
      angiogenesis_verdict = get_v("Angiogenesis promotion"),
      arm_flagged = isTRUE(res$arm_flag$geneA$flagged) || isTRUE(res$arm_flag$geneB$flagged)
    )
  }
}

out <- bind_rows(results)
# FIX: the same (gene_a, gene_b) pair can legitimately appear under more
# than one cancer_type (e.g. RB1+NF1 in both a GBM module and an OVCA
# module) -- pairs was correctly deduped on (gene_a, gene_b, cancer_type),
# but joining on (gene_a, gene_b) alone loses that distinction and matches
# each result row to EVERY cancer_type variant of that pair, producing a
# many-to-many join that silently duplicates/cross-wires rows. study_id is
# 1:1 with cancer_type (via STUDY_ID_MAP) and is present on both sides, so
# adding it to the join key resolves the ambiguity correctly.
out <- pairs %>% select(gene_a, gene_b, cancer_type, module_id, network_version, study_id) %>%
  right_join(out, by = c("gene_a", "gene_b", "study_id"))

n_errored <- sum(!is.na(out$error))
if (n_errored > 0) {
  message(sprintf("%d / %d pairs failed to run (missing genes / coverage gaps) -- excluded from scoring.",
                  n_errored, nrow(out)))
}
scored <- out %>% filter(is.na(error))

# ------------------------------------------------------------
# 4. SCORING -- Term 3 (mutual_exclusivity) ONLY.
#    Ground truth here is positive-only: every pair below is drawn from a
#    published significant ME module, so ground_truth_positive == TRUE for
#    all rows. This yields RECALL/SENSITIVITY only -- see note 1 at top.
# ------------------------------------------------------------
cat("\n=== Term 3 (mutual_exclusivity) verdict distribution on known ME-module pairs ===\n")
print(table(scored$mutual_exclusivity_verdict, useNA = "ifany"))

n_total <- nrow(scored)
n_hit   <- sum(scored$mutual_exclusivity_verdict == "MUTUALLY EXCLUSIVE", na.rm = TRUE)
n_na    <- sum(scored$mutual_exclusivity_verdict == "N/A (not testable)", na.rm = TRUE)

cat(sprintf("\nRECALL (of testable pairs): %d / %d = %.3f\n",
            n_hit, n_total - n_na, n_hit / max(n_total - n_na, 1)))
cat(sprintf("Untestable (insufficient mutation variation): %d / %d\n", n_na, n_total))
cat("NOTE: No precision/specificity figure is reported -- this file has no negative examples.\n")
cat("      A pair NOT called MUTUALLY EXCLUSIVE here is a likely miss, but some misses are\n")
cat("      expected from pairwise decomposition of >2-gene modules (see note 2 at top).\n")

cat("\n=== Breakdown by cancer type ===\n")
scored %>%
  group_by(cancer_type) %>%
  summarise(n = n(),
            n_hit = sum(mutual_exclusivity_verdict == "MUTUALLY EXCLUSIVE", na.rm = TRUE),
            recall = n_hit / n, .groups = "drop") %>%
  print()

write.csv(out, OUT_CSV, row.names = FALSE)
cat(sprintf("\nFull per-pair results (all 6 term verdicts, unscored beyond Term 3) written to %s\n", OUT_CSV))





study_list <- get_study_list()

study_list %>% filter(grepl("gbm|glioblastoma", studyId, ignore.case = TRUE) |
                        grepl("glioblastoma", name, ignore.case = TRUE)) %>%
  print(n = Inf)

study_list %>% filter(grepl("^ov_|ovarian", studyId, ignore.case = TRUE) |
                        grepl("ovarian", name, ignore.case = TRUE)) %>%
  print(n = Inf)

for (sid in c("gbm_tcga_pub", "ov_tcga_pub")) {
  cat("\n--- Testing:", sid, "---\n")
  res <- tryCatch({
    sl <- get_default_sample_list(sid)
    roster <- get_full_sample_roster(sl)
    cat(sprintf("OK -- sample list '%s' has %d samples\n", sl, length(roster)))
  }, error = function(e) cat("FAILED:", conditionMessage(e), "\n"))
}

test_result <- analyze_gene_pair("RB1", "CDK4", "gbm_tcga_pub")
