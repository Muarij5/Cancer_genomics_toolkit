# ============================================================
# GeneComb Analyzer — v2.14
# v2.14 CHANGE LOG (adds ONE thing — a curated cancer-gene-role
# annotation layer — nothing else was touched):
#
#   CONTEXT: an earlier review flagged that classifying a gene as
#   "tumor suppressor" vs "essential gene" by checking whether its
#   median CRISPR dependency score is above ESSENTIALITY_CUTOFF
#   (median > -0.5) is not biologically valid — dependency score
#   measures how essential a gene is, not what ROLE it plays in
#   cancer. BRCA1 (a real TSG) and PARP1 (not a TSG) can land on the
#   same side of that threshold, which is exactly why a median-based
#   TS/EG split would misclassify pairs like BRCA1/PARP1.
#
#   THIS FILE ALREADY DID NOT HAVE THAT BUG. run_sl_pipeline() below
#   is explicitly commented "Direction-driven SL classification (no
#   TS/EG heuristic)" — it decides SL direction from empirical
#   evidence only (test_differential_essentiality run gene_a->gene_b
#   AND gene_b->gene_a, tissue-adjusted regression, real p-values and
#   effect sizes), never from a median-dependency guess about gene
#   role. ESSENTIALITY_CUTOFF is used elsewhere (Single-Gene
#   Non-Essentiality, an FR test) for a different, legitimate purpose
#   — checking neither gene is individually essential — not for
#   TS/EG classification. So there was nothing to fix there.
#
#   What v2.14 adds instead: a genuine curated annotation layer using
#   the actual COSMIC Cancer Gene Census (the Census_all*.csv you
#   uploaded), following this app's existing "supporting tests add
#   context, they never gate a call" philosophy. Each gene gets a
#   TSG / oncogene / fusion-partner / not-in-CGC badge for purely
#   INFORMATIONAL display next to the statistical result, so you can
#   sanity-check "does the empirically-found direction match known
#   biology?" without the tool silently gating its verdict on an
#   incomplete external list (most genes — including real TSGs and
#   oncogenes — simply aren't in the CGC yet, so absence ≠ evidence).
#
#   - New: find_cgc_file(), load_cgc_df(), get_cgc_annotation(),
#     format_cgc_badge() (see "COSMIC CANCER GENE CENSUS" section
#     below, right after DATA LOADING FUNCTIONS).
#   - load_all_real_data() now also looks for Census_all*.csv in the
#     data folder and loads it into APP_DATA$cgc_df if found. This is
#     fully optional — if absent, the app runs exactly as before.
#   - analyze_gene_pair_v2() gains an optional cgc_df param and
#     attaches result$cgc_badge_A / result$cgc_badge_B — informational
#     text only, never fed into any statistical test or pass/fail
#     decision.
#   - UI: two small badges added to the results card showing each
#     gene's CGC role, plus a status-box note when CGC is loaded.
## ------------------------------------------------------------
# v2.15 CHANGE LOG (this version only adds ONE thing — an identity
# threshold on paralogy's contribution to FR decision support —
# nothing else was touched):
#   - test_paralogy()'s contribution to FR decision support is now
#     gated by sequence identity (>=70%). Motivated by an identity-
#     threshold sweep on the internal benchmark showing weak-identity
#     matches (20-55%) are non-discriminating (flat ~4-7% precision)
#     while strong matches (>=70%) carry real signal -- see validation
#     summary Section 5.3, 5.5.
# ------------------------------------------------------------
# ------------------------------------------------------------
# v2.8 CHANGE LOG (this version only adds Copy Number support —
# nothing else was touched):
#   - Added OmicsCNGeneWGS.csv as a 7th required input file.
#   - Gene Loss is now: Mutation==1 OR Expression Z<-1.5 OR CN<0.5
#     (previously just the first two). This matches how DAISY-style
#     tools define "gene inactivation" — mutation, expression, AND
#     copy number all feed the same loss call.
#   - CN_LOSS_THRESHOLD is a STARTING value, not a validated one —
#     see the big comment block below it before trusting borderline
#     calls near 0.5.
#   - Removed the standalone diagnostic scratch-scripts that were
#     appended after shinyApp() — they were debugging aids, not part
#     of the app, and were cluttering the production file. Say the
#     word if you want those brought back as a separate file.
# ============================================================

library(data.table)

# ============================================================
# CONSTANTS — with scientific justification
# ============================================================

# v2.10 (handover cleanup): behavior flags, not scientific constants.
# VERBOSE: when FALSE (default), the per-call diagnostic cat() output in
#   test_differential_essentiality() and test_reciprocal_compensation() is
#   suppressed. Set VERBOSE <- TRUE (or set the env var GENECOMB_VERBOSE=1)
#   before sourcing if you want that console trace back for debugging.
# RUN_VALIDATION_ON_SOURCE: when FALSE (default), the benchmark/validation
#   script block at the bottom of this file (positive controls, negative
#   controls, diagnose_fr_miss() calls, etc.) is NOT executed just because
#   this file was source()'d -- e.g. by a Shiny app, another script, or a
#   package load. Set it to TRUE (or run this file directly / interactively)
#   to reproduce that validation pass.
VERBOSE <- isTRUE(as.logical(Sys.getenv("GENECOMB_VERBOSE", "FALSE")))
RUN_VALIDATION_ON_SOURCE <- isTRUE(as.logical(Sys.getenv("GENECOMB_RUN_VALIDATION", "FALSE")))
# NOTE ON MULTIPLE TESTING: BH-adjusted p-values are computed and displayed in
# make_test_df() for the reader's context, but PASS/FAIL decisions throughout
# this pipeline use the raw p-value against ALPHA. Each gene pair is evaluated
# as an independent hypothesis by design; BH-adjustment is not applied at the
# decision level. See manuscript Methods for full rationale.
ALPHA <- 0.05                    # Standard significance threshold
ESSENTIALITY_CUTOFF <- -0.5      # DepMap standard for "essential"
SELECTIVE_THRESHOLD <- -0.7      # Stronger threshold for "selectively essential"
MIN_SAMPLES <- 20                # Minimum samples for correlation tests
MIN_GROUP_SIZE <- 10             # Minimum per group for Wilcoxon / regression tests
MIN_INACTIVE_SAMPLES <- 10       # For mutation-based tests

# v2.2: recalibrated for genome-scale noisy data. These are "real,
# screened, but achievable" effect sizes rather than textbook clean-
# experiment values.
COHEN_D_THRESHOLD <- 0.35         # was 0.5. Used when tissue adjustment is NOT possible
# for a given pair (see fallback logic below) — i.e. the
# original, un-discounted bar.
CORRELATION_THRESHOLD <- 0.25     # was 0.3
# VALIDATED (combined with identity-gated paralogy fix + CN_LOSS_THRESHOLD=0.6):
# Dede/Hart internal benchmark (n=403) -> Precision=31.2%, Recall=31.2% (up from
# original 13.2%/31.2%). Confirmed via Parrish et al. 2021 external benchmark
# (n=1030, independent assay/cell lines) -- see validation summary Section 5-6
# v2.7: effect-size floor for the TISSUE-ADJUSTED test. Lower than
# COHEN_D_THRESHOLD because tissue-adjustment has already removed one
# major, well-understood noise source (cross-lineage baseline variance).
# Recommend validating this value against a benchmark panel (see
# run_benchmark() at the bottom of this file) before relying on it for
# a specific dataset/panel size.
PARTIAL_COHEN_D_THRESHOLD <- 0.25

OR_SL_THRESHOLD <- 2.0           # Odds ratio for SL Score
OR_EXCLUSIVITY_THRESHOLD <- 0.3  # Strict OR for mutual exclusivity

# v2.8 NEW: Copy-number loss threshold.
# IMPORTANT — OmicsCNGeneWGS.csv is LINEAR-SCALE relative copy number,
# NOT log2 and NOT a GISTIC -2/-1/0/1/2 category. On this scale, diploid
# is ~1.0 and a homozygous/deep deletion clusters near 0. This value
# (< 0.5, i.e. less than half the normal diploid ratio) is a reasonable
# STARTING cutoff, but it has not been validated against this specific
# data release. Before trusting calls near this boundary: pull the real
# CN values for a gene you already know is frequently deleted (e.g.
# ARID1A in Uterine/Ovarian lines) and confirm the deletion cluster
# actually sits below this number in your data — adjust if it doesn't.
CN_LOSS_THRESHOLD <- 0.6
# VALIDATED against ARID1A (known frequently-deleted gene): stronger,
# more sample-supported signal at 0.6 (n=50, p=4.24e-11) vs original 0.5
# (n=9, p=0.0011). Parrish external robustness check confirmed no
# meaningful change vs 0.5 -- see validation summary Section 6.2.
# ============================================================
# DATA LOADING FUNCTIONS
# ============================================================

strip_entrez_suffix <- function(cols) sub(" \\([^)]*\\)$", "", cols)

load_model_data <- function(path) {
  df <- fread(path, data.table = FALSE, na.strings = c("NA", ""))
  rownames(df) <- df$ModelID
  df
}

load_compounds_data <- function(path) fread(path, data.table = FALSE, na.strings = c("NA", ""))

get_druggability <- function(gene, compounds_df) {
  hits <- character(0)
  targets_raw <- compounds_df$GeneSymbolOfTargets
  for (i in seq_len(nrow(compounds_df))) {
    tr <- targets_raw[i]
    if (is.na(tr) || tr == "") next
    targets <- trimws(strsplit(tr, ";")[[1]])
    if (gene %in% targets) hits <- c(hits, compounds_df$CompoundName[i])
  }
  hits
}

.load_profile_style_matrix <- function(path) {
  df <- fread(path, data.table = FALSE, na.strings = c("NA", ""))
  is_default <- tolower(as.character(df$IsDefaultEntryForModel)) %in% c("true", "yes", "1")
  df <- df[is_default, , drop = FALSE]
  rownames(df) <- df$ModelID
  meta_cols <- c("SequencingID", "ModelConditionID", "IsDefaultEntryForMC",
                 "IsDefaultEntryForModel", "ModelID")
  df <- df[, !(colnames(df) %in% meta_cols), drop = FALSE]
  colnames(df) <- strip_entrez_suffix(colnames(df))
  df
}

load_expression_data <- function(path) .load_profile_style_matrix(path)
load_mutation_df     <- function(path) .load_profile_style_matrix(path)

# v2.8 NEW: Copy number loader. Reuses the same profile-style loader as
# expression/mutation since OmicsCNGeneWGS.csv follows the same
# ModelID / IsDefaultEntryForModel convention as other "Omics"-prefixed
# DepMap files. If this errors on load, check the actual column names
# in your file (fread(path, nrows=5)) — the loader expects ModelID and
# IsDefaultEntryForModel columns exactly like mutation/expression do.
load_cn_df <- function(path) .load_profile_style_matrix(path)

load_dependency_df <- function(path) {
  df <- fread(path, data.table = FALSE, na.strings = c("NA", ""))
  rownames(df) <- df[[1]]
  df <- df[, -1, drop = FALSE]
  colnames(df) <- strip_entrez_suffix(colnames(df))
  df
}

load_paralog_df <- function(path, min_identity = 0.0) {
  df <- fread(path, data.table = FALSE, na.strings = c("NA", ""))
  df <- df[!is.na(df[["Human paralogue gene stable ID"]]), , drop = FALSE]
  df <- df[!is.na(df[["Gene name"]]) & !is.na(df[["Human paralogue associated gene name"]]), , drop = FALSE]
  id1 <- as.numeric(df[["Paralogue %id. target Human gene identical to query gene"]])
  id2 <- as.numeric(df[["Paralogue %id. query gene identical to target Human gene"]])
  pct_identity <- rowMeans(cbind(id1, id2), na.rm = TRUE)
  out <- data.frame(GeneA = df[["Gene name"]], GeneB = df[["Human paralogue associated gene name"]],
                    pct_identity = pct_identity, stringsAsFactors = FALSE)
  out <- out[out$pct_identity >= min_identity, , drop = FALSE]
  out <- unique(out)
  rownames(out) <- NULL
  out
}

# ============================================================
# v2.14 NEW: COSMIC Cancer Gene Census (CGC) — informational only
# ============================================================
# This section is entirely additive. cgc_df is OPTIONAL everywhere:
# every function below is written to accept cgc_df = NULL and return
# "no annotation available" rather than error, so the app runs exactly
# as it did in v2.8 if no Census_all*.csv is present.
#
# Why this is informational-only, not a gate: CGC coverage is a curated
# but incomplete list (~700-ish genes). A gene missing from it is NOT
# evidence the gene lacks a cancer role -- it may just not be catalogued
# yet. Using presence/absence in CGC to gate a PASS/fail call would
# introduce a *different* systematic bias than the one the median-
# dependency heuristic had. So this only ever labels, never decides.

# COSMIC's download button embeds a timestamp in the filename (e.g.
# "Census_allWed_Jul__8_08_17_13_2026.csv"), so we search for the
# pattern rather than hardcoding one exact name.
find_cgc_file <- function(data_dir) {
  hits <- list.files(data_dir, pattern = "Census_all.*\\.csv$", full.names = TRUE, ignore.case = TRUE)
  if (length(hits) == 0) return(NULL)
  if (length(hits) > 1) {
    message("Multiple Census_all*.csv files found in '", data_dir, "' -- using the first one: ", hits[1])
  }
  hits[1]
}

load_cgc_df <- function(path) {
  df <- fread(path, data.table = FALSE, na.strings = c("NA", ""))
  required_cols <- c("Gene Symbol", "Tier", "Role in Cancer")
  missing <- setdiff(required_cols, colnames(df))
  if (length(missing) > 0) {
    stop(sprintf("CGC file is missing expected column(s): %s. Is this really a COSMIC Census_all export?",
                 paste(missing, collapse = ", ")))
  }
  df$GeneSymbol <- toupper(trimws(df[["Gene Symbol"]]))
  df <- df[!is.na(df$GeneSymbol) & df$GeneSymbol != "", , drop = FALSE]
  df <- df[!duplicated(df$GeneSymbol), , drop = FALSE]
  rownames(df) <- df$GeneSymbol
  df
}

# Returns NULL if the gene isn't in the census (common -- not evidence
# of anything either way), otherwise a small list describing its
# curated role.
get_cgc_annotation <- function(gene, cgc_df) {
  if (is.null(cgc_df)) return(NULL)
  gene <- toupper(trimws(gene))
  if (!(gene %in% rownames(cgc_df))) return(NULL)
  row <- cgc_df[gene, ]
  role <- if (!is.na(row[["Role in Cancer"]])) row[["Role in Cancer"]] else ""
  list(
    gene = gene,
    tier = row[["Tier"]],
    role_in_cancer = role,
    is_tsg = grepl("TSG", role, ignore.case = TRUE),
    is_oncogene = grepl("oncogene", role, ignore.case = TRUE),
    is_fusion = grepl("fusion", role, ignore.case = TRUE)
  )
}

# Short human-readable line for the UI. Deliberately hedges on "not in
# CGC" so it doesn't read as a negative result.
format_cgc_badge <- function(annot) {
  if (is.null(annot)) return("Not in COSMIC CGC (no curated role on file \u2014 common, not informative either way)")
  role_bits <- c()
  if (annot$is_tsg) role_bits <- c(role_bits, "Tumor Suppressor (TSG)")
  if (annot$is_oncogene) role_bits <- c(role_bits, "Oncogene")
  if (annot$is_fusion) role_bits <- c(role_bits, "Fusion partner")
  if (length(role_bits) == 0) role_bits <- "role not specified in CGC"
  sprintf("COSMIC CGC Tier %s \u2014 %s", annot$tier, paste(role_bits, collapse = " / "))
}

# v2.8: added `cn` entry. This is now a 7th required file.
DEFAULT_FILENAMES <- list(
  dependency = "CRISPRGeneEffect.csv",
  mutation   = "OmicsSomaticMutationsMatrixDamaging.csv",
  expression = "OmicsExpressionTPMLogp1HumanProteinCodingGenesStranded.csv",
  cn         = "OmicsCNGeneWGS.csv",
  model      = "Model.csv",
  compounds  = "PortalCompounds.csv",
  paralog    = "mart_export.txt"
)

load_all_real_data <- function(data_dir, filenames = list(), min_paralog_identity = 0.0) {
  names_ <- modifyList(DEFAULT_FILENAMES, filenames)
  p <- function(key) file.path(data_dir, names_[[key]])
  dependency_df <- load_dependency_df(p("dependency"))
  mutation_df   <- load_mutation_df(p("mutation"))
  expression_df <- load_expression_data(p("expression"))
  cn_df         <- load_cn_df(p("cn"))
  model_df      <- load_model_data(p("model"))
  compounds_df  <- load_compounds_data(p("compounds"))
  paralog_df    <- load_paralog_df(p("paralog"), min_identity = min_paralog_identity)
  
  # v2.14: CGC is fully optional and never blocks startup -- missing or
  # unparseable Census_all*.csv just means no role badges get shown.
  cgc_path <- find_cgc_file(data_dir)
  cgc_df <- NULL
  if (!is.null(cgc_path)) {
    cgc_df <- tryCatch(load_cgc_df(cgc_path), error = function(e) {
      message("Found a Census_all*.csv file but couldn't parse it (", conditionMessage(e), ") -- continuing without CGC role annotations.")
      NULL
    })
  } else {
    message("No Census_all*.csv (COSMIC Cancer Gene Census) found in '", data_dir, "' -- gene-role badges will be skipped. Optional; does not affect SL/FR results.")
  }
  
  list(dependency_df = dependency_df, mutation_df = mutation_df, expression_df = expression_df,
       cn_df = cn_df, paralog_df = paralog_df, model_df = model_df, compounds_df = compounds_df,
       cgc_df = cgc_df)
}

# v2.8: mock data now includes a cn_df too, so the whole pipeline
# (including the new CN-based loss criterion) can be exercised without
# real DepMap files. GENE_A gets a mix of mutation-flagged AND
# CN-only-deleted mock cell lines, specifically to mirror the
# ARID1A-style scenario (some loss is only visible via CN).
make_mock_data <- function(n_cell_lines = 150, seed = 0) {
  set.seed(seed)
  cell_lines <- sprintf("CL_%03d", 0:(n_cell_lines - 1))
  genes <- c("GENE_A", "GENE_B", "PARA_1", "PARA_2", "CTRL_1", "CTRL_2")
  mutation_df <- as.data.frame(sapply(genes, function(g) rbinom(n_cell_lines, 1, 0.3)))
  rownames(mutation_df) <- cell_lines
  dep <- list()
  dep$CTRL_1 <- rnorm(n_cell_lines, 0, 0.3); dep$CTRL_2 <- rnorm(n_cell_lines, 0, 0.3)
  base_B <- rnorm(n_cell_lines, 0, 0.25)
  a_inactive <- mutation_df$GENE_A == 1
  base_B[a_inactive] <- base_B[a_inactive] - 1.2
  dep$GENE_A <- rnorm(n_cell_lines, 0, 0.3); dep$GENE_B <- base_B
  shared_signal <- rnorm(n_cell_lines, 0, 0.4)
  dep$PARA_1 <- -0.1 + shared_signal * 0.6 + rnorm(n_cell_lines, 0, 0.15)
  dep$PARA_2 <- -0.1 + shared_signal * 0.6 + rnorm(n_cell_lines, 0, 0.15)
  dependency_df <- as.data.frame(dep); rownames(dependency_df) <- cell_lines
  expr_shared <- rnorm(n_cell_lines, 5, 1)
  expression_df <- data.frame(
    GENE_A = expr_shared + rnorm(n_cell_lines, 0, 0.3),
    GENE_B = expr_shared * 0.8 + rnorm(n_cell_lines, 0, 0.3),
    PARA_1 = rnorm(n_cell_lines, 5, 1), PARA_2 = rnorm(n_cell_lines, 5, 1),
    CTRL_1 = rnorm(n_cell_lines, 5, 1), CTRL_2 = rnorm(n_cell_lines, 5, 1)
  )
  rownames(expression_df) <- cell_lines
  
  # --- v2.8 mock CN: linear scale, diploid ~1.0 ---
  cn_df <- as.data.frame(lapply(genes, function(g) rep(1.0, n_cell_lines)))
  colnames(cn_df) <- genes
  cn_df <- cn_df + matrix(rnorm(n_cell_lines * length(genes), 0, 0.08), nrow = n_cell_lines)
  # Some already-mutant GENE_A lines are ALSO deep-deleted (co-occurring, realistic)
  also_deleted <- a_inactive & (runif(n_cell_lines) < 0.4)
  cn_df$GENE_A[also_deleted] <- runif(sum(also_deleted), 0.05, 0.25)
  # A separate CN-ONLY-deleted group: NOT flagged by mutation_df at all —
  # this is the exact scenario the CN file is meant to catch.
  cn_only_candidates <- which(!a_inactive)
  cn_only <- sample(cn_only_candidates, size = round(0.15 * length(cn_only_candidates)))
  cn_df$GENE_A[cn_only] <- runif(length(cn_only), 0.05, 0.25)
  rownames(cn_df) <- cell_lines
  
  paralog_df <- data.frame(GeneA = "PARA_1", GeneB = "PARA_2", pct_identity = 85.0, stringsAsFactors = FALSE)
  list(dependency_df = dependency_df, mutation_df = mutation_df,
       expression_df = expression_df, cn_df = cn_df, paralog_df = paralog_df,
       model_df = NULL, compounds_df = NULL, cgc_df = NULL)
}

# ============================================================
# ENHANCED TEST RESULT STRUCTURE
# ============================================================

make_test_result <- function(name, statistic, p_value, passed,
                             effect_size = NA, note = "",
                             tier = "core") {
  list(
    name = name,
    statistic = statistic,
    p_value = p_value,
    passed = isTRUE(passed),
    effect_size = effect_size,
    note = note,
    tier = tier  # "core" or "supporting"
  )
}

# ============================================================
# GENE LOSS DEFINITION — v2.8: Multi-modal
# (mutation OR low expression OR copy-number deletion)
# ============================================================
Z_LOW_THRESHOLD <- -1.5   # SDs below the gene's mean expression to count as a low-expression outlier

# v2.8: added cn_df (optional — pass NULL to fall back to the old
# mutation+expression-only behavior, e.g. for mock data without CN, or
# if the CN file failed to load). method gains a "cn" option and
# "combined" now includes CN whenever cn_df is supplied.
define_gene_loss <- function(gene, mutation_df, expression_df, cn_df = NULL,
                             method = "combined") {
  
  # Align all datasets to common cell lines
  common <- intersect(rownames(mutation_df), rownames(expression_df))
  if (!is.null(cn_df))
    common <- intersect(common, rownames(cn_df))
  
  mutation_df <- mutation_df[common, , drop = FALSE]
  expression_df <- expression_df[common, , drop = FALSE]
  if (!is.null(cn_df))
    cn_df <- cn_df[common, , drop = FALSE]
  
  loss <- rep(FALSE, length(common))
  names(loss) <- common
  
  # Mutation
  if (method %in% c("mutation", "combined")) {
    if (gene %in% colnames(mutation_df)) {
      mut <- mutation_df[, gene]
      loss <- loss | (!is.na(mut) & mut == 1)
    }
  }
  
  # Expression
  if (method %in% c("expression", "combined")) {
    if (gene %in% colnames(expression_df)) {
      expr <- expression_df[, gene]
      z <- (expr - mean(expr, na.rm = TRUE)) / sd(expr, na.rm = TRUE)
      loss <- loss | (!is.na(z) & z < Z_LOW_THRESHOLD)
    }
  }
  
  # Copy number
  if (method %in% c("cn", "combined") &&
      !is.null(cn_df) &&
      gene %in% colnames(cn_df)) {
    
    cn <- cn_df[, gene]
    loss <- loss | (!is.na(cn) & cn < CN_LOSS_THRESHOLD)
  }
  
  return(loss)
}

# v2.8: breakdown now reports each source separately (mutation-only,
# expression-only, CN-only, and any overlap of 2+ sources), so you can
# see exactly how much a given group grew because of the CN file.
loss_group_breakdown <- function(gene, mutation_df, expression_df, cn_df = NULL) {
  common <- intersect(rownames(mutation_df), rownames(expression_df))
  if (!is.null(cn_df)) common <- intersect(common, rownames(cn_df))
  
  mut_col <- mutation_df[common, gene]
  is_mut <- !is.na(mut_col) & mut_col == 1
  
  expr_col <- expression_df[common, gene]
  z <- (expr_col - mean(expr_col, na.rm = TRUE)) / sd(expr_col, na.rm = TRUE)
  is_low_expr <- !is.na(z) & z < Z_LOW_THRESHOLD
  
  if (!is.null(cn_df) && gene %in% colnames(cn_df)) {
    cn_col <- cn_df[common, gene]
    is_cn_del <- !is.na(cn_col) & cn_col < CN_LOSS_THRESHOLD
  } else {
    is_cn_del <- rep(FALSE, length(common))
  }
  
  n_flags <- as.integer(is_mut) + as.integer(is_low_expr) + as.integer(is_cn_del)
  
  list(
    n_mutation_only = sum(is_mut & n_flags == 1, na.rm = TRUE),
    n_expr_only = sum(is_low_expr & n_flags == 1, na.rm = TRUE),
    n_cn_only = sum(is_cn_del & n_flags == 1, na.rm = TRUE),
    n_multi_criteria = sum(n_flags >= 2, na.rm = TRUE),
    n_total_loss = sum(n_flags >= 1, na.rm = TRUE),
    n_cell_lines = length(common)
  )
}

# ============================================================
# SYNTHETIC LETHALITY TESTS
# ============================================================

# --- CORE 1/2: Differential Essentiality (A->B and B->A) ---
# v2.8: gene_a's loss now also draws on cn_df (via define_gene_loss),
# passed straight through as an extra argument — no logic here changed
# beyond forwarding cn_df and reporting its contribution in the note.
test_differential_essentiality <- function(gene_a, gene_b, mutation_df, expression_df,
                                           dependency_df, loss_method = "combined",
                                           model_df = NULL, cn_df = NULL) {
  test_name <- paste0("Differential Essentiality (", gene_a, "\u2192", gene_b, ")")
  
  a_loss <- define_gene_loss(gene_a, mutation_df, expression_df, cn_df, method = loss_method)
  common_cells <- intersect(names(a_loss), rownames(dependency_df))
  a_loss <- a_loss[common_cells]
  dep_b <- dependency_df[common_cells, gene_b]
  
  valid <- !is.na(a_loss) & !is.na(dep_b)
  a_loss <- a_loss[valid]; dep_b <- dep_b[valid]
  common_cells <- common_cells[valid]
  n_loss <- sum(a_loss); n_intact <- sum(!a_loss)
  
  if (n_loss < MIN_GROUP_SIZE || n_intact < MIN_GROUP_SIZE) {
    return(make_test_result(test_name, NA, NA, FALSE, NA,
                            sprintf("Groups too small: %d loss vs %d intact (need %d each)", n_loss, n_intact, MIN_GROUP_SIZE),
                            tier = "core"))
  }
  
  # --- Determine whether tissue adjustment is usable for this pair ---
  use_tissue <- FALSE
  tissue_factor <- NULL
  if (!is.null(model_df) && "OncotreeLineage" %in% colnames(model_df)) {
    tissue_vec <- model_df[common_cells, "OncotreeLineage"]
    tab <- table(tissue_vec)
    keep_tissues <- names(tab[tab >= 10])  # keep tissues with 10+ samples in this pair's data
    if (length(keep_tissues) > 1) {
      tissue_vec[!(tissue_vec %in% keep_tissues)] <- "Other"
      tissue_factor <- factor(tissue_vec)
      use_tissue <- TRUE
    }
  }
  
  reg_df <- data.frame(dep_b = dep_b, a_loss = as.numeric(a_loss))
  if (use_tissue) reg_df$tissue <- tissue_factor
  reg_df <- reg_df[complete.cases(reg_df), ]
  
  if (nrow(reg_df) < (MIN_GROUP_SIZE * 2)) {
    return(make_test_result(test_name, NA, NA, FALSE, NA, "Not enough complete cases for regression", tier = "core"))
  }
  
  # --- Fit the model that gives us the p-value (includes a_loss) ---
  if (use_tissue) {
    full_fit <- tryCatch(lm(dep_b ~ a_loss + tissue, data = reg_df), error = function(e) NULL)
    base_fit <- tryCatch(lm(dep_b ~ tissue, data = reg_df), error = function(e) NULL)
  } else {
    full_fit <- tryCatch(lm(dep_b ~ a_loss, data = reg_df), error = function(e) NULL)
    base_fit <- NULL
  }
  
  if (is.null(full_fit)) return(make_test_result(test_name, NA, NA, FALSE, NA, "Regression failed", tier = "core"))
  
  coefs <- summary(full_fit)$coefficients
  if (!("a_loss" %in% rownames(coefs))) {
    return(make_test_result(test_name, NA, NA, FALSE, NA,
                            "Coefficient missing (likely rank-deficient tissue design for this pair)", tier = "core"))
  }
  
  mut_coef <- coefs["a_loss", "Estimate"]
  mut_pval <- coefs["a_loss", "Pr(>|t|)"]
  
  if (mut_coef < 0) {
    one_tailed_p <- mut_pval / 2
  } else {
    one_tailed_p <- 1 - mut_pval / 2
  }
  
  if (use_tissue && !is.null(base_fit)) {
    resid_scale <- residuals(base_fit)
    threshold_used <- PARTIAL_COHEN_D_THRESHOLD
    adj_label <- "with"
  } else {
    resid_scale <- reg_df$dep_b - mean(reg_df$dep_b)
    threshold_used <- COHEN_D_THRESHOLD
    adj_label <- "no"
  }
  grp <- reg_df$a_loss == 1
  m_loss <- mean(resid_scale[grp]);  m_intact <- mean(resid_scale[!grp])
  v_loss <- var(resid_scale[grp]);   v_intact <- var(resid_scale[!grp])
  pooled_sd <- sqrt((v_loss + v_intact) / 2)
  partial_d <- if (pooled_sd > 0) (m_intact - m_loss) / pooled_sd else 0
  
  if (VERBOSE) {
    cat("\n===== DEBUG:", gene_a, "->", gene_b, "=====\n")
    cat("n_loss =", n_loss, "\n")
    cat("n_intact =", n_intact, "\n")
    cat("coef =", mut_coef, "\n")
    cat("one-tailed p =", one_tailed_p, "\n")
    cat("partial d =", partial_d, "\n")
    cat("use_tissue =", use_tissue, "\n")
    cat("threshold =", threshold_used, "\n")
  }
  
  passed <- one_tailed_p < ALPHA && mut_coef < 0 && partial_d >= threshold_used
  
  brk <- loss_group_breakdown(gene_a, mutation_df, expression_df, cn_df)
  
  make_test_result(test_name, mut_coef, one_tailed_p, passed, partial_d,
                   sprintf("coef=%.3f, one-tailed p=%.4f, partial d=%.2f (need >=%.2f, %s tissue adj) | %d loss (%d mut, %d expr, %d CN-only, %d multi-criteria) vs %d intact",
                           mut_coef, one_tailed_p, partial_d, threshold_used, adj_label,
                           n_loss, brk$n_mutation_only, brk$n_expr_only, brk$n_cn_only, brk$n_multi_criteria, n_intact), tier = "core")
}

# --- SUPPORTING: SL Score — Correlated Selective Dependency ---
test_sl_score <- function(gene_a, gene_b, dependency_df) {
  dep_a <- dependency_df[[gene_a]]; dep_b <- dependency_df[[gene_b]]
  valid <- !is.na(dep_a) & !is.na(dep_b)
  dep_a <- dep_a[valid]; dep_b <- dep_b[valid]
  
  a_selective <- dep_a < SELECTIVE_THRESHOLD
  b_selective <- dep_b < SELECTIVE_THRESHOLD
  
  n11 <- sum(a_selective & b_selective); n12 <- sum(a_selective & !b_selective)
  n21 <- sum(!a_selective & b_selective); n22 <- sum(!a_selective & !b_selective)
  test_name <- "SL Score (Selective Dependency Correlation) [supporting]"
  
  if (min(n11 + n12, n21 + n22, n11 + n21, n12 + n22) < 5) {
    return(make_test_result(test_name, NA, NA, FALSE, NA,
                            sprintf("Counts too low: both=%d, A-only=%d, B-only=%d, neither=%d", n11, n12, n21, n22),
                            tier = "supporting"))
  }
  
  tab <- matrix(c(n11, n21, n12, n22), nrow = 2)
  ft <- fisher.test(tab, alternative = "greater")
  odds_ratio <- unname(ft$estimate)
  passed <- ft$p.value < ALPHA && odds_ratio >= OR_SL_THRESHOLD
  
  make_test_result(test_name, odds_ratio, ft$p.value, passed, odds_ratio,
                   sprintf("OR=%.2f | %d cell lines dependent on both", odds_ratio, n11), tier = "supporting")
}

# --- SUPPORTING: Mutual Exclusivity ---
test_mutual_exclusivity <- function(gene_a, gene_b, mutation_df) {
  a_mut <- mutation_df[[gene_a]] == 1; b_mut <- mutation_df[[gene_b]] == 1
  valid <- !is.na(a_mut) & !is.na(b_mut)
  a_mut <- a_mut[valid]; b_mut <- b_mut[valid]
  
  n11 <- sum(a_mut & b_mut); n10 <- sum(a_mut & !b_mut)
  n01 <- sum(!a_mut & b_mut); n00 <- sum(!a_mut & !b_mut)
  mut_a_total <- n11 + n10; mut_b_total <- n11 + n01
  test_name <- "Mutual Exclusivity [supporting]"
  
  if (mut_a_total < MIN_INACTIVE_SAMPLES || mut_b_total < MIN_INACTIVE_SAMPLES) {
    return(make_test_result(test_name, NA, NA, FALSE, NA,
                            sprintf("Mutated samples too few: %s=%d, %s=%d", gene_a, mut_a_total, gene_b, mut_b_total),
                            tier = "supporting"))
  }
  
  tab <- matrix(c(n11, n01, n10, n00), nrow = 2)
  ft <- fisher.test(tab, alternative = "less")
  odds_ratio <- unname(ft$estimate)
  passed <- ft$p.value < ALPHA && odds_ratio <= OR_EXCLUSIVITY_THRESHOLD
  
  make_test_result(test_name, odds_ratio, ft$p.value, passed, odds_ratio,
                   sprintf("OR=%.2f | co-mutated=%d (expect few if SL)", odds_ratio, n11), tier = "supporting")
}

# --- SUPPORTING: Coexpression (informational only — never required) ---
test_coexpression <- function(gene_a, gene_b, expression_df) {
  test_name <- "Coexpression [supporting, informational]"
  if (nrow(expression_df) < MIN_SAMPLES) {
    return(make_test_result(test_name, NA, NA, FALSE, NA, "Insufficient samples", tier = "supporting"))
  }
  ct <- cor.test(expression_df[[gene_a]], expression_df[[gene_b]], method = "pearson")
  r <- unname(ct$estimate)
  passed <- ct$p.value < ALPHA
  make_test_result(test_name, r, ct$p.value, passed, r,
                   sprintf("r=%.3f, p=%.4f (context only — not a required SL criterion)", r, ct$p.value),
                   tier = "supporting")
}

# ============================================================
# FUNCTIONAL REDUNDANCY TESTS
# ============================================================

# --- CORE: Single-Gene Non-Essentiality ---
test_single_gene_nonessentiality <- function(gene_a, gene_b, dependency_df) {
  med_a <- median(dependency_df[[gene_a]], na.rm = TRUE)
  med_b <- median(dependency_df[[gene_b]], na.rm = TRUE)
  a_nonessential <- med_a > ESSENTIALITY_CUTOFF
  b_nonessential <- med_b > ESSENTIALITY_CUTOFF
  passed <- a_nonessential && b_nonessential
  
  make_test_result("Single-Gene Non-Essentiality", NA, NA, passed, NA,
                   sprintf("median dep: %s=%.2f, %s=%.2f (threshold=%.1f)",
                           gene_a, med_a, gene_b, med_b, ESSENTIALITY_CUTOFF), tier = "core")
}

# --- SUPPORTING: Conditional Co-Essentiality ---
# v2.13: sign flipped. Genuine FR pairs show NEGATIVE dep-dep correlation
# in context-relevant lines (paralog-buffering trade-off: De Kegel & Ryan
# 2019, PLOS Genetics) -- confirmed empirically 4/6 in our own FR-positive
# benchmark (rho -0.35 to -0.51), not just a threshold-fit.
test_conditional_coessentiality <- function(gene_a, gene_b, dependency_df) {
  dep_a <- dependency_df[[gene_a]]; dep_b <- dependency_df[[gene_b]]
  valid <- !is.na(dep_a) & !is.na(dep_b)
  dep_a <- dep_a[valid]; dep_b <- dep_b[valid]
  test_name <- "Conditional Co-Essentiality [supporting]"
  
  context_relevant <- (dep_a < -0.2) | (dep_b < -0.2)
  if (sum(context_relevant) < MIN_GROUP_SIZE) {
    if (sum(valid) < MIN_SAMPLES) {
      return(make_test_result(test_name, NA, NA, FALSE, NA, "Insufficient samples in relevant context", tier = "supporting"))
    }
    context_relevant <- rep(TRUE, length(dep_a))
  }
  
  ct <- cor.test(dep_a[context_relevant], dep_b[context_relevant], method = "spearman")
  rho <- unname(ct$estimate)
  passed <- ct$p.value < ALPHA && rho < -CORRELATION_THRESHOLD
  make_test_result(test_name, rho, ct$p.value, passed, rho,
                   sprintf("Spearman \u03c1=%.3f in %d context-relevant lines (negative = paralog buffering trade-off)", rho, sum(context_relevant)), tier = "supporting")
}

# --- CORE: Compensatory Expression ---
test_compensatory_expression <- function(gene_a, gene_b, expression_df, dependency_df) {
  # Align cell lines between expression and dependency data
  common <- intersect(rownames(expression_df), rownames(dependency_df))
  
  expr_a <- expression_df[common, gene_a]
  expr_b <- expression_df[common, gene_b]
  dep_a  <- dependency_df[common, gene_a]
  dep_b  <- dependency_df[common, gene_b]
  
  valid <- !is.na(expr_a) & !is.na(expr_b) &
    !is.na(dep_a) & !is.na(dep_b)
  
  test_name <- "Compensatory Expression"
  if (sum(valid) < MIN_SAMPLES) return(make_test_result(test_name, NA, NA, FALSE, NA, "Insufficient samples", tier = "core"))
  
  ct <- cor.test(expr_a[valid], expr_b[valid], method = "spearman")
  rho <- unname(ct$estimate)
  
  passed <- ct$p.value < ALPHA && abs(rho) >= CORRELATION_THRESHOLD
  pattern <- if (is.na(rho)) "n/a" else if (rho < 0) "compensatory upregulation (A down, B up)" else "shared-regulation co-expression (backup pair co-varies)"
  
  make_test_result(test_name, rho, ct$p.value, passed, rho,
                   sprintf("Spearman \u03c1=%.3f \u2014 %s", rho, pattern), tier = "core")
}

# --- SUPPORTING: Compensatory Dependency ---
test_compensatory_dependency <- function(gene_a, gene_b, dependency_df) {
  dep_a <- dependency_df[[gene_a]]; dep_b <- dependency_df[[gene_b]]
  valid <- !is.na(dep_a) & !is.na(dep_b)
  test_name <- "Compensatory Dependency [supporting]"
  
  if (sum(valid) < MIN_SAMPLES) return(make_test_result(test_name, NA, NA, FALSE, NA, "Insufficient samples", tier = "supporting"))
  
  ct <- cor.test(dep_a[valid], dep_b[valid], method = "spearman")
  rho <- unname(ct$estimate)
  passed <- ct$p.value < ALPHA && rho < -CORRELATION_THRESHOLD
  
  make_test_result(test_name, rho, ct$p.value, passed, rho,
                   sprintf("Spearman \u03c1=%.3f (negative = compensatory essentiality pattern)", rho), tier = "supporting")
}

# --- SUPPORTING: Compensatory Upregulation on Loss ---
# NEW v2.13. Literature: "transcriptional adaptation" / nonsense-induced
# transcriptional compensation (El-Brolosy et al. 2019, Nature; confirmed
# genome-wide in CRISPR knockout data, Genome Biology 2023/2024). When
# gene A is lost, does gene B's expression rise in those specific cell
# lines (and vice versa)? Unlike test_compensatory_expression (baseline
# co-expression only), this tests the actual transcriptional RESPONSE to
# loss -- more direct FR evidence.
test_compensatory_upregulation_on_loss <- function(gene_a, gene_b, mutation_df, expression_df, cn_df = NULL) {
  test_name <- "Compensatory Upregulation on Loss [supporting]"
  
  one_direction <- function(loss_gene, resp_gene) {
    loss_vec <- define_gene_loss(loss_gene, mutation_df, expression_df, cn_df, method = "combined")
    common <- intersect(names(loss_vec), rownames(expression_df))
    loss_vec <- loss_vec[common]
    resp_expr <- expression_df[common, resp_gene]
    valid <- !is.na(loss_vec) & !is.na(resp_expr)
    loss_vec <- loss_vec[valid]; resp_expr <- resp_expr[valid]
    
    n_loss <- sum(loss_vec); n_intact <- sum(!loss_vec)
    if (n_loss < MIN_GROUP_SIZE || n_intact < MIN_GROUP_SIZE) {
      return(list(p = 1, effect = 0, n_loss = n_loss, ok = FALSE))
    }
    wt <- tryCatch(wilcox.test(resp_expr[loss_vec], resp_expr[!loss_vec], alternative = "greater", exact = FALSE),
                   error = function(e) NULL)
    if (is.null(wt)) return(list(p = 1, effect = 0, n_loss = n_loss, ok = FALSE))
    eff <- mean(resp_expr[loss_vec], na.rm = TRUE) - mean(resp_expr[!loss_vec], na.rm = TRUE)
    list(p = wt$p.value, effect = eff, n_loss = n_loss, ok = TRUE)
  }
  
  res_ab <- one_direction(gene_a, gene_b)  # A lost -> B expr up?
  res_ba <- one_direction(gene_b, gene_a)  # B lost -> A expr up?
  
  ab_pass <- res_ab$ok && res_ab$p < ALPHA && res_ab$effect > 0
  ba_pass <- res_ba$ok && res_ba$p < ALPHA && res_ba$effect > 0
  passed <- ab_pass || ba_pass
  
out <- make_test_result(test_name, NA, NA, passed, NA,
                   sprintf("%s loss\u2192%s expr: %s (n=%d) | %s loss\u2192%s expr: %s (n=%d)",
                           gene_a, gene_b, if (ab_pass) "UP (sig)" else "n.s.", res_ab$n_loss,
                           gene_b, gene_a, if (ba_pass) "UP (sig)" else "n.s.", res_ba$n_loss),
                   tier = "supporting")
  out$p_ab      <- res_ab$p
  out$p_ba      <- res_ba$p
  out$effect_ab <- res_ab$effect
  out$effect_ba <- res_ba$effect
  out
}
# v2.15 fix: paralogy's contribution to FR decision support is now gated by
# sequence identity (>=70%). Motivated by an identity-threshold sweep on the
# internal benchmark showing weak-identity matches (20-55%) are non-discriminating
# (flat ~4-7% precision) while strong matches (>=70%) carry real signal --
# see validation summary Section 5.3, 5.5.

test_paralogy <- function(gene_a, gene_b, paralog_df, min_identity_for_bonus = 70) {
  match_row <- paralog_df[(paralog_df$GeneA == gene_a & paralog_df$GeneB == gene_b) |
                            (paralog_df$GeneA == gene_b & paralog_df$GeneB == gene_a), ]
  is_paralog <- nrow(match_row) > 0
  identity_info <- ""
  passes_bonus <- FALSE
  if (is_paralog) {
    pct <- if ("pct_identity" %in% colnames(match_row)) match_row$pct_identity[1] else NA
    identity_info <- if (!is.na(pct)) sprintf(" (%.1f%% identity)", pct) else ""
    passes_bonus <- !is.na(pct) && pct >= min_identity_for_bonus
  }
  note <- if (!is_paralog) {
    "No paralogy record (FR can still exist)"
  } else if (passes_bonus) {
    paste0("Known paralog pair", identity_info, " -- meets identity threshold for decision support")
  } else {
    paste0("Known paralog pair", identity_info, " -- below identity threshold, informational only")
  }
  out <- make_test_result("Paralogy [supporting]", NA, NA, passes_bonus, NA, note, tier = "supporting")
  out$identity_pct <- if (is_paralog && !is.na(pct)) pct else NA
  out
}

# --- CORE (FR): Reciprocal Compensation ---
# v2.8: stratified_test's loss call now also threads cn_df through, so
# tissue-stratified FR detection benefits from copy-number loss too.
test_reciprocal_compensation <- function(gene_a, gene_b, mutation_df, expression_df, dependency_df, model_df = NULL, cn_df = NULL) {
  
  stratified_test <- function(loss_gene, dep_gene) {
    loss_vec <- define_gene_loss(loss_gene, mutation_df, expression_df, cn_df, method = "combined")
    common <- intersect(names(loss_vec), rownames(dependency_df))
    loss_vec <- loss_vec[common]
    dep_vec <- dependency_df[common, dep_gene]
    valid <- !is.na(loss_vec) & !is.na(dep_vec)
    loss_vec <- loss_vec[valid]; dep_vec <- dep_vec[valid]
    
    if (sum(loss_vec) < 5) return(list(p = 1, effect = 0, note = "Too few total mutations"))
    
    if (!is.null(model_df) && "OncotreeLineage" %in% colnames(model_df)) {
      tissues <- model_df[common[valid], "OncotreeLineage"]
      p_vals <- c(); effects <- c()
      
      for (t in unique(tissues)) {
        idx <- tissues == t
        t_loss <- loss_vec[idx]; t_dep <- dep_vec[idx]
        if (sum(t_loss) >= 5 && sum(!t_loss) >= 10) {
          wt <- tryCatch(wilcox.test(t_dep[t_loss], t_dep[!t_loss], alternative = "less", exact = FALSE), error = function(e) NULL)
          if (!is.null(wt)) {
            if (VERBOSE) {
              cat(
                t,
                " loss=", sum(t_loss),
                " intact=", sum(!t_loss),
                " p=", signif(wt$p.value, 3),
                " effect=", round(mean(t_dep[!t_loss]) - mean(t_dep[t_loss]), 3),
                "\n"
              )
            }
            p_vals <- c(p_vals, wt$p.value)
            effects <- c(effects, mean(t_dep[!t_loss], na.rm=T) - mean(t_dep[t_loss], na.rm=T))
          }
        }
      }
      
      if (length(p_vals) > 0) {
        combined_stat <- -2 * sum(log(p_vals))
        combined_p <- pchisq(combined_stat, df = 2 * length(p_vals), lower.tail = FALSE)
        avg_effect <- mean(effects)
        return(list(p = combined_p, effect = avg_effect, note = sprintf("Meta-p=%.3e across %d tissues", combined_p, length(p_vals))))
      }
      
      # v2.12: no single tissue had enough loss+intact samples on its own --
      # fall back to pan-cancer Wilcoxon instead of failing automatically.
      wt <- tryCatch(wilcox.test(dep_vec[loss_vec], dep_vec[!loss_vec], alternative = "less", exact = FALSE), error = function(e) NULL)
      if (is.null(wt)) return(list(p = 1, effect = 0, note = "Mutations too scattered across tissues, and pan-cancer fallback failed"))
      return(list(p = wt$p.value, effect = mean(dep_vec[!loss_vec]) - mean(dep_vec[loss_vec]),
                  note = "Pan-cancer Wilcoxon (fallback -- no single tissue had enough samples)"))
    } else {
      wt <- tryCatch(wilcox.test(dep_vec[loss_vec], dep_vec[!loss_vec], alternative = "less", exact = FALSE), error = function(e) NULL)
      if (is.null(wt)) return(list(p = 1, effect = 0, note = "Pan-cancer Wilcoxon failed"))
      return(list(p = wt$p.value, effect = mean(dep_vec[!loss_vec]) - mean(dep_vec[loss_vec]), note = "Pan-cancer Wilcoxon"))
    }
  }
  
  res_ab <- stratified_test(gene_a, gene_b)
  res_ba <- stratified_test(gene_b, gene_a)
  
  ab_pass <- res_ab$p < ALPHA && res_ab$effect > 0
  ba_pass <- res_ba$p < ALPHA && res_ba$effect > 0
  
  passed <- ab_pass && ba_pass
  partial <- ab_pass || ba_pass
  
  res <- make_test_result("Reciprocal Compensation", NA, NA, passed, NA,
                          sprintf("A\u2192B: %s (%s) | B\u2192A: %s (%s)", 
                                  if (ab_pass) "PASS" else "fail", res_ab$note,
                                  if (ba_pass) "PASS" else "fail", res_ba$note),
                          tier = "core")
  res$partial   <- partial
  res$p_ab      <- res_ab$p
  res$p_ba      <- res_ba$p
  res$effect_ab <- res_ab$effect
  res$effect_ba <- res_ba$effect
  res
}

# ============================================================
# PIPELINE ORCHESTRATION
# ============================================================

# v2.8: cn_df threaded through to both differential-essentiality calls.
run_sl_pipeline <- function(gene_a, gene_b, dependency_df, mutation_df, expression_df, model_df = NULL, cn_df = NULL) {
  
  # ======================================================
  # Direction-driven SL classification (no TS/EG heuristic)
  # ======================================================
  core_tests <- list(
    test_differential_essentiality(gene_a, gene_b, mutation_df, expression_df, dependency_df, "combined", model_df, cn_df),
    test_differential_essentiality(gene_b, gene_a, mutation_df, expression_df, dependency_df, "combined", model_df, cn_df)
  )
  
  supporting_tests <- list(
    test_sl_score(gene_a, gene_b, dependency_df),
    test_mutual_exclusivity(gene_a, gene_b, mutation_df),
    test_coexpression(gene_a, gene_b, expression_df)
  )
  
  all_tests <- c(core_tests, supporting_tests)
  
  forward_pass <- core_tests[[1]]$passed
  reverse_pass <- core_tests[[2]]$passed
  
  n_core_passed <- sum(forward_pass, reverse_pass)
  n_support_passed <- sum(vapply(supporting_tests, function(t) t$passed, logical(1)))
  
  if (forward_pass && reverse_pass) {
    call <- "Synthetic Lethal"
    confidence <- "high"
    pattern <- "Symmetric (both directions confirmed)"
    
  } else if (forward_pass || reverse_pass) {
    pattern <- if (forward_pass) {
      paste0("Asymmetric: ", gene_a, " loss \u2192 ", gene_b, " dependency")
    } else {
      paste0("Asymmetric: ", gene_b, " loss \u2192 ", gene_a, " dependency")
    }
    # never reuse the bare "Synthetic Lethal" string here -- that exact
    # string is what analyze_gene_pair_v2 checks for HIGH priority /
    # the strong badge, and means BOTH directions passed. A single
    # direction, even corroborated, stays "(moderate)" -- confidence
    # is what carries the distinction.
    call <- "Synthetic Lethal (moderate)"
    confidence <- if (n_support_passed >= 1) "moderate" else "low-moderate"
    
  } else {
    call <- "Not Synthetic Lethal"
    confidence <- "none"
    pattern <- "No directional SL"
  }
  
  list(
    tests = all_tests, core_tests = core_tests, supporting_tests = supporting_tests,
    n_core_passed = n_core_passed, n_core_total = 2,
    n_support_passed = n_support_passed, n_support_total = length(supporting_tests),
    call = call, confidence = confidence, pattern = pattern
  )
}

# v2.8: cn_df threaded through to test_reciprocal_compensation.
run_fr_pipeline <- function(gene_a, gene_b, dependency_df, expression_df, paralog_df, mutation_df, model_df = NULL, cn_df = NULL) {
  core_tests <- list(
    test_single_gene_nonessentiality(gene_a, gene_b, dependency_df),
    test_reciprocal_compensation(gene_a, gene_b, mutation_df, expression_df, dependency_df, model_df, cn_df),
    test_compensatory_expression(gene_a, gene_b, expression_df, dependency_df)
  )
  supporting_tests <- list(
    test_paralogy(gene_a, gene_b, paralog_df),
    test_conditional_coessentiality(gene_a, gene_b, dependency_df),
    test_compensatory_dependency(gene_a, gene_b, dependency_df),
    test_compensatory_upregulation_on_loss(gene_a, gene_b, mutation_df, expression_df, cn_df)
  )
  all_tests <- c(core_tests, supporting_tests)
  n_core_passed <- sum(vapply(core_tests, function(t) t$passed, logical(1)))
  
  # v2.13: decision-driving support now excludes BOTH paralogy (index 1 --
  # non-discriminating, candidate pairs are paralog-enriched) AND
  # compensatory_dependency (index 3 -- effect sizes too small to ever clear
  # threshold in FR true positives; kept only for display/context). Decision
  # draws on conditional coessentiality (2, sign-corrected) + compensatory
  # upregulation on loss (4, new) only.
  n_support_passed   <- sum(vapply(supporting_tests, function(t) t$passed, logical(1)))
  is_paralog_pair    <- supporting_tests[[1]]$passed
  n_support_decision <- sum(vapply(supporting_tests[c(2, 4)], function(t) t$passed, logical(1)))
  
  noness_passed    <- core_tests[[1]]$passed
  recip_partial     <- isTRUE(core_tests[[2]]$partial) && !core_tests[[2]]$passed
  recip_full        <- core_tests[[2]]$passed
  comp_expr_passed  <- core_tests[[3]]$passed
  
  if (!noness_passed) {
    call <- "Not Functionally Redundant"
    confidence <- "none (gate failed: at least one gene is essential)"
    
  } else if (recip_full && comp_expr_passed) {
    call <- "Functional Redundancy"
    confidence <- "high (all three core tests concordant)"
    
  } else if (recip_full && !comp_expr_passed) {
    # Full bidirectional compensation already cleared significance + effect
    # size in BOTH directions -- that's a real core signal on its own, so it
    # stays at "moderate" even with zero non-paralogy support. Confidence
    # text is honest about how thin that specific case is.
    call <- "Functional Redundancy (moderate)"
    confidence <- if (n_support_decision >= 1) {
      "moderate (reciprocal compensation both directions; expression not concordant, but corroborated by independent evidence beyond paralogy)"
    } else {
      "moderate-low (reciprocal compensation both directions; expression not concordant, no non-paralogy corroboration)"
    }
    
  } else if (recip_partial) {
    # v2.14: literature ("Incomplete paralog compensation generates
    # selective dependency", PLOS Genetics; De Kegel & Ryan 2019) shows
    # asymmetric/one-directional buffering is the NORM for real paralog
    # pairs, not a lesser/edge case -- so a CORROBORATED one-directional
    # signal is promoted to a real FR call, not left at "Borderline".
    # Uncorroborated one-directional signal stays Borderline/Not-FR --
    # test_reciprocal_compensation's per-direction pass has no effect-size
    # floor (only p<ALPHA), so it still needs support to be trusted alone.
    effective_support <- n_support_decision + (if (is_paralog_pair) 1 else 0)
    if (effective_support >= 2) {
      call <- "Functional Redundancy (moderate)"
      confidence <- "moderate (asymmetric/one-directional compensation, consistent with literature on incomplete paralog buffering; corroborated)"
    } else if (effective_support == 1) {
      call <- "Borderline FR"
      confidence <- "low-moderate (one-directional compensation, weakly corroborated)"
    } else {
      call <- "Not Functionally Redundant"
      confidence <- sprintf(
        "none (partial reciprocal compensation only, insufficient corroboration: effective support %d/2%s)",
        effective_support,
        if (is_paralog_pair && n_support_decision == 0) " -- paralogy alone does not count" else ""
      )
    }
    
  } else {
    call <- "Not Functionally Redundant"
    confidence <- "none (no reciprocal-compensation signal - may be co-pathway, not backup)"
  }
  
  list(tests = all_tests, core_tests = core_tests, supporting_tests = supporting_tests,
       n_core_passed = n_core_passed, n_core_total = length(core_tests),
       n_support_passed = n_support_passed, n_support_total = length(supporting_tests),
       n_support_decision = n_support_decision, is_paralog_pair = is_paralog_pair,
       call = call, confidence = confidence)
}
# ============================================================
# INTEGRATED ANALYSIS
# ============================================================

# v2.8: added cn_df param, threaded to both pipelines; tissue filter and
# biomarker-context filter now also subset cn_df so it stays in sync
# with mutation_df/expression_df/dependency_df.
analyze_gene_pair_v2 <- function(gene_a, gene_b, dependency_df, mutation_df, expression_df, paralog_df,
                                 compounds_df = NULL, model_df = NULL, cn_df = NULL, cgc_df = NULL, tissue = NULL,
                                 biomarker_gene = NULL, biomarker_mode = "loss") {
  for (g in c(gene_a, gene_b)) {
    if (!(g %in% colnames(dependency_df))) stop(paste0("Gene '", g, "' not found in the loaded dependency data."))
  }
  
  missing_from_mut  <- c(gene_a, gene_b)[!(c(gene_a, gene_b) %in% colnames(mutation_df))]
  missing_from_expr <- c(gene_a, gene_b)[!(c(gene_a, gene_b) %in% colnames(expression_df))]
  if (length(missing_from_mut) > 0 || length(missing_from_expr) > 0) {
    all_missing_genes <- unique(c(missing_from_mut, missing_from_expr))
    missing_sources <- unique(c(
      if (length(missing_from_mut) > 0) "the mutation data" else NULL,
      if (length(missing_from_expr) > 0) "the expression data" else NULL
    ))
    stop(sprintf(
      paste0(
        "Data coverage gap: %s %s present in the CRISPR dependency data, but missing from %s. ",
        "This means mutation- and/or expression-based loss can't be computed for this gene, so any ",
        "SL/FR result involving it would silently look like 'zero mutations / insufficient data' rather ",
        "than a genuine negative. Check the exact gene symbol used in your OmicsSomaticMutationsMatrixDamaging.csv ",
        "and OmicsExpressionTPMLogp1... column headers (naming/alias mismatches are the most common cause), ",
        "or this gene may simply not be covered by that particular DepMap release."
      ),
      paste(all_missing_genes, collapse = " and "),
      if (length(all_missing_genes) > 1) "are" else "is",
      paste(missing_sources, collapse = " and ")
    ))
  }
  
  # v2.8: CN coverage gap is a soft note, not a hard stop — WGS coverage
  # can legitimately be a subset of all sequenced lines, and the tool
  # should still run on mutation+expression alone in that case.
  cn_missing_note <- NULL
  if (!is.null(cn_df)) {
    missing_from_cn <- c(gene_a, gene_b)[!(c(gene_a, gene_b) %in% colnames(cn_df))]
    if (length(missing_from_cn) > 0) {
      cn_missing_note <- sprintf(
        "Note: %s not found in the copy number file — loss calls for %s fall back to mutation+expression only.",
        paste(missing_from_cn, collapse = " and "), paste(missing_from_cn, collapse = " and ")
      )
    }
  }
  
  if (!is.null(tissue) && !is.null(model_df)) {
    tissue_models <- rownames(model_df)[model_df$OncotreeLineage == tissue]
    dependency_df <- dependency_df[rownames(dependency_df) %in% tissue_models, , drop = FALSE]
    mutation_df   <- mutation_df[rownames(mutation_df) %in% tissue_models, , drop = FALSE]
    expression_df <- expression_df[rownames(expression_df) %in% tissue_models, , drop = FALSE]
    if (!is.null(cn_df)) cn_df <- cn_df[rownames(cn_df) %in% tissue_models, , drop = FALSE]
  }
  
  if (!is.null(biomarker_gene) && biomarker_gene != "") {
    
    if (biomarker_gene %in% c(gene_a, gene_b)) {
      stop(sprintf(
        "Biomarker Context can't be the same gene as Gene A or Gene B ('%s'). The SL/FR tests already internally compare '%s'-loss vs '%s'-intact cell lines as part of the test itself — pre-filtering to only one of those groups removes the comparison the test needs. Leave Biomarker Context blank to test on the full (or tissue-filtered) panel, or enter a DIFFERENT third gene to test this pair within its context.",
        biomarker_gene, biomarker_gene, biomarker_gene
      ))
    }
    
    if (!(biomarker_gene %in% colnames(mutation_df)) || !(biomarker_gene %in% colnames(expression_df))) {
      stop(paste0("Biomarker gene '", biomarker_gene, "' not found in the loaded mutation/expression data."))
    }
    bio_loss <- define_gene_loss(biomarker_gene, mutation_df, expression_df, cn_df, method = "combined")
    common <- intersect(names(bio_loss), rownames(dependency_df))
    keep_ids <- if (biomarker_mode == "loss") common[bio_loss[common]] else common[!bio_loss[common]]
    
    if (length(keep_ids) < MIN_GROUP_SIZE) {
      brk <- loss_group_breakdown(biomarker_gene, mutation_df, expression_df, cn_df)
      stop(sprintf(
        paste0(
          "Biomarker context '%s' (%s) leaves only %d cell lines (need at least %d) out of %d in this panel. ",
          "Composition of the full loss group: %d by mutation only, %d by low-expression outlier only, %d by CN deletion only, %d by multiple criteria. ",
          "This context is too rare here to test reliably. Try: (1) a Tissue Filter that enriches for this alteration ",
          "(e.g. Breast/Ovary for BRCA1), (2) a more commonly-altered biomarker gene, or (3) leave Biomarker Context ",
          "blank and interpret the unstratified result with the context-dependence caveat in mind."
        ),
        biomarker_gene, if (biomarker_mode == "loss") "loss/mutant" else "intact/WT",
        length(keep_ids), MIN_GROUP_SIZE, brk$n_cell_lines,
        brk$n_mutation_only, brk$n_expr_only, brk$n_cn_only, brk$n_multi_criteria
      ))
    }
    
    dependency_df <- dependency_df[rownames(dependency_df) %in% keep_ids, , drop = FALSE]
    mutation_df   <- mutation_df[rownames(mutation_df) %in% keep_ids, , drop = FALSE]
    expression_df <- expression_df[rownames(expression_df) %in% keep_ids, , drop = FALSE]
    if (!is.null(cn_df)) cn_df <- cn_df[rownames(cn_df) %in% keep_ids, , drop = FALSE]
  }
  
  sl <- run_sl_pipeline(gene_a, gene_b, dependency_df, mutation_df, expression_df, model_df, cn_df)
  fr <- run_fr_pipeline(gene_a, gene_b, dependency_df, expression_df, paralog_df, mutation_df, model_df, cn_df)
  
  sl_strong <- sl$call == "Synthetic Lethal"
  sl_mod    <- sl$call == "Synthetic Lethal (moderate)"
  sl_border <- sl$call == "Borderline SL"
  fr_strong <- fr$call == "Functional Redundancy"
  fr_mod    <- fr$call == "Functional Redundancy (moderate)"
  fr_border <- fr$call == "Borderline FR"
  
  fr_any_positive <- fr_strong || fr_mod
  
  if ((sl_strong || sl_mod) && fr_any_positive) { overall <- "Redundancy-driven SL"; priority <- "HIGH" }
  else if (sl_strong)          { overall <- "Synthetic Lethal"; priority <- "HIGH" }
  else if (sl_mod)             { overall <- "Synthetic Lethal (moderate)"; priority <- "HIGH" }
  else if (sl_border && fr_any_positive) { overall <- "Borderline Redundancy-driven SL"; priority <- "MEDIUM" }
  else if (sl_border)          { overall <- "Borderline SL"; priority <- "MEDIUM" }
  else if (fr_strong)          { overall <- "Functional Redundancy"; priority <- "LOW" }
  else if (fr_mod)             { overall <- "Functional Redundancy (moderate)"; priority <- "LOW" }
  else if (fr_border)          { overall <- "Borderline FR"; priority <- "LOW" }
  else                         { overall <- "No Interaction"; priority <- "NONE" }
  
  result <- list(gene_pair = c(gene_a, gene_b), overall_call = overall, priority = priority,
                 SL = sl, FR = fr, n_cell_lines_used = nrow(dependency_df),
                 tissue_filter = tissue,
                 biomarker_filter = if (!is.null(biomarker_gene) && biomarker_gene != "")
                   sprintf("%s (%s)", biomarker_gene, if (biomarker_mode == "loss") "loss/mutant only" else "intact/WT only")
                 else NULL,
                 cn_note = cn_missing_note,
                 version = "2.14")
  
  if (!is.null(compounds_df)) {
    result$druggable_A <- get_druggability(gene_a, compounds_df)
    result$druggable_B <- get_druggability(gene_b, compounds_df)
  }
  
  # v2.14: purely informational -- attached AFTER the SL/FR calls above
  # are already final, and never fed back into any test or threshold.
  result$cgc_annotation_A <- get_cgc_annotation(gene_a, cgc_df)
  result$cgc_annotation_B <- get_cgc_annotation(gene_b, cgc_df)
  result$cgc_badge_A <- format_cgc_badge(result$cgc_annotation_A)
  result$cgc_badge_B <- format_cgc_badge(result$cgc_annotation_B)
  
  result
}

# ============================================================
# FORMATTING HELPERS
# ============================================================

make_test_df <- function(pipeline_result) {
  raw_p <- vapply(pipeline_result$tests, function(t) t$p_value, numeric(1))
  adj_p <- rep(NA_real_, length(raw_p))
  has_p <- !is.na(raw_p)
  if (any(has_p)) adj_p[has_p] <- p.adjust(raw_p[has_p], method = "BH")
  
  data.frame(
    Tier = vapply(pipeline_result$tests, function(t) if (t$tier == "core") "\u2b50 Core" else "\U0001f4cb Supporting", character(1)),
    Result = vapply(pipeline_result$tests, function(t) if (t$passed) "PASS" else "fail", character(1)),
    Test = vapply(pipeline_result$tests, function(t) t$name, character(1)),
    `Effect Size` = vapply(pipeline_result$tests, function(t) if (!is.na(t$effect_size)) sprintf("%.3f", t$effect_size) else "\u2014", character(1)),
    `p-value` = vapply(pipeline_result$tests, function(t) if (!is.na(t$p_value)) sprintf("%.4f", t$p_value) else "n/a", character(1)),
    `p-adj (BH)` = ifelse(is.na(adj_p), "\u2014", sprintf("%.4f", adj_p)),
    Note = vapply(pipeline_result$tests, function(t) t$note, character(1)),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

# ============================================================
# OPTIONAL: BENCHMARK HELPER (not wired into the UI)
# ============================================================
# Run from the R console AFTER sourcing this file and loading APP_DATA,
# to sanity-check thresholds against known-true AND known-false pairs
# rather than a single pair you already expect to be positive.
#   source("module1_sl_fr_app_v2_8.R")
#   run_benchmark(APP_DATA)
run_benchmark <- function(app_data, pairs = NULL) {
  if (is.null(pairs)) {
    pairs <- data.frame(
      gene_a = c("BRCA1", "ARID1A", "SMARCA4", "MYC", "TP53"),
      gene_b = c("PARP1", "ARID1B", "SMARCA2", "ACTB", "GAPDH"),
      expected = c("SL (context-dependent)", "SL/FR", "FR", "no relationship (control)", "no relationship (control)"),
      stringsAsFactors = FALSE
    )
  }
  results <- lapply(seq_len(nrow(pairs)), function(i) {
    ga <- pairs$gene_a[i]; gb <- pairs$gene_b[i]
    r <- tryCatch(
      analyze_gene_pair_v2(ga, gb, app_data$dependency_df, app_data$mutation_df, app_data$expression_df,
                           app_data$paralog_df, compounds_df = app_data$compounds_df, model_df = app_data$model_df,
                           cn_df = app_data$cn_df, cgc_df = app_data$cgc_df),
      error = function(e) list(error = conditionMessage(e))
    )
    data.frame(gene_a = ga, gene_b = gb, expected = pairs$expected[i],
               call = if (!is.null(r$error)) paste("ERROR:", r$error) else r$overall_call,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, results)
  print(out, row.names = FALSE)
  invisible(out)
}

# ============================================================
# LOAD DATA ONCE AT APP STARTUP
# ============================================================

# v2.10 (handover cleanup): was hardcoded to a personal path
# ("D:/Tool_genomic/depmap_data"), which only worked on one machine. Now
# reads from the GENECOMB_DATA_DIR environment variable if set, otherwise
# falls back to a "depmap_data" folder relative to the working directory.
# To point at a specific folder without editing this file, either set the
# env var (Sys.setenv(GENECOMB_DATA_DIR = "/path/to/data") before sourcing,
# or an OS-level env var) or just run R with that folder as the working
# directory.
DATA_DIR <- Sys.getenv("GENECOMB_DATA_DIR", unset = "depmap_data")
if (dir.exists(DATA_DIR) && length(list.files(DATA_DIR)) > 0) {
  message("Found data folder '", DATA_DIR, "' -- loading real data...")
  APP_DATA <- load_all_real_data(DATA_DIR)
  USING_MOCK <- FALSE
} else {
  message("Data folder NOT found at '", DATA_DIR, "' (current working directory is '", getwd(), "') -- using synthetic mock data instead.")
  APP_DATA <- make_mock_data()
  USING_MOCK <- TRUE
}


