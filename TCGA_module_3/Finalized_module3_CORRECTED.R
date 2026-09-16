# =====================================================================
# TCGA Gene-Pair Multi-Term Analysis -- STANDALONE TERMINAL APP (v5)
# -----------------------------------------------------------------------
# One self-contained file. No Shiny, no browser, no Bioconductor.
# Run it from a terminal:
#     Rscript gene_pair_terminal_app_v5.R
#
# DATA SOURCE: cBioPortal REST API (https://www.cbioportal.org/api)
#
# UPDATE LOG (v5 vs v4):
#   UPDATE 12 (subtype specificity burden adjustment):
#     - PROBLEM: When falling back to tumor stage (or any categorical
#       surrogate), Term 1 ran a raw Chi-square test. Later-stage tumors
#       naturally carry more mutations, so a "significant" stage
#       association could just reflect general mutation burden, not a
#       specific link to the gene pair.
#     - FIX: Term 1 now fits a logistic regression (Likelihood Ratio
#       Test) comparing a model with the subtype/stage variable vs.
#       a model with ONLY burden covariates. The headline p-value
#       follows the same suppression artifact protection (max of raw
#       and adjusted) as the other terms.
#   UPDATE 13 (double-hit arm-level co-deletion flag):
#     - PROBLEM: For genes on the same chromosomal arm (e.g., PBRM1 and
#       BAP1 on 3p), an arm-level loss affects BOTH genes. A mutated
#       sample will automatically test positive for "deletion" just from
#       the background arm event. This makes the "double-hit" rate
#       identical to the plain mutation rate, falsely implying true
#       biallelic inactivation when it's just a passenger arm loss.
#     - FIX: Term 4 now compares the double-hit rate to the plain
#       mutation rate. If they are nearly identical (within 2%), it
#       explicitly flags: "NOTE: Double-hit rate mirrors mutation
#       rate—likely arm-level co-deletion artifact."
#   UPDATE 14 (subtype attribute priority + cleanup, found during final
#   validation pass -- see tool_validation_summary.md / Complete Summary):
#     - PROBLEM: PAM50_SUBTYPE (the real molecular-subtype call) was listed
#       AFTER PAM50_CALL in Term 1's attribute search order, so on cohorts
#       carrying both (e.g. brca_tcga_pub) the tool silently used the less
#       specific PAM50_CALL instead. This was documented as a known,
#       unfixed limitation in an earlier validation draft; it is fixed here.
#     - FIX: PAM50_SUBTYPE now checked before PAM50_CALL.
#     - Also removed: a hardcoded local setwd() path, and leftover
#       development debug print statements (including one hardcoded to
#       print GLI1-specific values regardless of which genes were run).
#   (Carried over from v4 -- see prior versions for detail:)
#     - Suppression artifact protection (adjusted p cannot be more
#       significant than raw p).
#     - Mutation-level exclusivity for co-localized genes.
#     - Patient-level clinical data fetching.
# =====================================================================
# NOTE: no setwd() here on purpose -- a hardcoded local path (previously
# "D:/Tool_genomic/TCGA_module") breaks the script for anyone but the
# original author. This script does not read/write any local files itself
# (all data comes from the cBioPortal REST API), so no working directory
# is required. If you want output written somewhere specific, set your
# working directory in your own session before sourcing this file.
# ---- 0. Install & load packages (CRAN only -- no Bioconductor) ----
pkgs <- c("httr", "jsonlite", "dplyr", "tidyr", "tibble", "purrr", "msigdbr")
new_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(new_pkgs) > 0) install.packages(new_pkgs)

suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(msigdbr)
})
select   <- dplyr::select
filter   <- dplyr::filter
mutate   <- dplyr::mutate
arrange  <- dplyr::arrange
summarise <- dplyr::summarise
BASE_URL <- "https://www.cbioportal.org/api"

CANCER_GENE_LIST_PATH <- Sys.getenv("TCGA_CANCER_GENE_LIST", unset = "cancerGeneList.tsv")
if (file.exists(CANCER_GENE_LIST_PATH)) {
  .cgl <- read.delim(CANCER_GENE_LIST_PATH, stringsAsFactors = FALSE, check.names = FALSE)
  TUMOR_SUPPRESSORS <- toupper(.cgl$`Hugo Symbol`[.cgl$`Gene Type` %in% c("TSG", "ONCOGENE_AND_TSG")])
  ONCOGENES <- toupper(.cgl$`Hugo Symbol`[.cgl$`Gene Type` %in% c("ONCOGENE", "ONCOGENE_AND_TSG")])
  message(sprintf("Loaded cancer gene list: %d tumor suppressors, %d oncogenes.",
                  length(TUMOR_SUPPRESSORS), length(ONCOGENES)))
} else {
  warning("Cancer gene list file not found -- gene-role-aware CNA filtering will be disabled (all genes use broad mut+del+amp definition).")
  TUMOR_SUPPRESSORS <- character(0)
  ONCOGENES <- character(0)
}

check_gene_role <- function(gene) {
  g <- toupper(gene)
  if (g %in% TUMOR_SUPPRESSORS) return("tumor_suppressor")
  if (g %in% ONCOGENES) return("oncogene")
  message(sprintf("    !! '%s' not in cancer gene list -- using broad altered definition (mut+del+amp).", gene))
  return("unknown")
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ==================== 1. LOW-LEVEL API HELPERS ====================

HTTP_TIMEOUT_SEC <- 60
HTTP_MAX_RETRIES <- 3

with_retry <- function(expr_fn, max_tries = HTTP_MAX_RETRIES, label = "request") {
  last_err <- NULL
  for (attempt in seq_len(max_tries)) {
    outcome <- tryCatch(list(ok = TRUE, value = expr_fn()),
                        error = function(e) list(ok = FALSE, err = e))
    if (outcome$ok) return(outcome$value)
    last_err <- outcome$err
    message(sprintf("    !! %s attempt %d/%d failed: %s",
                    label, attempt, max_tries, conditionMessage(last_err)))
    if (attempt < max_tries) Sys.sleep(2 * attempt)
  }
  stop(last_err)
}

cbio_get <- function(path, query = list()) {
  with_retry(function() {
    resp <- httr::GET(paste0(BASE_URL, path), query = query, httr::accept_json(),
                      httr::timeout(HTTP_TIMEOUT_SEC))
    if (httr::http_error(resp)) {
      stop(sprintf("GET %s failed [%s]: %s", path, httr::status_code(resp),
                   httr::content(resp, as = "text", encoding = "UTF-8")))
    }
    jsonlite::fromJSON(httr::content(resp, as = "text", encoding = "UTF-8"), flatten = TRUE)
  }, label = paste("GET", path))
}

cbio_post <- function(path, body, query = list()) {
  with_retry(function() {
    resp <- httr::POST(paste0(BASE_URL, path), query = query,
                       body = jsonlite::toJSON(body, auto_unbox = TRUE), encode = "raw",
                       httr::content_type_json(), httr::accept_json(),
                       httr::timeout(HTTP_TIMEOUT_SEC))
    if (httr::http_error(resp)) {
      stop(sprintf("POST %s failed [%s]: %s", path, httr::status_code(resp),
                   httr::content(resp, as = "text", encoding = "UTF-8")))
    }
    txt <- httr::content(resp, as = "text", encoding = "UTF-8")
    if (identical(txt, "[]") || identical(txt, "")) return(tibble())
    jsonlite::fromJSON(txt, flatten = TRUE)
  }, label = paste("POST", path))
}

check_network_connectivity <- function() {
  ok <- tryCatch({
    resp <- httr::GET(paste0(BASE_URL, "/studies"), query = list(pageSize = 1),
                      httr::accept_json(), httr::timeout(10))
    !httr::http_error(resp)
  }, error = function(e) {
    message("\n!! Cannot reach cBioPortal (", BASE_URL, "):")
    message("   ", conditionMessage(e))
    message("   This looks like a DNS/network/proxy/firewall issue.")
    FALSE
  })
  ok
}

get_study_list <- function() {
  cbio_get("/studies") %>% as_tibble() %>% dplyr::select(studyId, name, cancerTypeId)
}

get_molecular_profiles <- function(study_id) {
  cbio_get(paste0("/studies/", study_id, "/molecular-profiles")) %>% as_tibble()
}

get_default_sample_list <- function(study_id) {
  lists <- cbio_get(paste0("/studies/", study_id, "/sample-lists")) %>% as_tibble()
  hit <- lists$sampleListId[lists$category == "all_cases_in_study"][1]
  if (is.na(hit) || length(hit) == 0) paste0(study_id, "_all") else hit
}

get_full_sample_roster <- function(sample_list_id) {
  res <- cbio_get(paste0("/sample-lists/", sample_list_id))
  ids <- res$sampleIds
  if (is.null(ids) || length(ids) == 0) {
    stop(sprintf("Sample list '%s' returned no sampleIds.", sample_list_id))
  }
  unique(as.character(unlist(ids)))
}

get_entrez_ids <- function(hugo_symbols) {
  res <- cbio_post("/genes/fetch", body = as.list(hugo_symbols),
                   query = list(geneIdType = "HUGO_GENE_SYMBOL")) %>% as_tibble()
  if (nrow(res) == 0) stop("Could not resolve Entrez IDs for: ", paste(hugo_symbols, collapse = ", "))
  res
}

attach_hugo_symbol <- function(df, gene_lookup) {
  if (nrow(df) == 0) return(df)
  if (!"entrezGeneId" %in% names(df)) {
    stop("Response is missing 'entrezGeneId' -- API response shape has changed.")
  }
  df <- df %>% select(-any_of(c("hugoGeneSymbol", "gene.hugoGeneSymbol")))
  lookup <- gene_lookup %>% select(entrezGeneId, hugoGeneSymbol) %>%
    mutate(entrezGeneId = as.integer(entrezGeneId))
  df %>% mutate(entrezGeneId = as.integer(entrezGeneId)) %>%
    left_join(lookup, by = "entrezGeneId")
}

# ---- Burden covariates ----
BURDEN_ATTR_CANDIDATES <- list(
  mutation_count          = c("MUTATION_COUNT", "TOTAL_MUTATION_COUNT", "TMB_NONSYNONYMOUS"),
  fraction_genome_altered = c("FRACTION_GENOME_ALTERED", "FRACTION_GENOME_CHANGED")
)

get_burden_covariates <- function(clinical) {
  out <- tibble(sampleId = unique(clinical$sampleId))
  found <- character(0)
  for (covar in names(BURDEN_ATTR_CANDIDATES)) {
    attr_id <- NULL
    for (cand in BURDEN_ATTR_CANDIDATES[[covar]]) {
      if (cand %in% clinical$clinicalAttributeId) { attr_id <- cand; break }
    }
    if (!is.null(attr_id)) {
      vals <- clinical %>% filter(clinicalAttributeId == attr_id) %>%
        transmute(sampleId, !!covar := suppressWarnings(as.numeric(value)))
      out <- out %>% left_join(vals, by = "sampleId")
      found <- c(found, sprintf("%s (from %s)", covar, attr_id))
    } else {
      out[[covar]] <- NA_real_
    }
  }
  attr(out, "found") <- found
  out
}

# ---- Pull mutation + CNA + expression data for the 2 genes ----
fetch_gene_pair_data <- function(gene1, gene2, study_id) {
  genes <- c(gene1, gene2)
  gene_lookup <- get_entrez_ids(genes)
  
  # get_entrez_ids only errors if NOTHING resolved -- it silently drops any
  # individual symbol the API didn't recognize. For the two target genes that
  # must be caught explicitly, or the pipeline will quietly analyze a
  # nonexistent gene as if it has zero alterations everywhere.
  resolved_symbols  <- toupper(gene_lookup$hugoGeneSymbol)
  requested_symbols <- toupper(genes)
  missing_symbols   <- setdiff(requested_symbols, resolved_symbols)
  if (length(missing_symbols) > 0) {
    stop(sprintf("Could not resolve the following gene symbol(s) -- check spelling: %s",
                 paste(missing_symbols, collapse = ", ")))
  }
  
  entrez_ids <- gene_lookup$entrezGeneId
  
  message(sprintf("[1/6] Resolved genes: %s",
                  paste(gene_lookup$hugoGeneSymbol, "=", gene_lookup$entrezGeneId, collapse = ", ")))
  # Check physical co-localization directly from cytoband, rather than
  # inferring it indirectly from CNA call patterns. This catches BOTH
  # arm-level co-deletion (e.g. PBRM1/BAP1 on 3p) and focal co-amplification
  # in a shared amplicon (e.g. CDK12/ERBB2 on 17q12) with one criterion.
  # cBioPortal's gene endpoints do NOT return cytoband or coordinates
  # (confirmed empirically -- only entrezGeneId/hugoGeneSymbol/type).
  # MyGene.info provides exact genomic coordinates, which is both more
  # reliable and more precise than a cytoband string: two genes sharing
  # an "arm" can be 50+ Mb apart, but a real shared CNA event (co-deletion
  # or co-amplification) puts them within a few Mb of each other.
  COLOCALIZATION_DISTANCE_BP <- 10000000  # 10 Mb -- generous for arm-level/amplicon events
  
  colocalization <- tryCatch({
    get_position <- function(entrez_id) {
      resp <- httr::GET(sprintf("https://mygene.info/v3/gene/%s", entrez_id),
                        query = list(fields = "genomic_pos"), httr::accept_json(),
                        httr::timeout(15))
      if (httr::http_error(resp)) return(NULL)
      info <- jsonlite::fromJSON(httr::content(resp, as = "text", encoding = "UTF-8"), flatten = TRUE)
      gp <- info$genomic_pos
      if (is.null(gp)) return(NULL)
      if (is.data.frame(gp)) gp <- as.list(gp[1, ])
      list(chr = as.character(gp$chr), start = as.numeric(gp$start), end = as.numeric(gp$end))
    }
    
    p1 <- get_position(gene_lookup$entrezGeneId[gene_lookup$hugoGeneSymbol == genes[1]][1])
    p2 <- get_position(gene_lookup$entrezGeneId[gene_lookup$hugoGeneSymbol == genes[2]][1])
    
    if (is.null(p1) || is.null(p2)) {
      list(same_arm = FALSE, gene1_pos = NA, gene2_pos = NA, distance_bp = NA)
    } else {
      same_chr <- identical(p1$chr, p2$chr)
      dist_bp <- if (same_chr) abs(p1$start - p2$start) else NA_real_
      close <- same_chr && !is.na(dist_bp) && dist_bp <= COLOCALIZATION_DISTANCE_BP
      list(same_arm = close,
           gene1_pos = sprintf("chr%s:%d", p1$chr, p1$start),
           gene2_pos = sprintf("chr%s:%d", p2$chr, p2$start),
           distance_bp = dist_bp)
    }
  }, error = function(e) {
    message(sprintf("    NOTE: MyGene.info co-localization lookup failed (%s) -- arm-level/amplicon flag defaults to FALSE for this pair.", conditionMessage(e)))
    list(same_arm = FALSE, gene1_pos = NA, gene2_pos = NA, distance_bp = NA)
  })
  
  if (isTRUE(colocalization$same_arm)) {
    message(sprintf("    !! CO-LOCALIZATION: %s and %s are %.1f Mb apart on the same chromosome.",
                    genes[1], genes[2], colocalization$distance_bp / 1e6))
  }
  profiles <- get_molecular_profiles(study_id)
  sample_list_id <- get_default_sample_list(study_id)
  message(sprintf("[2/6] Study '%s' has %d molecular profiles. Sample list: %s",
                  study_id, nrow(profiles), sample_list_id))
  
  full_sample_ids <- tryCatch(
    get_full_sample_roster(sample_list_id),
    error = function(e) {
      message(sprintf("    !! Could not fetch full sample roster (%s) -- falling back to inferring from mutation/CNA data.", conditionMessage(e)))
      NULL
    }
  )
  message(sprintf("    Sample roster: %s",
                  if (!is.null(full_sample_ids)) {
                    sprintf("%d samples confirmed from sample list '%s'", length(full_sample_ids), sample_list_id)
                  } else {
                    "UNAVAILABLE -- using fallback inference"
                  }))
  mut_id  <- profiles$molecularProfileId[profiles$molecularAlterationType == "MUTATION_EXTENDED"][1]
  expr_id <- profiles$molecularProfileId[profiles$molecularAlterationType == "MRNA_EXPRESSION" &
                                           grepl("zscore", profiles$molecularProfileId, ignore.case = TRUE) &
                                           !grepl("mirna", profiles$molecularProfileId, ignore.case = TRUE)][1]
  if (is.na(expr_id)) message("    !! No mRNA z-score profile found for this study (only miRNA or none) -- Terms 5/6 will be N/A")
  
  mut_df <- if (!is.na(mut_id)) {
    cbio_post(paste0("/molecular-profiles/", mut_id, "/mutations/fetch"),
              body = list(entrezGeneIds = as.list(entrez_ids), sampleListId = sample_list_id),
              query = list(projection = "SUMMARY")) %>% as_tibble() %>% attach_hugo_symbol(gene_lookup)
  } else tibble()
  message(sprintf("[4/6] Mutation rows returned: %d", nrow(mut_df)))
  
  cna_candidates <- profiles$molecularProfileId[profiles$molecularAlterationType == "COPY_NUMBER_ALTERATION" &
                                                  grepl("DISCRETE", profiles$datatype, ignore.case = TRUE)]
  cna_id <- NA
  cna_df <- tibble()
  for (candidate in cna_candidates) {
    trial <- cbio_post(paste0("/molecular-profiles/", candidate, "/molecular-data/fetch"),
                       body = list(entrezGeneIds = as.list(entrez_ids), sampleListId = sample_list_id),
                       query = list(projection = "SUMMARY")) %>% as_tibble() %>% attach_hugo_symbol(gene_lookup)
    covered_genes <- unique(trial$hugoGeneSymbol)
    if (all(toupper(genes) %in% toupper(covered_genes))) {
      cna_id <- candidate
      cna_df <- trial
      break
    } else {
      message(sprintf("    [CNA fallback] Profile '%s' missing gene(s): %s -- trying next candidate",
                      candidate, paste(setdiff(toupper(genes), toupper(covered_genes)), collapse=", ")))
    }
  }
  if (is.na(cna_id)) message("    !! No CNA profile fully covers both genes -- using last attempted (partial coverage)")
  
  message(sprintf("[3/6] Picked profiles -> mutations: %s | CNA: %s | expression z-score: %s",
                  mut_id, cna_id, expr_id))
  message(sprintf("    CNA rows returned: %d", nrow(cna_df)))
  if (nrow(cna_df) > 0) {
    bad_vals <- setdiff(unique(cna_df$value), c(-2,-1,0,1,2))
    if (length(bad_vals) > 0) {
      warning(sprintf("!! CNA profile '%s' has non-GISTIC values: %s -- amp/del thresholds will be WRONG",
                      cna_id, paste(head(sort(bad_vals), 10), collapse = ", ")))
    }
  }
  hallmark_symbols <- unique(unlist(Filter(
    Negate(is.null),
    lapply(c("HALLMARK_E2F_TARGETS", "HALLMARK_ANGIOGENESIS"), get_hallmark_geneset)
  )))
  expr_gene_lookup <- gene_lookup
  expr_entrez_ids  <- entrez_ids
  if (length(hallmark_symbols) > 0) {
    hallmark_lookup <- tryCatch(get_entrez_ids(hallmark_symbols), error = function(e) tibble())
    if (nrow(hallmark_lookup) > 0) {
      expr_gene_lookup <- bind_rows(gene_lookup, hallmark_lookup) %>% distinct(entrezGeneId, .keep_all = TRUE)
      expr_entrez_ids  <- unique(c(entrez_ids, hallmark_lookup$entrezGeneId))
    }
  }
  
  message(sprintf("    Expression fetch will include %d genes (2 target + %d Hallmark pathway genes)",
                  length(expr_entrez_ids), length(expr_entrez_ids) - length(entrez_ids)))
  
  expr_df <- if (!is.na(expr_id)) {
    cbio_post(paste0("/molecular-profiles/", expr_id, "/molecular-data/fetch"),
              body = list(entrezGeneIds = as.list(expr_entrez_ids), sampleListId = sample_list_id),
              query = list(projection = "SUMMARY")) %>% as_tibble() %>% attach_hugo_symbol(expr_gene_lookup)
  } else tibble()
  message(sprintf("    Expression (z-score) rows returned: %d", nrow(expr_df)))
  
  clinical_sample <- cbio_get(paste0("/studies/", study_id, "/clinical-data"),
                              query = list(clinicalDataType = "SAMPLE", projection = "SUMMARY")) %>% as_tibble()
  
  clinical_patient <- tryCatch(
    cbio_get(paste0("/studies/", study_id, "/clinical-data"),
             query = list(clinicalDataType = "PATIENT", projection = "SUMMARY")) %>% as_tibble(),
    error = function(e) tibble()
  )
  
  if (nrow(clinical_patient) > 0) {
    samples <- cbio_get(paste0("/studies/", study_id, "/samples")) %>% as_tibble() %>%
      select(sampleId, patientId)
    patient_mapped <- clinical_patient %>%
      inner_join(samples, by = "patientId") %>%
      select(-patientId)
    clinical <- bind_rows(clinical_sample, patient_mapped)
  } else {
    clinical <- clinical_sample
  }
  message(sprintf("[5/6] Clinical rows returned: %d (sample) + %d (patient-mapped)",
                  nrow(clinical_sample), nrow(clinical) - nrow(clinical_sample)))
  
  burden_df <- get_burden_covariates(clinical)
  found_burden <- attr(burden_df, "found")
  message(sprintf("    Burden covariates for adjustment: %s",
                  if (length(found_burden) > 0) paste(found_burden, collapse = ", ") else "NONE available"))
  
  list(mut_df = mut_df, cna_df = cna_df, expr_df = expr_df, clinical = clinical, burden_df = burden_df,
       gene_lookup = gene_lookup, mut_id = mut_id, cna_id = cna_id, expr_id = expr_id,
       sample_list_id = sample_list_id, full_sample_ids = full_sample_ids, colocalization = colocalization)
}

# ---- Build one row-per-sample status table for the 2 genes ----
# ---- Build one row-per-sample status table for the 2 genes ----
build_status_table <- function(fetched, gene1, gene2) {
  mut_df <- fetched$mut_df
  cna_df <- fetched$cna_df
  
  if (!is.null(fetched$full_sample_ids)) {
    all_samples <- fetched$full_sample_ids
    roster_source <- sprintf("study sample list '%s'", fetched$sample_list_id)
  } else {
    all_samples <- unique(c(mut_df$sampleId, cna_df$sampleId))
    roster_source <- "FALLBACK: inferred from mutation/CNA rows"
  }
  if (length(all_samples) == 0) stop("No samples returned for this gene pair / study combination.")
  
  status <- tibble(sampleId = all_samples) %>%
    mutate(
      geneA_mut = nrow(mut_df) > 0 & sampleId %in% mut_df$sampleId[mut_df$hugoGeneSymbol == gene1],
      geneB_mut = nrow(mut_df) > 0 & sampleId %in% mut_df$sampleId[mut_df$hugoGeneSymbol == gene2],
      geneA_del = nrow(cna_df) > 0 & sampleId %in% cna_df$sampleId[cna_df$hugoGeneSymbol == gene1 & cna_df$value <= -1],
      geneB_del = nrow(cna_df) > 0 & sampleId %in% cna_df$sampleId[cna_df$hugoGeneSymbol == gene2 & cna_df$value <= -1],
      geneA_amp = nrow(cna_df) > 0 & sampleId %in% cna_df$sampleId[cna_df$hugoGeneSymbol == gene1 & cna_df$value >= 1],
      geneB_amp = nrow(cna_df) > 0 & sampleId %in% cna_df$sampleId[cna_df$hugoGeneSymbol == gene2 & cna_df$value >= 1],
      geneA_deep_del = nrow(cna_df) > 0 & sampleId %in% cna_df$sampleId[cna_df$hugoGeneSymbol == gene1 & cna_df$value <= -2],
      geneB_deep_del = nrow(cna_df) > 0 & sampleId %in% cna_df$sampleId[cna_df$hugoGeneSymbol == gene2 & cna_df$value <= -2],
      # GISTIC value 2 = high-level/focal amplification, 1 = low-level/broad gain.
      # Mirrors deep_del's distinction between focal and arm-level events.
      geneA_deep_amp = nrow(cna_df) > 0 & sampleId %in% cna_df$sampleId[cna_df$hugoGeneSymbol == gene1 & cna_df$value >= 2],
      geneB_deep_amp = nrow(cna_df) > 0 & sampleId %in% cna_df$sampleId[cna_df$hugoGeneSymbol == gene2 & cna_df$value >= 2],
      geneA_altered = if (check_gene_role(gene1) == "tumor_suppressor") geneA_mut | geneA_del else geneA_mut | geneA_del | geneA_amp,
      geneB_altered = if (check_gene_role(gene2) == "tumor_suppressor") geneB_mut | geneB_del else geneB_mut | geneB_del | geneB_amp,
      either_altered = geneA_altered | geneB_altered,
      both_altered   = geneA_altered & geneB_altered,
      # STRICT now excludes shallow/broad events on BOTH sides (del and amp),
      # keeping only mutation, deep/focal deletion, and deep/focal amplification.
      geneA_altered_strict = geneA_mut | geneA_deep_del | geneA_deep_amp,
      geneB_altered_strict = geneB_mut | geneB_deep_del | geneB_deep_amp,
      either_altered_strict = geneA_altered_strict | geneB_altered_strict
    ) %>%
    left_join(fetched$burden_df, by = "sampleId")
  
  message(sprintf("[6/6] Status table built: %d samples (roster: %s) | %s altered: %d | %s altered: %d",
                  nrow(status), roster_source, gene1, sum(status$geneA_altered), gene2, sum(status$geneB_altered)))
  status
}
# ---- Detect arm-level co-deletion signature ----
detect_arm_level <- function(status, gene1, gene2) {
  check_one <- function(broad_col, focal_col, label, event_type) {
    broad_rate <- mean(status[[broad_col]])
    focal_rate <- mean(status[[focal_col]])
    flagged <- broad_rate > 0.5 && (focal_rate / max(broad_rate, 1e-9)) < 0.1
    list(gene = label, event_type = event_type,
         broad_rate = broad_rate, focal_rate = focal_rate, flagged = flagged,
         del_rate = broad_rate, deep_rate = focal_rate)
  }
  
  del_A <- check_one("geneA_del", "geneA_deep_del", gene1, "deletion")
  del_B <- check_one("geneB_del", "geneB_deep_del", gene2, "deletion")
  amp_A <- check_one("geneA_amp", "geneA_deep_amp", gene1, "amplification")
  amp_B <- check_one("geneB_amp", "geneB_deep_amp", gene2, "amplification")
  
  geneA <- if (del_A$flagged) del_A else if (amp_A$flagged) amp_A else del_A
  geneB <- if (del_B$flagged) del_B else if (amp_B$flagged) amp_B else del_B
  
  list(geneA = geneA, geneB = geneB)
}
# ==================== 2. ONE FUNCTION PER TERM ====================

safe_fisher <- function(status, alternative = "two.sided") {
  tab <- table(status$geneA_altered, status$geneB_altered)
  n_total <- nrow(status)
  size_note <- sprintf("[n=%d total | geneA altered: %d/%d | geneB altered: %d/%d]",
                       n_total, sum(status$geneA_altered), n_total, sum(status$geneB_altered), n_total)
  if (nrow(tab) < 2 || ncol(tab) < 2) {
    return(list(table = tab, odds_ratio = NA_real_, p_value = NA_real_, size_note = size_note))
  }
  test <- fisher.test(tab, alternative = alternative)
  list(table = tab, odds_ratio = unname(test$estimate), p_value = test$p.value, size_note = size_note)
}

burden_adjusted_association <- function(status, geneA_col = "geneA_altered", geneB_col = "geneB_altered") {
  covars <- c("mutation_count", "fraction_genome_altered")
  usable <- covars[vapply(covars, function(cv) {
    v <- status[[cv]]
    !is.null(v) && sum(!is.na(v)) >= 0.5 * nrow(status) && length(unique(na.omit(v))) > 1
  }, logical(1))]
  
  if (length(usable) == 0) return(list(available = FALSE, note = "No burden covariate available."))
  
  model_data <- status %>% select(all_of(c(geneA_col, geneB_col)), all_of(usable)) %>% tidyr::drop_na()
  if (nrow(model_data) < 20 || length(unique(model_data[[geneA_col]])) < 2 || length(unique(model_data[[geneB_col]])) < 2) {
    return(list(available = FALSE, note = sprintf("Only %d samples had complete burden data.", nrow(model_data))))
  }
  
  for (cv in usable) model_data[[cv]] <- as.numeric(scale(model_data[[cv]]))
  form <- as.formula(paste(geneB_col, "~", geneA_col, "+", paste(usable, collapse = " + ")))
  fit <- tryCatch(glm(form, data = model_data, family = binomial()), error = function(e) NULL)
  if (is.null(fit)) return(list(available = FALSE, note = "Adjusted logistic model failed to fit."))
  
  coefs <- summary(fit)$coefficients
  coef_row <- grep(paste0("^", geneA_col), rownames(coefs), value = TRUE)[1]
  if (is.na(coef_row)) return(list(available = FALSE, note = "Could not extract coefficient."))
  
  est <- coefs[coef_row, "Estimate"]
  se  <- coefs[coef_row, "Std. Error"]
  if (!is.finite(est) || !is.finite(se) || abs(est) > 10 || se > 10) {
    return(list(available = FALSE,
                note = "Adjusted model shows signs of separation (small/skewed sample) -- coefficient unreliable."))
  }
  
  z_stat <- est / se
  p_one_sided_less <- pnorm(z_stat)  # P(Z < z) -- tests OR < 1 (exclusivity direction)
  list(available = TRUE,
       odds_ratio = unname(exp(est)),
       p_value = unname(coefs[coef_row, "Pr(>|z|)"]),
       p_value_one_sided_less = unname(p_one_sided_less),
       covariates_used = usable,
       n_used = nrow(model_data))
}
term_combination_signal <- function(status) {
  raw <- safe_fisher(status)
  adj <- burden_adjusted_association(status)
  raw_txt <- if (!is.na(raw$p_value)) sprintf("Raw: OR=%.2f, p=%.4g", raw$odds_ratio, raw$p_value) else "Raw: not testable"
  
  if (isTRUE(adj$available)) {
    adj_txt <- sprintf("Burden-adjusted: OR=%.2f, p=%.4g (controlling for %s; n=%d)",
                       adj$odds_ratio, adj$p_value, paste(adj$covariates_used, collapse = " + "), adj$n_used)
    
    p_for_summary <- if (!is.na(raw$p_value) && adj$p_value < raw$p_value) raw$p_value else adj$p_value
    
    call <- if (!is.na(p_for_summary) && p_for_summary < 0.05 && adj$odds_ratio > 1) {
      if (!is.na(raw$p_value) && adj$p_value < raw$p_value) {
        "Co-occurring signal (raw only -- adjusted model showed suppression artifact)"
      } else {
        "Co-occurring / combination signal (holds after burden adjustment)"
      }
    } else "No significant combination signal"
    interp <- sprintf("%s. %s | %s %s", call, raw_txt, adj_txt, raw$size_note %||% "")
  } else {
    call <- if (!is.na(raw$p_value) && raw$odds_ratio > 1 && raw$p_value < 0.05) {
      "Co-occurring, UNADJUSTED result only -- interpret cautiously"
    } else "No significant combination signal (unadjusted)"
    interp <- sprintf("%s. %s | Adjusted test unavailable: %s %s", call, raw_txt, adj$note, raw$size_note %||% "")
    p_for_summary <- raw$p_value
  }
  
  verdict <- if (isTRUE(adj$available)) {
    if (!is.na(p_for_summary) && p_for_summary < 0.05 && adj$odds_ratio > 1) "SIGNIFICANT (co-occurring)" else "NOT SIGNIFICANT"
  } else if (grepl("separation", adj$note %||% "")) {
    "UNTRUSTWORTHY (too few samples)"
  } else if (!is.na(raw$p_value) && raw$odds_ratio > 1 && raw$p_value < 0.05) {
    "SIGNIFICANT (unadjusted only)"
  } else "NOT SIGNIFICANT"
  
  list(table = raw$table, odds_ratio_raw = raw$odds_ratio, p_value_raw = raw$p_value,
       p_value = p_for_summary, interpretation = interp, verdict = verdict)
}

term_mutual_exclusivity <- function(status, arm_flag, status_strict = NULL) {
  n_total <- nrow(status)
  
  # ---- Evidence 1: mutation-only. Unconfounded by CNA co-localization, but
  # blind to any gene whose primary alteration mechanism is CNA (amp/del). ----
  mut_tab <- table(status$geneA_mut, status$geneB_mut)
  size_note <- sprintf("[n=%d total | geneA mutated: %d/%d | geneB mutated: %d/%d]",
                       n_total, sum(status$geneA_mut), n_total, sum(status$geneB_mut), n_total)
  mut_test <- if (nrow(mut_tab) < 2 || ncol(mut_tab) < 2) {
    list(p_value = NA_real_, odds_ratio = NA_real_)
  } else {
    ft <- fisher.test(mut_tab, alternative = "less")
    list(p_value = ft$p.value, odds_ratio = unname(ft$estimate))
  }
  mut_adj <- burden_adjusted_association(status, "geneA_mut", "geneB_mut")
  mut_adj_p_use <- if (isTRUE(mut_adj$available) && !is.null(mut_adj$p_value_one_sided_less)) mut_adj$p_value_one_sided_less else mut_adj$p_value
  mut_p <- if (isTRUE(mut_adj$available)) {
    if (!is.na(mut_test$p_value) && mut_adj_p_use < mut_test$p_value) mut_test$p_value else mut_adj_p_use
  } else mut_test$p_value
  mut_or <- if (isTRUE(mut_adj$available)) mut_adj$odds_ratio else mut_test$odds_ratio
  mut_txt <- if (!is.na(mut_test$p_value)) {
    sprintf("Mutation-only: raw OR=%.2g p=%.4g%s", mut_test$odds_ratio, mut_test$p_value,
            if (isTRUE(mut_adj$available)) sprintf(" | adjusted OR=%.2f p=%.4g (n=%d)", mut_adj$odds_ratio, mut_adj$p_value, mut_adj$n_used) else "")
  } else "Mutation-only: not testable (insufficient mutation variation)"
  
  # ---- Evidence 2: altered-level (mut + CNA). THIS is the test with actual
  # power for CNA-driven pathway members (e.g. CDK4 amp, RB1/CDKN2A del) that
  # mutation-only structurally misses. Reuses the exact same machinery Term 2
  # already uses (safe_fisher / burden_adjusted_association), just consumed
  # here too instead of only there. ----
  broad_raw <- safe_fisher(status, alternative = "less")
  broad_adj <- burden_adjusted_association(status)  # default cols = geneA_altered/geneB_altered
  
  use_strict <- !is.null(status_strict) && (isTRUE(arm_flag$geneA$flagged) || isTRUE(arm_flag$geneB$flagged))
  strict_raw <- if (use_strict) safe_fisher(status_strict, alternative = "less") else NULL
  strict_adj <- if (use_strict) burden_adjusted_association(status_strict) else NULL
  
  pick_tier <- function(raw, adj, label) {
    if (is.null(raw) || is.na(raw$p_value)) return(NULL)
    if (isTRUE(adj$available)) {
      # same suppression-artifact protection convention as the rest of the
      # pipeline: adjustment is never allowed to manufacture significance
      # raw didn't already support.
      adj_p_use <- if (!is.null(adj$p_value_one_sided_less)) adj$p_value_one_sided_less else adj$p_value
      p <- if (!is.na(raw$p_value) && adj_p_use < raw$p_value) raw$p_value else adj_p_use
      or <- adj$odds_ratio
      txt <- sprintf("%s: raw OR=%.2f p=%.4g | adjusted OR=%.2f p=%.4g (n=%d)",
                     label, raw$odds_ratio, raw$p_value, adj$odds_ratio, adj$p_value, adj$n_used)
    } else {
      p <- raw$p_value
      or <- raw$odds_ratio
      txt <- sprintf("%s: raw OR=%.2f p=%.4g (adjustment unavailable: %s)", label, raw$odds_ratio, raw$p_value, adj$note %||% "n/a")
    }
    list(p = p, or = or, txt = txt)
  }
  
  # If the pair is arm-flagged (shared arm-level deletion or shared amplicon),
  # the STRICT definition (mut + deep/focal CNA only) is the trustworthy one --
  # the broad definition would just be re-detecting the shared CNA event, not
  # gene-specific exclusivity. Fall back to broad only if strict has no power.
  altered_tier <- if (use_strict) {
    strict_result <- pick_tier(strict_raw, strict_adj,
                               "STRICT altered-level (mut+deep/focal CNA, arm-level/amplicon events excluded)")
    if (!is.null(strict_result)) strict_result else pick_tier(broad_raw, broad_adj, "Altered-level (broad -- strict fallback had no power)")
  } else {
    pick_tier(broad_raw, broad_adj, "Altered-level (mut+CNA)")
  }
  
  # ---- Final verdict. Altered-level is now PRIMARY (it's the test with power
  # for CNA-driven genes); mutation-only is kept as a transparent secondary
  # line of evidence and flagged explicitly when the two disagree, instead of
  # being silently discarded either way. ----
  if (!is.null(altered_tier) && !is.na(altered_tier$p) && !is.na(altered_tier$or)) {
    p_for_summary  <- altered_tier$p
    or_for_summary <- altered_tier$or
    evidence_used  <- "altered-level"
  } else if (!is.na(mut_p) && !is.na(mut_or)) {
    p_for_summary  <- mut_p
    or_for_summary <- mut_or
    evidence_used  <- "mutation-only (altered-level not testable)"
  } else {
    p_for_summary  <- NA_real_
    or_for_summary <- NA_real_
    evidence_used  <- "none"
  }
  
  call <- if (is.na(p_for_summary) || is.na(or_for_summary)) {
    "No exclusivity detected (insufficient data at any level)"
  } else if (or_for_summary < 1 && p_for_summary < 0.05) {
    sprintf("Mutually exclusive pattern (%s)", evidence_used)
  } else if (or_for_summary > 1 && p_for_summary < 0.05) {
    sprintf("Co-occurring, not exclusive (%s)", evidence_used)
  } else {
    sprintf("No exclusivity detected (%s)", evidence_used)
  }
  
  discord <- !is.na(mut_or) && !is.null(altered_tier) && !is.na(altered_tier$or) &&
    ((mut_or < 1) != (altered_tier$or < 1))
  discord_note <- if (isTRUE(discord)) {
    " -- NOTE: mutation-only and altered-level tests disagree on direction; expected when a gene's dominant alteration mechanism is CNA rather than point mutation -- trusting the altered-level result."
  } else ""
  
  near_miss_note <- if (!is.na(p_for_summary) && !is.na(or_for_summary) &&
                        or_for_summary < 1 && p_for_summary >= 0.05 && p_for_summary < 0.10) {
    " [NOTE: direction consistent with exclusivity (OR<1); did not reach significance at this sample size.]"
  } else ""
  
  scope_note <- " [SCOPE: pairwise test only -- cannot detect exclusivity that requires pooling >2 genes in a shared pathway module (e.g. MEMo-style module tests); a 'NOT SIGNIFICANT' here does not rule out pathway-level exclusivity documented elsewhere in the literature.]"
  interp <- sprintf("%s. %s | %s%s %s%s%s",
                    call, mut_txt,
                    if (!is.null(altered_tier)) altered_tier$txt else "Altered-level: not testable",
                    discord_note, size_note, scope_note, near_miss_note)
  
  verdict <- if (is.na(p_for_summary) || is.na(or_for_summary)) {
    "N/A (not testable)"
  } else if (or_for_summary < 1 && p_for_summary < 0.05) {
    "MUTUALLY EXCLUSIVE"
  } else if (or_for_summary > 1 && p_for_summary < 0.05) {
    "SIGNIFICANT (co-occurring)"
  } else "NOT SIGNIFICANT"
  
  list(
    p_value = p_for_summary, odds_ratio = or_for_summary, evidence_used = evidence_used,
    p_value_mutation = mut_p, odds_ratio_mutation = mut_or,
    interpretation = interp, verdict = verdict
  )
}
# UPDATE 12: Added burden adjustment via Likelihood Ratio Test (LRT)
# UPDATE 14: During validation on brca_tcga_pub, PAM50_SUBTYPE (the actual
# molecular-subtype call) was found AFTER PAM50_CALL in this list, so on
# cohorts carrying both attributes the loop picked the less specific one
# first without erroring. PAM50_SUBTYPE is now checked before PAM50_CALL.
# (PAM50_CALL is kept as a fallback for cohorts that only have that field.)
term_subtype_specificity <- function(status, clinical, subtype_attr = NULL, test_col = "either_altered"){
  possible_attrs <- c(subtype_attr, "SUBTYPE", "HISTOLOGICAL_SUBTYPE",
                      "BREAST_CANCER_SUBTYPE", "PAM50_SUBTYPE", "PAM50_CALL", "GBM_SUBTYPE",
                      "LUNG_CANCER_SUBTYPE", "COADREAD_SUBTYPE")
  possible_attrs <- possible_attrs[!is.na(possible_attrs) & possible_attrs != ""]
  
  found_attr <- NULL
  for (attr_name in possible_attrs) {
    if (attr_name %in% clinical$clinicalAttributeId) {
      vals <- clinical$value[clinical$clinicalAttributeId == attr_name]
      if (length(unique(na.omit(vals))) >= 2) {
        found_attr <- attr_name
        break
      }
    }
  }
  
  if (is.null(found_attr)) {
    all_attrs <- unique(clinical$clinicalAttributeId)
    stage_matches <- all_attrs[grepl("STAGE", all_attrs, ignore.case = TRUE)]
    grade_matches <- all_attrs[grepl("GRADE", all_attrs, ignore.case = TRUE)]
    for (attr in c(stage_matches, grade_matches)) {
      vals <- clinical$value[clinical$clinicalAttributeId == attr]
      if (length(unique(na.omit(vals))) >= 2) {
        found_attr <- attr
        message(sprintf("    [Subtype Fallback] Using clinical surrogate: '%s'", attr))
        break
      }
    }
  }
  
  if (is.null(found_attr)) {
    avail <- clinical %>%
      group_by(clinicalAttributeId) %>%
      summarise(n_unique = n_distinct(value), .groups = "drop") %>%
      filter(n_unique >= 2 & n_unique <= 20) %>%
      pull(clinicalAttributeId)
    msg <- if (length(avail) > 0) {
      sprintf("No subtype attribute found. Available categorical attributes: %s", paste(head(avail, 10), collapse = ", "))
    } else "No subtype attribute found and no suitable categorical attributes available"
    return(list(result = msg))
  }
  
  subtype_wide <- clinical %>%
    filter(clinicalAttributeId == found_attr) %>%
    select(sampleId, subtype = value)
  
  # Join status (which includes burden covariates) with the subtype data
  merged <- status %>% left_join(subtype_wide, by = "sampleId")
  tab <- table(merged$subtype, merged[[test_col]])
  
  n_total <- nrow(merged)
  n_altered <- sum(merged[[test_col]], na.rm = TRUE)
  n_unaltered <- n_total - n_altered
  size_note <- sprintf("[n=%d total | either-gene altered: %d | unaltered: %d | categories: %d]",
                       n_total, n_altered, n_unaltered, nrow(tab))
  
  if (nrow(tab) < 2 || ncol(tab) < 2) {
    return(list(table = tab, p_value = NA_real_,
                interpretation = paste0("Insufficient variation to test. ", size_note)))
  }
  
  # Raw Chi-square test
  test <- tryCatch(chisq.test(tab), error = function(e) NULL)
  p_raw <- if (!is.null(test)) test$p.value else NA_real_
  
  # UPDATE 12: Burden-adjusted LRT
  p_adj <- NA_real_
  adj_note <- ""
  burden_covs <- c("mutation_count", "fraction_genome_altered")
  usable_covs <- burden_covs[vapply(burden_covs, function(cv) {
    v <- merged[[cv]]
    !is.null(v) && sum(!is.na(v)) >= 0.5 * nrow(merged)
  }, logical(1))]
  
  if (length(usable_covs) > 0) {
    model_data <- merged %>% select(all_of(test_col), subtype, all_of(usable_covs)) %>% drop_na()
    if (nrow(model_data) > 30 && length(unique(model_data$subtype)) >= 2) {
      for (cv in usable_covs) model_data[[cv]] <- as.numeric(scale(model_data[[cv]]))
      adj_result <- tryCatch({
        null_form <- as.formula(paste(test_col, "~", paste(usable_covs, collapse = " + ")))
        fit_null <- glm(null_form, data = model_data, family = binomial())
        full_form <- as.formula(paste(test_col, "~ subtype +", paste(usable_covs, collapse = " + ")))
        fit_full <- glm(full_form, data = model_data, family = binomial())
        lrt <- anova(fit_null, fit_full, test = "LRT")
        p <- lrt$`Pr(>Chi)`[2]
        list(p_adj = p, note = sprintf(" (burden-adjusted p=%.4g)", p))
      }, error = function(e) {
        list(p_adj = NA_real_, note = sprintf(" (adjusted model failed to fit: %s)", conditionMessage(e)))
      })
      p_adj <- adj_result$p_adj
      adj_note <- adj_result$note
    } else {
      adj_note <- " (insufficient data for adjustment)"
    }
  } else {
    adj_note <- " (no burden covariates available)"
  }
  
  # Suppression artifact protection
  p_for_summary <- if (!is.na(p_adj) && !is.na(p_raw)) max(p_raw, p_adj) else p_raw
  top_subtype_note <- ""
  if (!is.na(p_for_summary) && p_for_summary < 0.05 && nrow(tab) >= 2 && ncol(tab) == 2) {
    subtype_n <- rowSums(tab)
    reliable <- subtype_n >= 10
    if (any(reliable)) {
      altered_col <- if ("TRUE" %in% colnames(tab)) "TRUE" else colnames(tab)[2]
      prop_altered_by_subtype <- (tab[, altered_col] / rowSums(tab))[reliable]
      top_subtype <- names(which.max(prop_altered_by_subtype))
      top_pct <- round(100 * max(prop_altered_by_subtype), 1)
      top_n <- subtype_n[top_subtype]
      top_subtype_note <- sprintf(" | Most-associated subtype: '%s' (%.1f%% altered, n=%d in this subtype)", top_subtype, top_pct, top_n)
    } else {
      top_subtype_note <- " | (all subtype categories too small (n<10) for reliable most-associated-subtype reporting)"
    }
  }
  interp <- if (is.na(p_for_summary)) {
    paste0("Chi-squared test could not be computed. ", size_note)
  } else if (p_for_summary < 0.05) {
    sprintf("Significant subtype association (p=%.4g)%s, %d subtypes tested %s%s", p_for_summary, adj_note, nrow(tab), size_note, top_subtype_note)
  } else {
    sprintf("No significant subtype association (p=%.4g)%s, %d subtypes tested %s", p_for_summary, adj_note, nrow(tab), size_note)
  }
  verdict <- if (is.na(p_for_summary)) "N/A" else if (p_for_summary < 0.05) "SIGNIFICANT" else "NOT SIGNIFICANT"
  
  list(table = tab, p_value = p_for_summary, interpretation = interp, verdict = verdict)
}

# UPDATE 13: Added co-localization artifact flag
term_double_hit <- function(status, gene1, gene2, arm_flag) {
  n <- nrow(status)
  rateA <- sum(status$geneA_mut & status$geneA_del) / n
  rateB <- sum(status$geneB_mut & status$geneB_del) / n
  deepA <- sum(status$geneA_mut & status$geneA_deep_del) / n
  deepB <- sum(status$geneB_mut & status$geneB_deep_del) / n
  any_altA <- mean(status$geneA_altered)
  any_altB <- mean(status$geneB_altered)
  
  # Use the same arm-level detector as the rest of the pipeline (relative,
  # deep-vs-shallow ratio based) instead of a separate absolute-difference
  # heuristic, so Term 4's note agrees with the ARM-LEVEL CO-DELETION FLAG.
  note_A <- if (isTRUE(arm_flag$geneA$flagged)) {
    " NOTE: Double-hit rate for Gene A mirrors its mutation rate--likely arm-level co-deletion artifact rather than independent biallelic inactivation."
  } else ""
  note_B <- if (isTRUE(arm_flag$geneB$flagged)) {
    " NOTE: Same artifact likely affects Gene B."
  } else ""
  
  interp <- sprintf(
    paste0(
      "[NOTE: calculated separately per gene -- NOT a comparison between the two genes] ",
      "n=%d total samples. %s: %.1f%% double-hit (mut+LOH/del), %.1f%% deep-del-only, %.1f%% any alteration | ",
      "%s: %.1f%% double-hit, %.1f%% deep-del-only, %.1f%% any alteration. %s%s",
      "Caution: mutation+CNA data only -- promoter methylation silencing (a known third ",
      "inactivation route for tumor suppressors such as VHL) is not captured by this pipeline."
    ),
    n, gene1, 100 * rateA, 100 * deepA, 100 * any_altA,
    gene2, 100 * rateB, 100 * deepB, 100 * any_altB, note_A, note_B)
  
  verdict <- if (note_A != "" || note_B != "") "LIKELY FALSE (chromosome position)" else "DESCRIPTIVE ONLY"
  
  list(
    geneA_double_hit_rate = rateA, geneB_double_hit_rate = rateB,
    geneA_double_hit_rate_deep_del_only = deepA, geneB_double_hit_rate_deep_del_only = deepB,
    geneA_any_alteration_rate = any_altA, geneB_any_alteration_rate = any_altB,
    interpretation = interp, verdict = verdict
  )
}

# ---- Hallmark gene-set fetch: cached + retried ----
.hallmark_cache_dir <- file.path(tools::R_user_dir("gene_pair_app", "cache"))
if (!dir.exists(.hallmark_cache_dir)) dir.create(.hallmark_cache_dir, recursive = TRUE, showWarnings = FALSE)

get_hallmark_geneset <- function(hallmark_name, timeout_sec = 60, max_tries = 3) {
  cache_file <- file.path(.hallmark_cache_dir, paste0(hallmark_name, ".rds"))
  if (file.exists(cache_file)) return(readRDS(cache_file))
  
  old_timeout <- getOption("timeout")
  options(timeout = max(timeout_sec, old_timeout %||% 60))
  on.exit(options(timeout = old_timeout), add = TRUE)
  
  genes <- NULL
  for (attempt in seq_len(max_tries)) {
    genes <- tryCatch({
      msigdbr(species = "Homo sapiens", category = "H") %>%
        filter(gs_name == hallmark_name) %>% pull(gene_symbol)
    }, error = function(e) {
      message(sprintf("    !! [%s] MSigDB fetch attempt %d/%d failed: %s", hallmark_name, attempt, max_tries, conditionMessage(e)))
      NULL
    })
    if (!is.null(genes) && length(genes) > 0) break
    if (attempt < max_tries) Sys.sleep(2 * attempt)
  }
  if (!is.null(genes) && length(genes) > 0) saveRDS(genes, cache_file)
  genes
}

burden_adjusted_score_test <- function(merged, altered_col) {
  covars <- c("mutation_count", "fraction_genome_altered")
  usable <- covars[vapply(covars, function(cv) {
    v <- merged[[cv]]
    !is.null(v) && sum(!is.na(v)) >= 0.5 * nrow(merged) && length(unique(na.omit(v))) > 1
  }, logical(1))]
  
  if (length(usable) == 0) return(list(available = FALSE, note = "No burden covariate available."))
  
  model_data <- merged %>% select(score, all_of(altered_col), all_of(usable)) %>% tidyr::drop_na()
  if (nrow(model_data) < 20 || length(unique(model_data[[altered_col]])) < 2) {
    return(list(available = FALSE, note = "Too few complete-burden-data samples."))
  }
  
  for (cv in usable) model_data[[cv]] <- as.numeric(scale(model_data[[cv]]))
  form <- as.formula(paste("score ~", altered_col, "+", paste(usable, collapse = " + ")))
  fit <- tryCatch(lm(form, data = model_data), error = function(e) NULL)
  if (is.null(fit)) return(list(available = FALSE, note = "Adjusted linear model failed to fit."))
  
  coefs <- summary(fit)$coefficients
  coef_row <- grep(paste0("^", altered_col), rownames(coefs), value = TRUE)[1]
  if (is.na(coef_row)) return(list(available = FALSE, note = "Could not extract coefficient."))
  
  est <- coefs[coef_row, "Estimate"]
  se  <- coefs[coef_row, "Std. Error"]
  if (!is.finite(est) || !is.finite(se)) {
    return(list(available = FALSE, note = "Adjusted model produced a non-finite coefficient -- unreliable."))
  }
  
  list(available = TRUE,
       estimate = unname(est),
       p_value = unname(coefs[coef_row, "Pr(>|t|)"]),
       covariates_used = usable,
       n_used = nrow(model_data))
}
score_and_test <- function(expr_df, status, hallmark_name, altered_col) {
  if (nrow(expr_df) == 0) return(list(result = "Expression (z-score) data not available for this study", verdict = "N/A"))
  
  gene_set <- get_hallmark_geneset(hallmark_name)
  if (is.null(gene_set) || length(gene_set) == 0) {
    return(list(result = sprintf("Could not retrieve the '%s' gene set from MSigDB.", hallmark_name), verdict = "N/A"))
  }
  
  matched <- expr_df %>% filter(hugoGeneSymbol %in% gene_set)
  if (length(unique(matched$hugoGeneSymbol)) < 3) {
    return(list(result = "Too few gene-set genes matched in this study's expression data", verdict = "N/A"))
  }
  
  scores <- matched %>%
    group_by(sampleId) %>%
    summarise(score = mean(value, na.rm = TRUE), n_genes = n_distinct(hugoGeneSymbol), .groups = "drop")
  
  merged <- status %>% inner_join(scores, by = "sampleId")
  
  n_altered_scored <- sum(merged[[altered_col]], na.rm = TRUE)
  n_unaltered_scored <- nrow(merged) - n_altered_scored
  size_note <- sprintf("[n=%d scored | altered: %d | unaltered: %d]", nrow(merged), n_altered_scored, n_unaltered_scored)
  
  if (length(unique(merged[[altered_col]])) < 2) {
    return(list(scores = scores, n_samples_scored = nrow(scores), p_value = NA_real_,
                interpretation = paste0("Insufficient variation to test. ", size_note), verdict = "N/A"))
  }
  
  wtest <- tryCatch(wilcox.test(score ~ get(altered_col), data = merged), error = function(e) NULL)
  p_raw <- if (!is.null(wtest)) wtest$p.value else NA_real_
  
  med_altered <- median(merged$score[merged[[altered_col]]], na.rm = TRUE)
  med_unaltered <- median(merged$score[!merged[[altered_col]]], na.rm = TRUE)
  direction <- if (is.na(med_altered) || is.na(med_unaltered)) "" else
    if (med_altered > med_unaltered) "higher" else "lower"
  
  adj <- burden_adjusted_score_test(merged, altered_col)
  
  if (isTRUE(adj$available)) {
    adj_txt <- sprintf("Burden-adjusted: p=%.4g (controlling for %s; n=%d)",
                       adj$p_value, paste(adj$covariates_used, collapse = " + "), adj$n_used)
    
    p_for_summary <- if (!is.na(p_raw) && adj$p_value < p_raw) p_raw else adj$p_value
    
    call <- if (!is.na(p_for_summary) && p_for_summary < 0.05) {
      if (!is.na(p_raw) && adj$p_value < p_raw) {
        sprintf("Significant (%s in altered group, raw only -- adjusted model showed suppression artifact)", direction)
      } else {
        sprintf("Significant (%s in altered group, holds after burden adjustment)", if (adj$estimate > 0) "higher" else "lower")
      }
    } else "No significant difference after burden adjustment"
    interp <- sprintf("%s. Raw Wilcoxon: p=%.4g (%s in altered group) | %s %s",
                      call, p_raw, direction, adj_txt, size_note)
  } else {
    interp <- if (is.na(p_raw)) {
      paste0("Could not compute test (check group sizes). ", size_note)
    } else if (p_raw < 0.05) {
      sprintf("Significant difference (%s in altered group), UNADJUSTED result only -- interpret cautiously. Adjusted test unavailable: %s %s",
              direction, adj$note, size_note)
    } else {
      sprintf("No significant difference between altered/unaltered groups %s. Adjusted test unavailable: %s", size_note, adj$note)
    }
    p_for_summary <- p_raw
  }
  
  verdict <- if (is.na(p_for_summary)) "N/A" else if (p_for_summary < 0.05) sprintf("SIGNIFICANT (%s)", toupper(direction)) else "NOT SIGNIFICANT"
  
  list(scores = scores, n_samples_scored = nrow(scores), p_value = p_for_summary,
       p_value_raw = p_raw, interpretation = interp, verdict = verdict)
}

term_cell_cycle <- function(fetched, status) {
  n <- nrow(status)
  n_alt <- sum(status$either_altered)
  if (n_alt / n > 0.85 || (n - n_alt) < 30) {
    message(sprintf("    [Low contrast] %.1f%% altered, %d unaltered -- pathway terms may be uninformative", 100 * n_alt / n, n - n_alt))
  }
  score_and_test(fetched$expr_df, status, "HALLMARK_E2F_TARGETS", "either_altered")
}

term_angiogenesis <- function(fetched, status) {
  n <- nrow(status)
  n_alt <- sum(status$either_altered)
  if (n_alt / n > 0.85 || (n - n_alt) < 30) {
    message(sprintf("    [Low contrast] %.1f%% altered, %d unaltered -- pathway terms may be uninformative", 100 * n_alt / n, n - n_alt))
  }
  score_and_test(fetched$expr_df, status, "HALLMARK_ANGIOGENESIS", "either_altered")
}

# ==================== 3. MAIN WRAPPER FUNCTION ====================

safe_term <- function(expr, label) {
  tryCatch(expr, error = function(e) list(result = sprintf("Term failed: %s", conditionMessage(e))))
}
# ==================== SINGLE-GENE ANALYSIS (Terms 1, 4, 5, 6 only) ====================

fetch_single_gene_data <- function(gene1, study_id) {
  gene_lookup <- get_entrez_ids(c(gene1))
  if (nrow(gene_lookup) == 0 || !toupper(gene1) %in% toupper(gene_lookup$hugoGeneSymbol)) {
    stop(sprintf("Could not resolve gene symbol: %s", gene1))
  }
  entrez_ids <- gene_lookup$entrezGeneId
  
  profiles <- get_molecular_profiles(study_id)
  sample_list_id <- get_default_sample_list(study_id)
  full_sample_ids <- tryCatch(get_full_sample_roster(sample_list_id), error = function(e) NULL)
  
  mut_id  <- profiles$molecularProfileId[profiles$molecularAlterationType == "MUTATION_EXTENDED"][1]
  expr_id <- profiles$molecularProfileId[profiles$molecularAlterationType == "MRNA_EXPRESSION" &
                                           grepl("zscore", profiles$molecularProfileId, ignore.case = TRUE) &
                                           !grepl("mirna", profiles$molecularProfileId, ignore.case = TRUE)][1]
  
  mut_df <- if (!is.na(mut_id)) {
    cbio_post(paste0("/molecular-profiles/", mut_id, "/mutations/fetch"),
              body = list(entrezGeneIds = as.list(entrez_ids), sampleListId = sample_list_id),
              query = list(projection = "SUMMARY")) %>% as_tibble() %>% attach_hugo_symbol(gene_lookup)
  } else tibble()
  
  cna_candidates <- profiles$molecularProfileId[profiles$molecularAlterationType == "COPY_NUMBER_ALTERATION" &
                                                  grepl("DISCRETE", profiles$datatype, ignore.case = TRUE)]
  cna_id <- NA; cna_df <- tibble()
  for (candidate in cna_candidates) {
    trial <- cbio_post(paste0("/molecular-profiles/", candidate, "/molecular-data/fetch"),
                       body = list(entrezGeneIds = as.list(entrez_ids), sampleListId = sample_list_id),
                       query = list(projection = "SUMMARY")) %>% as_tibble() %>% attach_hugo_symbol(gene_lookup)
    if (toupper(gene1) %in% toupper(unique(trial$hugoGeneSymbol))) { cna_id <- candidate; cna_df <- trial; break }
  }
  message(sprintf("[single-gene] mutations: %s | CNA: %s | expression: %s", mut_id, cna_id, expr_id))
  
  hallmark_symbols <- unique(unlist(Filter(Negate(is.null),
                                           lapply(c("HALLMARK_E2F_TARGETS", "HALLMARK_ANGIOGENESIS"), get_hallmark_geneset))))
  expr_gene_lookup <- gene_lookup; expr_entrez_ids <- entrez_ids
  if (length(hallmark_symbols) > 0) {
    hallmark_lookup <- tryCatch(get_entrez_ids(hallmark_symbols), error = function(e) tibble())
    if (nrow(hallmark_lookup) > 0) {
      expr_gene_lookup <- bind_rows(gene_lookup, hallmark_lookup) %>% distinct(entrezGeneId, .keep_all = TRUE)
      expr_entrez_ids  <- unique(c(entrez_ids, hallmark_lookup$entrezGeneId))
    }
  }
  expr_df <- if (!is.na(expr_id)) {
    cbio_post(paste0("/molecular-profiles/", expr_id, "/molecular-data/fetch"),
              body = list(entrezGeneIds = as.list(expr_entrez_ids), sampleListId = sample_list_id),
              query = list(projection = "SUMMARY")) %>% as_tibble() %>% attach_hugo_symbol(expr_gene_lookup)
  } else tibble()
  
  clinical_sample <- cbio_get(paste0("/studies/", study_id, "/clinical-data"),
                              query = list(clinicalDataType = "SAMPLE", projection = "SUMMARY")) %>% as_tibble()
  clinical_patient <- tryCatch(cbio_get(paste0("/studies/", study_id, "/clinical-data"),
                                        query = list(clinicalDataType = "PATIENT", projection = "SUMMARY")) %>% as_tibble(), error = function(e) tibble())
  clinical <- if (nrow(clinical_patient) > 0) {
    samples <- cbio_get(paste0("/studies/", study_id, "/samples")) %>% as_tibble() %>% select(sampleId, patientId)
    bind_rows(clinical_sample, clinical_patient %>% inner_join(samples, by = "patientId") %>% select(-patientId))
  } else clinical_sample
  
  burden_df <- get_burden_covariates(clinical)
  
  list(mut_df = mut_df, cna_df = cna_df, expr_df = expr_df, clinical = clinical, burden_df = burden_df,
       gene_lookup = gene_lookup, sample_list_id = sample_list_id, full_sample_ids = full_sample_ids)
}

build_single_gene_status_table <- function(fetched, gene1) {
  mut_df <- fetched$mut_df; cna_df <- fetched$cna_df
  all_samples <- if (!is.null(fetched$full_sample_ids)) fetched$full_sample_ids else unique(c(mut_df$sampleId, cna_df$sampleId))
  if (length(all_samples) == 0) stop("No samples returned.")
  
  status <- tibble(sampleId = all_samples) %>%
    mutate(
      geneA_mut = nrow(mut_df) > 0 & sampleId %in% mut_df$sampleId[mut_df$hugoGeneSymbol == gene1],
      geneA_del = nrow(cna_df) > 0 & sampleId %in% cna_df$sampleId[cna_df$hugoGeneSymbol == gene1 & cna_df$value <= -1],
      geneA_amp = nrow(cna_df) > 0 & sampleId %in% cna_df$sampleId[cna_df$hugoGeneSymbol == gene1 & cna_df$value >= 1],
      geneA_deep_del = nrow(cna_df) > 0 & sampleId %in% cna_df$sampleId[cna_df$hugoGeneSymbol == gene1 & cna_df$value <= -2],
      geneA_altered = if (check_gene_role(gene1) == "tumor_suppressor") geneA_mut | geneA_del else geneA_mut | geneA_del | geneA_amp,
      either_altered = geneA_altered
    ) %>%
    left_join(fetched$burden_df, by = "sampleId")
  
  message(sprintf("[single-gene] status table: %d samples | %s altered: %d", nrow(status), gene1, sum(status$geneA_altered)))
  status
}

analyze_single_gene <- function(gene1, study_id) {
  fetched <- fetch_single_gene_data(gene1, study_id)
  status  <- build_single_gene_status_table(fetched, gene1)
  
  subtype_res <- safe_term(term_subtype_specificity(status, fetched$clinical, test_col = "geneA_altered"), "subtype")
  cell_cycle_res <- safe_term(score_and_test(fetched$expr_df, status, "HALLMARK_E2F_TARGETS", "geneA_altered"), "cell_cycle")
  angio_res <- safe_term(score_and_test(fetched$expr_df, status, "HALLMARK_ANGIOGENESIS", "geneA_altered"), "angiogenesis")
  
  n <- nrow(status)
  rateA <- sum(status$geneA_mut & status$geneA_del) / n
  deepA <- sum(status$geneA_mut & status$geneA_deep_del) / n
  any_altA <- mean(status$geneA_altered)
  double_hit_res <- list(
    interpretation = sprintf("n=%d. %s: %.1f%% double-hit (mut+LOH/del), %.1f%% deep-del-only, %.1f%% any alteration.",
                             n, gene1, 100*rateA, 100*deepA, 100*any_altA),
    verdict = "DESCRIPTIVE ONLY", geneA_double_hit_rate = rateA
  )
  
  cat(sprintf("\n===== SINGLE-GENE RESULTS: %s | Study: %s =====\n", gene1, study_id))
  cat(sprintf("1. Subtype specificity >>> %s\n   %s\n\n", subtype_res$verdict %||% "N/A", subtype_res$interpretation %||% subtype_res$result %||% ""))
  cat(sprintf("4. Double hit >>> %s\n   %s\n\n", double_hit_res$verdict, double_hit_res$interpretation))
  cat(sprintf("5. Cell-cycle >>> %s\n   %s\n\n", cell_cycle_res$verdict %||% "N/A", cell_cycle_res$interpretation %||% cell_cycle_res$result %||% ""))
  cat(sprintf("6. Angiogenesis >>> %s\n   %s\n\n", angio_res$verdict %||% "N/A", angio_res$interpretation %||% angio_res$result %||% ""))
  
  list(subtype = subtype_res, double_hit = double_hit_res, cell_cycle = cell_cycle_res, angiogenesis = angio_res)
}
analyze_gene_pair <- function(gene1, gene2, study_id) {
  fetched <- fetch_gene_pair_data(gene1, gene2, study_id)
  status  <- build_status_table(fetched, gene1, gene2)
  arm_flag <- detect_arm_level(status, gene1, gene2)
  # A gene pair is treated as flagged if EITHER the CNA-rate heuristic
  # (catches arm-level deletion) OR direct cytoband co-localization
  # (catches focal co-amplification/shared amplicons) trips.
  if (isTRUE(fetched$colocalization$same_arm)) {
    if (!arm_flag$geneA$flagged) arm_flag$geneA$flagged <- TRUE
    if (!arm_flag$geneB$flagged) arm_flag$geneB$flagged <- TRUE
    arm_flag$geneA$cytoband_flag <- TRUE
    arm_flag$geneB$cytoband_flag <- TRUE
    arm_flag$geneA$cytoband <- fetched$colocalization$gene1_pos
    arm_flag$geneB$cytoband <- fetched$colocalization$gene2_pos
  }  
  status_strict <- if (arm_flag$geneA$flagged || arm_flag$geneB$flagged) {
    status %>%
      mutate(
        geneA_altered = geneA_altered_strict,
        geneB_altered = geneB_altered_strict,
        either_altered = either_altered_strict,
        both_altered = geneA_altered_strict & geneB_altered_strict
      )
  } else NULL
  
  main <- list(
    tumor_subtype_specificity    = safe_term(term_subtype_specificity(status, fetched$clinical), "subtype"),
    biomarker_combination_signal = safe_term(term_combination_signal(status), "combination"),
    mutual_exclusivity           = safe_term(term_mutual_exclusivity(status, arm_flag, status_strict), "exclusivity"),
    tumor_suppressor_double_hit  = safe_term(term_double_hit(status, gene1, gene2, arm_flag), "double_hit"),
    cell_cycle_deregulation      = safe_term(term_cell_cycle(fetched, status), "cell_cycle"),
    angiogenesis_promotion       = safe_term(term_angiogenesis(fetched, status), "angiogenesis")
  )
  
  sensitivity <- NULL
  if (!is.null(status_strict)) {
    sensitivity <- list(
      tumor_subtype_specificity    = safe_term(term_subtype_specificity(status_strict, fetched$clinical), "subtype_strict"),
      biomarker_combination_signal = safe_term(term_combination_signal(status_strict), "combination_strict"),
      cell_cycle_deregulation      = safe_term(term_cell_cycle(fetched, status_strict), "cell_cycle_strict"),
      angiogenesis_promotion       = safe_term(term_angiogenesis(fetched, status_strict), "angiogenesis_strict")
    )
  }
  
  list(main = main, arm_flag = arm_flag, sensitivity = sensitivity)
}

# ==================== 4. SUMMARY TABLE ====================

summarize_results <- function(result) {
  pull_p <- function(r) if (!is.null(r$p_value)) r$p_value else NA_real_
  pull_txt <- function(r) {
    if (!is.null(r$interpretation)) return(r$interpretation)
    if (!is.null(r$result)) return(r$result)
    NA_character_
  }
  
  p_subtype <- pull_p(result$tumor_subtype_specificity)
  p_combination <- pull_p(result$biomarker_combination_signal)
  p_exclusivity <- pull_p(result$mutual_exclusivity)
  p_cell_cycle <- pull_p(result$cell_cycle_deregulation)
  p_angio <- pull_p(result$angiogenesis_promotion)
  
  raw_p <- c(subtype = p_subtype, pair_relationship = p_combination,
             exclusivity = p_exclusivity, cell_cycle = p_cell_cycle, angiogenesis = p_angio)
  adj_p <- p.adjust(raw_p, method = "BH")
  
  pull_verdict <- function(r) if (!is.null(r$verdict)) r$verdict else "N/A"
  
  tibble(
    term = c("Tumor subtype specificity", "Biomarker combination signal",
             "Mutually exclusive (same pathway)", "Tumor suppressor double hit (per-gene, not pair)",
             "Cell-cycle deregulation", "Angiogenesis promotion"),
    verdict = c(pull_verdict(result$tumor_subtype_specificity), pull_verdict(result$biomarker_combination_signal),
                pull_verdict(result$mutual_exclusivity), pull_verdict(result$tumor_suppressor_double_hit),
                pull_verdict(result$cell_cycle_deregulation), pull_verdict(result$angiogenesis_promotion)),
    p_value = c(p_subtype, p_combination, p_exclusivity, NA_real_, p_cell_cycle, p_angio),
    adj_p_value = c(unname(adj_p["subtype"]), unname(adj_p["pair_relationship"]),
                    unname(adj_p["exclusivity"]), NA_real_,
                    unname(adj_p["cell_cycle"]), unname(adj_p["angiogenesis"])),
    detail = c(pull_txt(result$tumor_subtype_specificity),
               pull_txt(result$biomarker_combination_signal),
               pull_txt(result$mutual_exclusivity),
               pull_txt(result$tumor_suppressor_double_hit),
               pull_txt(result$cell_cycle_deregulation),
               pull_txt(result$angiogenesis_promotion)),
    batch_adj_p_value = NA_real_
  )
}
# ==================== BATCH-LEVEL FDR (across multiple gene pairs) ====================
apply_batch_fdr <- function(results_list) {
  pair_ids <- names(results_list)
  terms <- results_list[[1]]$term
  
  for (term_name in terms) {
    raw_p <- vapply(results_list, function(tab) tab$p_value[tab$term == term_name], numeric(1))
    if (all(is.na(raw_p))) next
    batch_p <- p.adjust(raw_p, method = "BH")
    for (i in seq_along(pair_ids)) {
      results_list[[pair_ids[i]]]$batch_adj_p_value[results_list[[pair_ids[i]]]$term == term_name] <- batch_p[i]
    }
  }
  results_list
}
# ==================== 5. TERMINAL INTERACTION ====================

.stdin_con <- if (!interactive()) file("stdin", "r") else NULL

prompt_input <- function(msg, default = NULL) {
  full_msg <- if (!is.null(default)) paste0(msg, " (default: ", default, ") > ") else paste0(msg, " > ")
  if (interactive()) {
    ans <- readline(full_msg)
  } else {
    cat(full_msg)
    ans <- readLines(.stdin_con, n = 1)
  }
  ans <- sanitize_text_input(ans)
  if (identical(ans, "") && !is.null(default)) default else ans
}

sanitize_text_input <- function(x) {
  x <- trimws(x)
  x <- gsub("^[]'\":[]+", "", x)
  x <- gsub("[]'\":[]+$", "", x)
  trimws(x)
}

looks_like_valid_study_id <- function(x) {
  nzchar(x) && grepl("^[A-Za-z0-9_.\\-]+$", x)
}

select_study_interactive <- function(study_list, default_study_id = "kirc_tcga_pan_can_atlas_2018") {
  if (is.null(study_list) || nrow(study_list) == 0) {
    cat("\n(Study list unavailable -- falling back to manual entry.)\n")
    return(prompt_input("Enter TCGA study ID", default_study_id))
  }
  
  repeat {
    keyword <- prompt_input("\nSearch cancer type (e.g. 'kidney', 'breast', 'lung') -- or press Enter to see all", "")
    matches <- if (identical(keyword, "")) {
      study_list
    } else {
      study_list %>%
        filter(grepl(keyword, name, ignore.case = TRUE) |
                 grepl(keyword, cancerTypeId, ignore.case = TRUE) |
                 grepl(keyword, studyId, ignore.case = TRUE))
    }
    
    if (nrow(matches) == 0) { cat(sprintf("No studies matched '%s'.\n", keyword)); next }
    
    matches <- matches %>% arrange(name)
    n_show <- min(nrow(matches), 25)
    if (nrow(matches) > n_show) cat(sprintf("Found %d matches -- showing the first %d.\n", nrow(matches), n_show))
    
    cat("\n")
    for (i in seq_len(n_show)) cat(sprintf("  %2d. %s  [%s]\n", i, matches$name[i], matches$studyId[i]))
    cat("   0. (search again with a different keyword)\n\n")
    
    choice <- prompt_input(sprintf("Pick a study (1-%d, or 0 to search again)", n_show), "1")
    choice_num <- suppressWarnings(as.integer(choice))
    
    if (!is.na(choice_num) && choice_num == 0) next
    if (!is.na(choice_num) && choice_num >= 1 && choice_num <= n_show) return(matches$studyId[choice_num])
    cat("Not a valid choice -- try again.\n")
  }
}

print_results <- function(result, gene1, gene2, study_id) {
  cat("\n=====================================================\n")
  cat(sprintf(" RESULTS: %s + %s  |  Study: %s\n", gene1, gene2, study_id))
  cat("=====================================================\n\n")
  tab <- summarize_results(result$main)
  for (i in seq_len(nrow(tab))) {
    cat(sprintf("%d. %s  >>> %s\n", i, tab$term[i], tab$verdict[i]))
    line <- character(0)
    if (!is.na(tab$p_value[i]))     line <- c(line, sprintf("p-value: %.4g", tab$p_value[i]))
    if (!is.na(tab$adj_p_value[i])) line <- c(line, sprintf("FDR-adjusted p-value: %.4g", tab$adj_p_value[i]))
    if (length(line) > 0) cat(paste0("   ", paste(line, collapse = "   | "), "\n"))
    cat(sprintf("   %s\n\n", tab$detail[i]))
  }
  if (result$arm_flag$geneA$flagged || result$arm_flag$geneB$flagged) {
    cat("*** ARM-LEVEL CO-DELETION FLAG ***\n")
    for (g in list(result$arm_flag$geneA, result$arm_flag$geneB)) {
      if (g$flagged) {
        if (isTRUE(g$cytoband_flag) && !is.null(g$event_type)) {
          verb <- if (g$event_type == "amplification") "gain/amp" else "del"
          noun <- if (g$event_type == "amplification") "high-level amplification" else "deep/focal loss"
          cat(sprintf("  %s (%s): %.1f%% broad %s, only %.1f%% of those are focal, AND physically co-localized with the other gene -- signature of a shared CNA/amplicon event.\n",
                      g$gene, g$cytoband, 100 * g$broad_rate, verb, 100 * g$focal_rate))
        } else if (isTRUE(g$cytoband_flag)) {
          cat(sprintf("  %s (%s): physically co-localized with the other gene on the same chromosome arm -- co-occurrence may reflect a shared CNA event rather than independent biology.\n",
                      g$gene, g$cytoband))
        } else {
          verb <- if (g$event_type == "amplification") "gain/amp" else "del"
          noun <- if (g$event_type == "amplification") "high-level amplification" else "deep/focal loss"
          cat(sprintf("  %s: %.1f%% broad %s, only %.1f%% of those are focal -- signature of an arm-level/amplicon event, not %s.\n",
                      g$gene, 100 * g$broad_rate, verb, 100 * g$focal_rate, noun))
        }
      }
    }
    cat("  Terms 1, 2, 5, 6 use a broad 'altered' definition (mut/del/amp) that this event\n")
    cat("  inflates. Re-run below with a STRICT definition (mut/deep-del/amp only,\n")
    cat("  shallow/arm-level deletions excluded) as a sensitivity check.\n\n")
    
    strict_tab <- summarize_results(result$sensitivity)
    cat("  --- STRICT sensitivity results ---\n")
    for (i in seq_len(nrow(strict_tab))) {
      if (is.na(strict_tab$p_value[i]) && is.na(strict_tab$adj_p_value[i]) &&
          strict_tab$term[i] %in% c("Mutually exclusive (same pathway)",
                                    "Tumor suppressor double hit (per-gene, not pair)")) next
      cat(sprintf("  %s: p=%.4g\n", strict_tab$term[i],
                  ifelse(is.na(strict_tab$p_value[i]), NA, strict_tab$p_value[i])))
    }
    cat("  Compare these p-values to the main results above:\n")
    cat("    - If a term stops being significant under strict, its 'main' result was\n")
    cat("      likely driven by the arm-level event rather than gene-specific alteration.\n")
    cat("    - If a term BECOMES significant under strict, the arm-level event was adding\n")
    cat("      noise that masked a SIGNIFICANT in the broad 'altered' definition -- the\n")
    cat("      strict result is likely the more trustworthy one.\n\n")
  }
  cat("-----------------------------------------------------\n")
  cat("TEST TYPES & ADJUSTMENT LOGIC:\n")
  cat("  - Terms 2 & 3 directly compare gene A vs gene B (true pair relationship).\n")
  cat("  - Terms 1, 5 & 6 test 'either gene altered' vs a subtype/pathway score.\n")
  cat("  - Term 4 is calculated separately per gene and is NOT a pair comparison.\n")
  cat("  - To prevent suppression artifacts, the headline p-value is ALWAYS the\n")
  cat("    most conservative (largest) of the raw and burden-adjusted p-values.\n")
  cat("  - Term 1 uses LRT burden adjustment to isolate subtype from mutation count.\n")
  cat("  - Term 3 uses mutation-only exclusivity to handle co-localized genes.\n")
  cat("  - Term 4 flags arm-level co-deletion artifacts for double-hit rates.\n")
  cat("-----------------------------------------------------\n\n")
}

run_terminal_app <- function() {
  cat("\n===== TCGA Gene-Pair Multi-Term Analysis Tool (v5) =====\n")
  cat("(Data source: cBioPortal REST API -- no Bioconductor required)\n\n")
  
  cat("Checking connectivity to cBioPortal...\n")
  if (!check_network_connectivity()) {
    cat("\nContinuing anyway in case this was transient.\n\n")
  } else {
    cat("OK.\n\n")
  }
  
  cat("Loading list of available studies...\n")
  study_list <- tryCatch(get_study_list(), error = function(e) {
    cat(sprintf("!! Could not load study list (%s).\n", conditionMessage(e)))
    NULL
  })
  if (!is.null(study_list)) cat(sprintf("Loaded %d studies.\n", nrow(study_list)))
  
  repeat {
    gene1    <- toupper(prompt_input("Enter Gene 1 (Hugo symbol)", "VHL"))
    gene2    <- toupper(prompt_input("Enter Gene 2 (Hugo symbol)", "VEGFA"))
    study_id <- select_study_interactive(study_list)
    
    if (!looks_like_valid_study_id(study_id)) {
      cat(sprintf("\n!! '%s' doesn't look like a valid study ID.\n", study_id))
      next
    }
    
    result <- tryCatch(
      analyze_gene_pair(gene1, gene2, study_id),
      error = function(e) { cat("\n!! ERROR:", conditionMessage(e), "\n\n"); NULL }
    )
    
    if (!is.null(result)) print_results(result, gene1, gene2, study_id)
    
    again <- tolower(prompt_input("Analyze another gene pair? (y/n)", "y"))
    if (!identical(again, "y")) break
  }
  
  if (!is.null(.stdin_con)) close(.stdin_con)
  cat("\nDone.\n")
}

# Auto-run when this file is executed
#run_terminal_app()


