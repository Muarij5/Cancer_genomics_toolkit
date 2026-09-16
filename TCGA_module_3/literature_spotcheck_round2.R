# ============================================================
# literature_spotcheck_round2.R
# Adds a SECOND literature-grounded pair per term (n=2 per term
# total, combined with round 1). Uses analyze_single_gene() where
# applicable (new in the corrected file) instead of GAPDH-pairing.
# Run AFTER Finalized_module3_CORRECTED.R is sourced.
# ============================================================

# ------------------------------------------------------------
# Check 6 -- Term 1 (Subtype): IDH1 mutation defines Proneural GBM
# Literature: Verhaak et al. 2010, Cancer Cell -- IDH1 mutation is the
# single strongest defining feature of the Proneural GBM subtype.
# ------------------------------------------------------------
cat("\n\n################ CHECK 6: IDH1 -> Proneural subtype (GBM) ################\n")
r6 <- analyze_single_gene("IDH1", "gbm_tcga_pub")
cat("Subtype verdict:", r6$subtype$verdict %||% "N/A", "\n")
cat(r6$subtype$interpretation %||% r6$subtype$result %||% "", "\n")

# ------------------------------------------------------------
# Check 7 -- Term 2 (Combination): KRAS + STK11 co-occurrence in LUAD
# Literature: well-documented co-occurring alterations defining the
# STK11/KEAP1-mutant KRAS subtype of lung adenocarcinoma (Skoulidis
# et al. 2018, Cancer Discovery).
# ------------------------------------------------------------
cat("\n\n################ CHECK 7: KRAS+STK11 co-occurrence (LUAD) ################\n")
r7 <- analyze_gene_pair("KRAS", "STK11", "luad_tcga_pan_can_atlas_2018")
print(summarize_results(r7$main) %>% filter(term == "Biomarker combination signal"))

# ------------------------------------------------------------
# Check 8 -- Term 4 (Double-hit): TP53 near-universal biallelic
# inactivation in high-grade serous ovarian cancer.
# Literature: TCGA Ovarian 2011, Nature -- TP53 mutated/inactivated
# in ~96% of HGSOC, classically via mutation + LOH.
# ------------------------------------------------------------
cat("\n\n################ CHECK 8: TP53 double-hit (OVCA) ################\n")
r8 <- analyze_single_gene("TP53", "ov_tcga_pub")
cat("Double-hit:", r8$double_hit$interpretation, "\n")

# ------------------------------------------------------------
# Check 9 -- Term 5 (Cell-cycle): CDKN2A loss activates CDK4/6-RB-E2F axis
# Literature: CDKN2A (p16) loss de-represses CDK4/6, which
# phosphorylates RB1, releasing E2F -- classic, textbook cell-cycle
# deregulation mechanism, independent of the RB1-direct check in round 1.
# ------------------------------------------------------------
cat("\n\n################ CHECK 9: CDKN2A -> cell-cycle deregulation (GBM) ################\n")
r9 <- analyze_single_gene("CDKN2A", "gbm_tcga_pub")
cat("Cell-cycle:", r9$cell_cycle$verdict %||% "N/A", "\n")
cat(r9$cell_cycle$interpretation %||% r9$cell_cycle$result %||% "", "\n")

# ------------------------------------------------------------
# Check 10 -- Term 6 (Angiogenesis): EGFR amplification drives
# angiogenic signaling in GBM (independent of the VHL/HIF-axis check
# in round 1 -- different mechanism, EGFR->VEGF induction).
# Literature: well-established in GBM biology (e.g. Guo et al. 2003,
# Am J Pathol; standard GBM/EGFR/VEGF literature).
# ------------------------------------------------------------
cat("\n\n################ CHECK 10: EGFR -> angiogenesis (GBM) ################\n")
r10 <- analyze_single_gene("EGFR", "gbm_tcga_pub")
cat("Angiogenesis:", r10$angiogenesis$verdict %||% "N/A", "\n")
cat(r10$angiogenesis$interpretation %||% r10$angiogenesis$result %||% "", "\n")

cat("\n\n================ ALL 5 ADDITIONAL CALLS DONE (round 2) ================\n")
cat("Combined with round 1, every term now has 2 independent literature spot-checks.\n")


length(TUMOR_SUPPRESSORS)
length(ONCOGENES)
"TP53" %in% TUMOR_SUPPRESSORS