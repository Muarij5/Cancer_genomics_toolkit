# ============================================================
# literature_spotcheck_round3.R
# Third literature-grounded pair per term (n=3 per term total,
# combined with rounds 1 and 2). Run AFTER
# Finalized_module3_CORRECTED.R is sourced.
# ============================================================

# ------------------------------------------------------------
# Check 11 -- Term 1 (Subtype): TP53 mutation enriched in Basal-like
# breast cancer subtype.
# Literature: TCGA Breast 2012, Nature -- TP53 mutation rate is
# dramatically higher in Basal-like tumors (~80%) vs Luminal A (~12%).
# ------------------------------------------------------------
cat("\n\n################ CHECK 11: TP53 -> Basal-like subtype (Breast) ################\n")
r11 <- analyze_single_gene("TP53", "brca_tcga_pub")
cat("Subtype verdict:", r11$subtype$verdict %||% "N/A", "\n")
cat(r11$subtype$interpretation %||% r11$subtype$result %||% "", "\n")

# ------------------------------------------------------------
# Check 12 -- Term 2 (Combination): ARID1A + PIK3CA co-occurrence in
# endometrioid endometrial cancer (UCEC).
# Literature: TCGA UCEC 2013, Nature -- ARID1A and PIK3CA alterations
# are well-documented as frequently co-occurring in this cohort.
# ------------------------------------------------------------
cat("\n\n################ CHECK 12: ARID1A+PIK3CA co-occurrence (UCEC) ################\n")
r12 <- analyze_gene_pair("ARID1A", "PIK3CA", "ucec_tcga_pan_can_atlas_2018")
print(summarize_results(r12$main) %>% filter(term == "Biomarker combination signal"))

# ------------------------------------------------------------
# Check 13 -- Term 4 (Double-hit): PTEN near-universal biallelic loss
# in endometrioid endometrial cancer.
# Literature: TCGA UCEC 2013 -- PTEN is inactivated (mutation + LOH)
# in the large majority of endometrioid-type tumors, a classic
# double-hit tumor suppressor example independent of VHL/TP53 above.
# ------------------------------------------------------------
cat("\n\n################ CHECK 13: PTEN double-hit (UCEC) ################\n")
r13 <- analyze_single_gene("PTEN", "ucec_tcga_pan_can_atlas_2018")
cat("Double-hit:", r13$double_hit$interpretation, "\n")

# ------------------------------------------------------------
# Check 14 -- Term 5 (Cell-cycle): MYC amplification directly
# transactivates E2F1 and E2F target genes.
# Literature: classic MYC-E2F1 transcriptional axis (e.g. Leone et al.
# 1997, Nature); independent mechanism from the RB1/CDKN2A checks above.
# ------------------------------------------------------------
cat("\n\n################ CHECK 14: MYC -> cell-cycle/E2F targets (Breast) ################\n")
r14 <- analyze_single_gene("MYC", "brca_tcga_pub")
cat("Cell-cycle:", r14$cell_cycle$verdict %||% "N/A", "\n")
cat(r14$cell_cycle$interpretation %||% r14$cell_cycle$result %||% "", "\n")

# ------------------------------------------------------------
# Check 15 -- Term 6 (Angiogenesis): PTEN loss drives angiogenesis via
# the PI3K/AKT/HIF1A/VEGF axis.
# Literature: well-established PI3K/AKT-driven HIF1A/VEGF induction
# mechanism in GBM, distinct from the VHL (HIF-stabilization) and
# EGFR (direct VEGF induction) mechanisms already checked.
# ------------------------------------------------------------
cat("\n\n################ CHECK 15: PTEN -> angiogenesis (GBM) ################\n")
r15 <- analyze_single_gene("PTEN", "gbm_tcga_pub")
cat("Angiogenesis:", r15$angiogenesis$verdict %||% "N/A", "\n")
cat(r15$angiogenesis$interpretation %||% r15$angiogenesis$result %||% "", "\n")

cat("\n\n================ ALL 5 ADDITIONAL CALLS DONE (round 3) ================\n")
cat("Combined with rounds 1 and 2, every term now has 3 independent literature spot-checks.\n")














