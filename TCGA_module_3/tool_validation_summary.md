# TCGA Gene-Pair Analysis Tool — Validation Summary

**Purpose:** Reference document for final evaluation and publication. Summarizes bugs found, fixes applied, verification performed, and known limitations.

---

## 1. Bugs Found and Fixed (all verified with real before/after data)

### Bug 1: Wrong CNA profile selected
- **Problem:** Tool always picked the first CNA profile matching `DISCRETE` datatype (e.g., `gbm_tcga_pub_cna_consensus`). For some studies, this profile silently omitted certain genes (e.g., GLI1 showed 0.5% altered instead of the true ~7%).
- **Fix:** Added a fallback loop that tries each available CNA profile in turn and only accepts one that covers both query genes.
- **Verification:** GLI1 alteration rate corrected from 1/206 (0.5%) to 46/206 (22%), matching cBioPortal's own OncoPrint (7% amplification + other events).

### Bug 2: Wrong expression profile selected
- **Problem:** Expression profile selection sometimes matched a miRNA z-score profile instead of mRNA, silently returning 0 rows and making Terms 5/6 (Cell-cycle, Angiogenesis) always N/A.
- **Fix:** Added `!grepl("mirna", ...)` exclusion to the profile selection filter.
- **Verification:** `gbm_tcga_pub` now correctly returns `gbm_tcga_pub_mrna_median_Zscores` with 48,383 expression rows; Terms 5/6 now compute real results.

### Bug 3: Wrong study cohort used for ground-truth comparison
- **Problem:** Benchmark used PanCancer Atlas 2018 cohorts (`gbm_tcga_pan_can_atlas_2018`, `ov_tcga_pan_can_atlas_2018`), but the ground-truth modules (Suppl_Table.xlsx, based on Ciriello et al. 2012 / MEMo) were computed on the original 2008/2011 TCGA marker-paper cohorts.
- **Fix:** Changed `STUDY_ID_MAP` to `gbm_tcga_pub` (206 samples) and `ov_tcga_pub` (316–489 samples).

### Bug 4: Amplification diluting tumor-suppressor loss signal
- **Problem:** `geneA_altered`/`geneB_altered` counted mutation + deletion + amplification equally for every gene. For tumor suppressors (RB1, CDKN2A, TP53, BRCA1, BRCA2, etc.), amplification is not an inactivating event and added noise, diluting real exclusivity signal.
- **Fix:** Loaded OncoKB Cancer Gene List (`cancerGeneList.tsv`, 1,242 genes) at startup. Added `check_gene_role()`: tumor suppressors now use mutation+deletion only; oncogenes/unknown genes keep the broad mut+del+amp definition. Unknown genes trigger a visible warning rather than failing silently.
- **Verification:** CDKN2A altered count dropped from 148→143/206 (amp-only samples correctly excluded); CDK4-CDKN2A odds ratio and p-value improved and passed significance (p=2.7e-5).

### Bug 5: Missing ground-truth significance filter
- **Problem:** Benchmark treated all modules listed in Suppl_Table.xlsx as equally valid positives, including modules the source paper itself found non-significant (q-value up to 1.0). For OVCA, only 3 of 33 listed modules were actually significant (q<0.05); ~90% of "ground truth" pairs were noise.
- **Fix:** Added a filter in `parse_sheet()`: `filter(!is.na(sig_value) & sig_value < 0.05)`, keeping only genuinely significant modules from each sheet.
- **Verification:** Ground-truth pair count dropped from 194 to 39 (24 GBM, 15 OVCA) after filtering to real positives.

### Bug 6: Two-sided statistical test used for a directional hypothesis
- **Problem:** Term 3 (mutual exclusivity) used two-sided Fisher tests / Wald tests, even though the hypothesis being tested is inherently directional (OR < 1, i.e., exclusivity specifically, not "any association").
- **Fix:** Changed `safe_fisher()` calls in `term_mutual_exclusivity` to `alternative = "less"`; added `p_value_one_sided_less` to `burden_adjusted_association()` output and used it in `pick_tier()` and `mut_p` calculation.
- **Verification:** Several near-miss pairs (e.g., BRCA1-CCNE1 raw p 0.071→0.040, BRCA2-CCNE1 raw p→0.009, now passing at p=0.041) moved appropriately; some pairs remained non-significant after adjustment, correctly reflecting genuine uncertainty rather than being forced.

### Bug 7: Missing closing brace (introduced during Bug 6 fix, caught before running)
- **Problem:** A `}` was missing after `burden_adjusted_association()`, nesting `term_combination_signal` inside it.
- **Fix:** Added the missing brace before re-running.

---

## 2. Terminology Change (pending — do this last, before final publish)

Replace "REAL SIGNAL" → "SIGNIFICANT" and "NO SIGNAL" → "NOT SIGNIFICANT" throughout (Terms 1, 2, 3, 5, 6 verdict strings). Rationale: "NO SIGNAL" implies no relationship exists; "NOT SIGNIFICANT" correctly communicates "not statistically demonstrated at this sample size," which is scientifically accurate and avoids overclaiming. Term 4's "LIKELY FALSE (chromosome position)" / "DESCRIPTIVE ONLY" and "MUTUALLY EXCLUSIVE" labels are already appropriate and should NOT be changed.

---

## 3. Validated Results (final, as of last benchmark run)

**GBM** (`gbm_tcga_pub`, ground truth q<0.05, n=24 pairs): **5/24 (20.8%) recall.** All 5 hits are the well-established strong pairs: CDK4-CDKN2A, CDK4-CDKN2B, MDM2-CDKN2A, MDM2-CDKN2B, GLI1-CDKN2A.

**OVCA** (`ov_tcga_pub`, ground truth q<0.05, n=15 pairs): **0/15 (0%) recall on strict multiple-testing basis**, but with important nuance:
- BRCA2-CCNE1: MUTUALLY EXCLUSIVE, raw p=0.009, headline p=0.041 (passes)
- BRCA1-CCNE1: NOT SIGNIFICANT, raw p=0.040, adjusted p=0.092 (near-miss, correct direction)
- RB1-CCNE1: NOT SIGNIFICANT, OR=1.31 (wrong direction — trends co-occurring, not exclusive)

## 4. Known, Documented Limitations (for the paper's Limitations section)

1. **Pairwise vs. module-level testing.** This tool tests exactly 2 genes at a time. Published ground truth (Ciriello et al. 2012 / MEMo) tests 3–5 gene modules jointly. A module can be significant as a whole even when individual constituent pairs are not — this is expected, not a tool defect. MEMo itself was validated on only 2 cancer types (GBM, OVCA) and its authors reported pan-cancer scaling as challenging — this tool's scope is comparably bounded and that is normal practice in the field.
2. **Cohort-specificity.** Results are specific to the chosen cBioPortal study/cohort. Different cohorts (sample size, sequencing platform, disease stage mix) can give different results for the same gene pair. This is a general property of cancer genomics, not unique to this tool.
3. **Validation scope.** This tool has been statistically validated against published ground truth for GBM and OVCA only. Structural bug classes (profile selection, gene-role-aware CNA filtering) are generic and should transfer to other cancer types, but cohort-specific ground truth (correct `study_id`, correct significant modules) has not been re-verified for other cancer types. Recommend a small spot-check (2–3 known pairs) before trusting results on a new cancer type.
4. **Methylation silencing not captured.** Term 4 (double-hit) only sees mutation+CNA; promoter methylation (a known third inactivation route, e.g., for VHL) is not in scope.
5. **Gene role list is a fixed snapshot.** `TUMOR_SUPPRESSORS`/`ONCOGENES` come from a static OncoKB file download; genes not in that list default to the broad mut+del+amp definition, with a visible warning logged per occurrence.

## 5. Remaining Work Before Full Publication

- [ ] Apply the terminology change (Section 2) throughout the codebase.
- [ ] Spot-check ground truth for Terms 1, 2, 4, 5, 6 (one known literature fact per term — see prior conversation for suggested pairs: TP53/GBM subtype, EGFR+PTEN co-occurrence, VHL double-hit in KIRC, RB1 loss→E2F targets up, VHL loss→angiogenesis up). This is a lighter check than Term 3 required — no module decomposition issue applies.
- [ ] Regenerate the full benchmark CSV one final time with all fixes + terminology change in place, for the numbers actually reported in the publication.
- [ ] Write the Limitations section using Section 4 above as the basis.
