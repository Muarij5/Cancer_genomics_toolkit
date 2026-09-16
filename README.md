# Cancer Genomics Toolkit

A three-module suite for cancer genomics analysis: synthetic lethality prediction,
drug synergy scoring, and TCGA-based gene-pair mutual exclusivity analysis.
Developed independently; validation and benchmarking documented per module.

## Modules

### Module 1 — GeneComb Analyzer (`Module_1/`)
Predicts candidate synthetic lethal / functional redundancy gene pairs using
DepMap dependency data — combining differential essentiality, mutual exclusivity,
and co-expression signals, with an automatic arm-level/co-amplicon artifact
detector (checks focal-vs-broad CNA rates and genomic co-localisation ≤10 Mb)
to flag likely chromosomal-proximity artifacts rather than true biological signal.
Benchmarked against the Dede/Hart and Parrish (2021) synthetic lethality datasets.

**Note:** requires large DepMap release files (CRISPRGeneEffect, Omics matrices,
mart_export) not included in this repo due to size — see Setup below.

### Module 2 — CancerComb Analyzer (`Module_2/`)
A drug synergy scoring pipeline implementing Bliss and Loewe synergy models on dose-response data, with `(0,0)` control normalization and transparent exact-vs-approximate Loewe reporting. Benchmarked for conformance against SynergyFinder's bundled reference dataset. Callable via `run_pipeline(csv_path, unit, response_type)`.

### Module 3 — TCGA Gene-Pair Multi-Term Analysis (`TCGA_module_3/`)
Queries the live cBioPortal REST API to test pairwise mutual exclusivity and
co-occurrence across TCGA cohorts, with Benjamini-Hochberg correction applied
both within-pair and across the batch. Includes literature spot-checks against
published gene-pair claims (MEMo and related sources).

## Setup

Each module is a set of R scripts; no shared package structure yet (see Status).

- **Module_1**: requires DepMap release CSVs (CRISPRGeneEffect, Omics*, Model,
  PortalCompounds) placed in `Module_1/`, plus `mart_export.txt` from Ensembl —
  excluded from this repo via `.gitignore` due to size (~1.5 GB total).
- **Module_2**: run `final_module_2.R` — call `run_pipeline(csv_path, unit, response_type)` in R/RStudio.
- **Module_3**: run `Finalized_module3_CORRECTED.R` — queries cBioPortal live,
  no local data setup required beyond the included `cancerGeneList.tsv`
  (OncoKB static snapshot).

## Validation

Each module includes its own benchmarking scripts and result files
(`benchmark_comparison_results*.csv`, `external_benchmark_results*.csv`,
literature spot-check scripts) documenting a find-bug → fix → quantify
validation trail. Full validation summaries are in the `.docx` files per module.

## Status

This is a working research prototype, not production software:
- No shared package structure or automated tests yet
- No CI/CD or dependency pinning (`renv`) yet
- Statistical calibration (random-pair null baseline, precision reporting
  alongside recall) is an active area of ongoing refinement
- Module 1's large data dependency currently limits it to local/offline use

Developed as an independent project; feedback and suggestions welcome.

## License

For academic and research reference purposes.
