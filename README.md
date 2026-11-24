# OSCC Dendritic-Cell Atlas (`oscc-dc-atlas`)

Analysis code, configuration tables, and source data accompanying the study:

> **Maturation and regulatory expression in migratory dendritic cells of oral squamous cell carcinoma across patient-paired cohorts**

## Overview

We integrate 10 single-cell RNA-seq (scRNA-seq) libraries from 8 patients with oral squamous cell carcinoma (OSCC; 85,164 quality-controlled cells, including 1,372 curated dendritic cells across six subtypes: cDC1, cDC2, Langerhans-like [LC-like], migratory/mregDC [migDC], pDC, and proliferating DC) to dissect migratory DC transcriptional programs using donor-paired pseudobulk negative-binomial models (`edgeR` quasi-likelihood), self-contained rotation gene-set tests (`limma::fry`), competitive pathway enrichment (`CAMERA`), and patient-level tumour microenvironment (TME) covariation. Independent replication is performed across 10 primary oral-cavity tumours in `GSE181919` (`2,902` myeloid antigen-presenting cells) after recovering cDC2 cells from the broad myeloid compartment.

![External validation and cross-cohort meta-analysis](results/figures/Fig5_external_validation.png)

## Key findings

- **Donor-paired `migDC` vs `cDC2` differential expression**: Across 5 primary tumours meeting the 10-cell-per-subtype threshold (`4,151` tested genes), `933` genes reach `FDR < 0.05`, including *CD83* ($\log_2\text{FC} = +0.70$), *PDCD1LG2* ($+2.73$), and *CD200* ($+6.46$).
- **Coupled maturation, feedback, and regulatory programs**: Held-out `Regulatory expression` ($+2.74\ \log_2\text{CPM}$), `Maturation` ($+1.01$), and `NF-κB feedback` ($+1.06$) are consistently elevated in `migDC` across all 5 discovery patients (paired *t* `P < 0.05` and directional `limma::fry` `FDR < 0.05`), whereas broad `HALLMARK_INFLAMMATORY_RESPONSE` enrichment favours `cDC2` (`CAMERA` `q = 0.00671`).
- **Mucosal `LC-like` specificity and *IDO1* distribution**: *IDO1* is elevated in both `migDC` ($+3.20\ \log_2\text{FC}$) and `cDC1` ($+2.15$) relative to `cDC2`, whereas `LC-like` DCs resemble `cDC2` ($+0.24$, `FDR = 0.752`).
- **Patient-level TME covariation**: Across 8 primary tumours, retained-DC *CCL17* correlates with Treg *CCR4* expression (Spearman $\rho = 0.93$, `BH FDR = 0.0362`), alongside positive `CCL19`–`CCR7` covariation ($\rho = 0.88$, `FDR = 0.0809`).
- **Cross-cohort replication (`GSE181919`)**: Recovering `243/245` (`99.2%`) gated `cDC2` cells from the broad `Macrophages` cluster yields 6–7 paired validation donors ($\rho = 0.671\text{--}0.679$; `93.4%–94.8%` hit concordance; $\rho = 0.751$ and `96.8%` in the 4-patient CD163-negative hybrid model). Across combined discovery and validation cohorts, `Maturation`, `NF-κB feedback`, and `Regulatory expression` are concordantly upregulated in `11/11` (10-cell threshold, exact sign-flip $P = 9.77 \times 10^{-4}$) and `12/12` (5-cell threshold, $P = 4.88 \times 10^{-4}$) patients.

## Datasets

| GEO Accession | Cohort Role | Samples Included | Reference |
|---|---|---|---|
| [GSE172577](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE172577) | Discovery atlas | 6 primary oral-tongue tumours | Peng *et al.*, *Oral Oncol* 2021 (PMID: 34044317) |
| [GSE193766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193766) | Discovery atlas | 2 primary OSCC tumours + paired adjacent mucosa | Chen *et al.*, *Sci Adv* 2022 (PMID: 35613277) |
| [GSE181919](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181919) | Independent validation | 10 primary oral-cavity tumours (`CA`, `OC`) | Choi *et al.*, *Nat Commun* 2023 (PMID: 36828832) |

## Repository structure

```text
oscc-dc-atlas/
├── config/                  # Sample manifest, cell-type maps, gene programs, and DC gates
├── data/
│   ├── genesets/            # Cached MSigDB Hallmark, KEGG, and C3 TFT GMT files
│   ├── checksums.csv        # MD5 checksums for discovery input objects
│   └── README.md            # Input data download and setup instructions
├── docs/
│   ├── data_provenance.md   # GEO cohort details and sample accession mapping
│   ├── methods.md           # Full computational methods supplement
│   └── sessionInfo.txt      # R / Bioconductor package versions
├── R/                       # Numbered modular analysis pipeline (00_setup.R -> 14_manuscript_figures.R)
│   ├── pseudobulk.R         # Shared donor-paired aggregation, edgeR, and fry helpers
│   └── verify_primary.R     # Standalone primary-model reproduction from tracked source counts
├── results/
│   ├── figures/             # Main figures (Fig1–Fig5) and supplementary figures (PNG + PDF)
│   └── tables/              # Source data tables, DE results, program scores, and TME correlations
├── Makefile
├── run_all.sh
└── LICENSE
```

## Quick start

### 1. Reproduce the primary donor-paired analysis from tracked source counts

No raw downloads or Seurat objects are required to reproduce the primary `migDC` vs `cDC2` `edgeR` estimates using the tracked aggregate matrix (`results/tables/11_primary_pseudobulk_counts.tsv.gz`):

```bash
Rscript -e 'install.packages("BiocManager"); BiocManager::install(c("edgeR", "statmod"))'
Rscript R/verify_primary.R
```

### 2. Run the full pipeline

After placing the GEO inputs under `data/raw/` (see [`data/README.md`](data/README.md)):

```bash
./run_all.sh          # Discovery atlas, donor-paired models, and TME covariation (01–11, 13)
./run_all.sh 12 14    # Independent GSE181919 validation (12) and main figures (14)
```

Or selectively via `make`:

```bash
make all              # Run steps 01–11 and 13
make validate         # Run GSE181919 validation (step 12)
make figures          # Regenerate main composite figures (step 14)
```

## License

Code is released under the [MIT License](LICENSE). Generated figures and tables (`results/figures/`, `results/tables/`) are available under CC-BY-4.0. Underlying public sequencing datasets remain subject to the original GEO data usage terms.
