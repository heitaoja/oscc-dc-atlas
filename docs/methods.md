# Computational Methods

All workflow parameters are centralized in [`R/00_setup.R`](../R/00_setup.R) (`PARAMS`) with fixed random seed `20240601`.

## 1. Discovery single-cell atlas and annotation

1. **Input processing and gene harmonization**: Integer RNA-layer UMI counts from six primary oral-tongue tumours (`GSE172577`) and two matched primary tumour–adjacent mucosa pairs (`GSE193766`; 10 libraries across 8 donors) are mapped from Ensembl IDs to HGNC symbols using `org.Hs.eg.db` (`multiVals = "first"`). Counts from duplicated gene symbols are summed via sparse matrix aggregation, and the intersection of shared symbols across series (`20,814` genes) is retained.
2. **Quality control and doublet removal**: Per-library filtering removes cells with `< 200` or `> 6,000` detected genes, `< 500` UMIs, `> 15%` mitochondrial reads, `> 3%` haemoglobin reads, or `log10(nFeature_RNA) > median + 3 × MAD` within each library, followed by per-library doublet removal using `scDblFinder` (`85,164` quality-controlled cells retained).
3. **Normalization, integration, and clustering**: Counts are log-normalized (`scale.factor = 1e4`), and the top `2,000` variable genes (`vst`) are scaled and embedded via PCA (`50` PCs). Batch correction across the 10 libraries is performed with `Harmony` (`theta = 2`, first `30` PCs), followed by Leiden clustering (`resolution = 0.5`) and UMAP (`dims = 1:30`).
4. **Major lineage and DC subtype annotation**: Major cell types (`24` lineages) are assigned by integrating cluster markers (`FindAllMarkers`, Wilcoxon test) and `SingleR` voting across three reference atlases (`BlueprintEncodeData`, `MonacoImmuneData`, `HumanPrimaryCellAtlasData`) in [`config/annotation_major.csv`](../config/annotation_major.csv). Cells annotated as `DC` or `pDC` (`1,872` cells) are re-clustered (`1,500` variable features, `30` PCs, `Harmony` by library, Leiden `resolution = 0.4`), excluding two low-complexity/ambient epithelial clusters (`500` cells) to retain `1,372` DCs across six subtypes (`cDC1`, `cDC2`, `LC-like`, `migDC`, `pDC`, `Proliferating DC`; [`config/annotation_dc.csv`](../config/annotation_dc.csv)).
5. **Composition summaries**: Candidate and retained DC fractions use total quality-controlled cells per library as the denominator; subtype fractions use retained DCs per library. Matched tumour–adjacent pairs (`P7`, `P8`) are summarized descriptively without inferential testing.

## 2. Donor-paired pseudobulk differential expression

To avoid cell-level pseudoreplication, raw integer RNA counts are summed within each primary tumour library and DC subtype (`R/11_donor_pseudobulk.R`, `R/pseudobulk.R`). Eligible donors require at least `min_cells` cells in both compared subtypes within the same primary tumour, and each model requires at least `4` eligible donors.

- **Primary and lineage contrasts**: `migDC` vs `cDC2` (reference), `migDC` vs `LC-like`, and `LC-like` vs `cDC2` require `≥ 10` cells per state (`5` eligible donors); `migDC` vs `cDC1` and `cDC1` vs `cDC2` require `≥ 5` cells per state (`5` eligible donors).
- **Filtering and model fitting**: Expression filtering uses `edgeR::filterByExpr` grouped by `dc_subtype` with `min.count = 3` and `min.total.count = 10` tailored to sparse DC pseudobulk library sizes. Library sizes are recalculated and normalized via TMM, followed by robust dispersion estimation and quasi-likelihood negative-binomial regression (`~ patient + dc_subtype`) with `edgeR::glmQLFTest`. Benjamini–Hochberg (BH) FDR is controlled across all tested genes per contrast.
- **Sensitivity analyses**:
  1. *Leave-one-donor-out*: Re-fit after omitting each primary donor and compute shared-gene Spearman rank correlation ($\rho$).
  2. *Cell support thresholds*: Evaluate `3`, `5`, `15`, and `20` cells per state.
  3. *Normalization & filtering*: Re-fit with `TMMwsp` and under the standard bulk filter (`strict_filter_min10`: `min.count = 10`, `min.total.count = 15`).
  4. *Equal-cell subsampling*: Downsample the larger subtype to match the smaller subtype within each tumour (`40` repeats) and recalculate TMM-normalized paired scores.
  5. *Objective marker-gate reclassification*: Reassign cells within the `cDC1`/`cDC2`/`migDC` pool using mutually exclusive two-of-three detection gates (`cDC1`: *CLEC9A/XCR1/CADM1*; `cDC2`: *CD1C/CLEC10A/FCER1A*; `migDC`: *LAMP3/CCR7/FSCN1*; [`config/external_dc_gates.csv`](../config/external_dc_gates.csv)) at 5- and 10-cell minima, comparing non-assignment fold-change ranks.

## 3. Transcriptional program and pathway analyses

- **Held-out expression programs**: Four literature-curated gene panels (`Maturation`, `NF-kB feedback`, `Regulatory expression`, `MHC class II`; [`config/dc_programs.csv`](../config/dc_programs.csv)) exclude all explicit annotation and gate markers (`CD274` and `IDO1` are held out of the main discovery panels, and additionally evaluated in an expanded `Regulatory expression (incl. CD274/IDO1)` panel under the 9-gene marker-gate classification where neither gene participates in gating). Programs require `≥ 60%` target gene coverage post-filtering. Donor scores are computed as the mean within-patient `migDC` minus `cDC2` difference on TMM $\log_2\text{CPM}$ (`prior.count = 1`) and tested via two-sided exact sign-flip tests, one-sample paired *t* tests, and self-contained rotation gene-set tests (`limma::fry` on the filtered $\log_2\text{CPM}$ matrix with `~ patient + dc_subtype`).
- **Competitive pathway enrichment (`CAMERA`)**: `limma::camera` tests MSigDB Hallmark sets (`data/genesets/h.all.symbols.gmt`, `≥ 10` measurable members, default inter-gene correlation `0.01`) on the primary paired design after removing all explicit annotation and gate genes while preserving primary TMM effective library sizes.

## 4. Tumour microenvironment covariation

Across all quality-controlled atlas cells (`R/13_tme_covariation.R`), library-level proportions are computed relative to total immune cells (`T cells`, `Treg`, `NK cells`, `B cells`, `Plasma cells`, `Macrophages`, `DC`, `pDC`, `Mast cells`, `Neutrophils`), total `T/NK` cells (`T cells + Treg + NK cells`), and retained DCs. Within primary tumours, integer counts for retained DCs, `migDC` (`≥ 3` cells), `Treg`, `T cells`, and `Macrophages` are aggregated per library and TMM-normalized on each population's full transcriptome. Predefined ligand–receptor, composition, and paired-program axes are tested across primary tumours using two-sided Spearman correlation (`cor.test(..., exact = FALSE)`) with BH FDR correction across all testable axes.

## 5. Independent validation in GSE181919 and cross-cohort meta-analysis

Primary oral-cavity tumours (`tissue.type == "CA"`, `subsite == "OC"`, 10 patients, `10,397` cells) from `GSE181919` are extracted via streaming column selection (`R/12_external_validation.R`). Because coarse clustering in `GSE181919` assigns `243/245` (`99.2%`) `CD1C/CLEC10A/FCER1A`-gated `cDC2` cells to `Macrophages` (`2,214` cells) rather than `Dendritic.cells` (`688` cells), the full myeloid antigen-presenting-cell pool (`2,902` cells) is extracted and evaluated alongside `CD163`-negative (`CD163` UMI `== 0`) and hybrid (`Dendritic.cells` `migDC` vs `CD163`-negative `cDC2`) sensitivity models. Paired `edgeR` models, cross-cohort non-assignment Spearman rank correlations, and discovery-hit direction concordance are computed as in discovery. For shared-gene meta-analysis (`12_combined_program_summary.csv`), donor program scores are recomputed on the intersection of filtered genes across both cohorts and evaluated with combined sign-flip and paired *t* tests alongside Fisher's combined `limma::fry` *P* values ($\chi^2_4$).

## 6. Software environment

Analyses were run in R `4.6.2` with `Seurat 5.6.0`, `edgeR 4.10.7`, and `limma 3.68.5`. Complete session metadata are archived in [`docs/sessionInfo.txt`](sessionInfo.txt) and [`results/tables/11_sessionInfo.txt`](../results/tables/11_sessionInfo.txt).
