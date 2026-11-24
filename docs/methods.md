# Technical methods

## Discovery atlas and composition

The discovery atlas comprises six GSE172577 primary oral-tongue tumours and two GSE193766 primary tumour–adjacent pairs. The resulting 10 libraries represent eight donors. Sample accession is the join key for source series, patient and tissue. The GSE193766 assignments are GSM5820668 and GSM5820670 for adjacent tissue, and GSM5820669 and GSM5820671 for the matched tumours. Reading an atlas object refreshes these fields from the sample manifest. [Sample manifest](../config/samples.csv) · [Data provenance](data_provenance.md) · [Library composition](../results/tables/09_dc_fraction_by_sample.csv)

1. Import the integer RNA UMI layer of the harmonised sample snapshots. Map Ensembl IDs to symbols with `org.Hs.eg.db`, sum duplicated symbols and retain the shared feature universe across libraries.
2. Apply per-library QC using the saved parameters. The absolute feature bounds are 200–6,000, with a 500 UMI minimum, a mitochondrial fraction limit of 15% and a haemoglobin fraction limit of 3%. Apply the sample-adaptive feature cap in log space and use scDblFinder to remove doublets.
3. Log-normalise for feature selection and embedding, select variable genes, perform PCA, use library as the Harmony batch variable, and cluster with Leiden. The counts layer remains available for pseudobulk.
4. Annotate major cell types using cluster markers and SingleR. Recluster cells with broad DC or pDC labels and use the reviewed DC mapping to retain cDC1, cDC2, LC-like, migDC, pDC and proliferating DC. The LC-like annotation includes *CD207*, *CD1A* and *CD1E*.
5. Calculate candidate and retained DC fractions with all quality-controlled cells in each library as denominator. Use retained DCs as denominator for subtype proportions. Present patient-matched tumour–adjacent differences descriptively.

[QC parameter table](../results/tables/00_params_qc.csv) · [DC annotation](../config/annotation_dc.csv) · [Pooled DC composition](../results/tables/09_dc_subtype_pooled.csv) · [Paired tissue composition](../results/tables/09_paired_tissue_descriptive.csv)

The original cell-level analysis uses Seurat's Wilcoxon marker test and its Bonferroni-adjusted significance. Legacy gene-set ranking, pseudotime and transcription-factor scores are descriptive atlas outputs. Slingshot orders cells conditional on a selected root. [Cell-level results](../results/tables/05_de_migDC_vs_cDC2.csv) · [Pseudotime](../results/tables/07_pseudotime_by_cell.csv) · [Transcription-factor scores](../results/tables/08_tf_activity_by_subtype.csv)

## Stage 11 donor-paired differential expression

Use the Stage 05 RNA counts and reviewed subtype labels. Aggregate integer counts by library and subtype, then restrict fits to primary tumours. An eligible patient has both states at or above the contrast-specific cell minimum. A fit requires at least four eligible patients with one aggregate per state and tumour. The focal comparisons are migDC–cDC2, migDC–LC-like and LC-like–cDC2 at a 10-cell minimum per state, and migDC–cDC1 and cDC1–cDC2 at a five-cell minimum. [Donor eligibility](../results/tables/11_donor_eligibility.csv) · [Model settings](../results/tables/11_sensitivity.csv)

Use `edgeR::filterByExpr` with subtype groups, `min.count = 3` and `min.total.count = 10` for the sparse aggregates. Recompute library sizes after filtering and apply TMM normalisation. Use robust dispersion estimation and a robust negative-binomial quasi-likelihood fit with `~ patient + dc_subtype`. Test the subtype coefficient using the quasi-likelihood F test. The second named state is the reference. Apply BH adjustment over genes tested within each fit. [Shared model code](../R/pseudobulk.R) · [Primary results](../results/tables/11_de_migDC_vs_cDC2.csv) · [edgeR](https://bioconductor.org/packages/edgeR)

| Sensitivity | Implementation |
|---|---|
| Donor leverage | Refit after omitting each primary donor. Report shared-gene Spearman correlation and per-gene direction stability. |
| Subtype cell support | Test 3-, 5-, 15- and 20-cell thresholds under the same fitting minimum. Export eligibility even when a fit is unavailable. |
| Normalisation | Replace TMM with TMMwsp at the primary cell threshold. |
| Expression filtering | Keep the primary cell threshold and TMM in `strict_filter_min10`, with `min.count = 10` and `min.total.count = 15`. |
| Cell imbalance | Downsample each patient's two states to the smaller state size in 40 repeats. Recalculate TMM over the primary tested genes. Report variation across subsamples. |
| Marker-gate labels | Reclassify the original cDC1, cDC2 and migDC pool using the gates below and evaluate five- and 10-cell minima. |
| Discovery series | Summarise the program effects by source series within the integrated discovery atlas. |

[Leave-one-out](../results/tables/11_leave_one_donor_out.csv) · [Cell thresholds and count filters](../results/tables/11_sensitivity.csv) · [Equal-cell subsampling](../results/tables/11_balanced_cell_sensitivity.csv) · [Cohort summaries](../results/tables/11_program_by_cohort.csv)

The fixed gates are *LAMP3/CCR7/FSCN1* for migDC, *CD1C/CLEC10A/FCER1A* for cDC2 and *CLEC9A/XCR1/CADM1* for cDC1. A positive detection is any UMI count above zero. Assign the cell to a state when at least two marker genes are detected in exactly one gate. Record cells passing several gates as ambiguous and cells passing none as unassigned. Exclude both from gate-based fits. For concordance with the original fit, omit the union of explicit discovery assignment and gate genes. [Gate definitions](../config/external_dc_gates.csv) · [Label agreement](../results/tables/11_label_agreement.csv) · [Gate sensitivity](../results/tables/11_marker_gate_sensitivity.csv)

The saved count matrix and library map contain the eligible primary aggregates. `R/verify_primary.R` verifies their mapping and re-estimates the primary log₂ fold changes, P values and FDR. [Aggregate counts](../results/tables/11_primary_pseudobulk_counts.tsv.gz) · [Library map](../results/tables/11_primary_pseudobulk_libraries.csv) · [Primary verifier](../R/verify_primary.R)

## Program and pathway statistics

The program definitions contain four held-out gene lists for Maturation, NF-κB feedback, Regulatory expression and MHC class II. Explicit discovery assignment markers and gate genes are excluded from the lists. Program eligibility requires at least 60% of target genes to pass expression filtering. Calculate `edgeR::cpm(y, log = TRUE, prior.count = 1)` on the count-filtered, TMM-normalised libraries. Average the paired log₂ CPM differences across eligible program genes within each patient. The published summary reports patient counts, sign counts, mean effects and two-sided 95% *t* intervals. [Definitions](../config/dc_programs.csv) · [Donor effects](../results/tables/11_program_donor_effects.csv) · [Program summary](../results/tables/11_program_summary.csv)

The exact two-sided sign-flip test compares the absolute observed sum of patient differences against the full signed distribution. Its BH-adjusted result is the `fdr` column. The two-sided paired *t* test is a one-sample test of the same patient differences against zero, with outputs `p_paired_t` and `fdr_paired_t`. The self-contained `limma::fry` analysis operates on the filtered TMM log₂ CPM matrix using `model.matrix(~ patient + dc_subtype, y$samples)` and the final design coefficient. It reports the direction of the overall signal and directional and mixed-direction P values with separate BH adjustment. Each test family is evaluated within its contrast. The t and fry tests address related hypotheses using different aggregations. [Discovery fry](../results/tables/11_program_fry.csv) · [limma](https://bioconductor.org/packages/limma)

In the objective-gate fits, add *CD274* and *IDO1* to the Regulatory expression panel and test it alongside the four held-out programs. The nine gate genes exclude both additions. Discovery's original annotation used these genes, so this expanded panel is reserved for the gate-based discovery and external analyses. [Discovery gate-based programs](../results/tables/11_marker_gate_programs.csv) · [External program summary](../results/tables/12_program_summary.csv)

CAMERA uses the primary paired design and the count-filtered background after excluding the union of discovery assignment and gate genes. It retains the primary TMM effective sizes and uses locally cached MSigDB Hallmark sets with at least 10 assayed genes. The analysis uses the DGEList method's default intergene correlation of 0.01. BH adjusts across tested sets. Direction refers to enrichment relative to the rest of the tested genes. Figure 3C uses primary TMM effective sizes to show descriptive pairs for selected low-count genes that fall below the DE filter. [Hallmark results](../results/tables/11_paired_hallmark_camera.csv) · [CAMERA](https://pubmed.ncbi.nlm.nih.gov/22638577/) · [Figure 3](../results/figures/Fig3_programs_and_inflammation_paradox.png)

## Stage 12 independent GSE181919 validation

Select primary oral-cavity samples using `tissue.type == "CA"` and `subsite == "OC"`. The myeloid APC pool includes both `Dendritic.cells` and `Macrophages`. Extract its UMI columns from the deposited gzipped counts in row blocks, using a unique match to the barcode metadata after harmonising a terminal dot or hyphen before a numeric suffix. Validate gene identifiers, selected barcodes and the nonnegative integer counts. Store the selected sparse matrix with its barcodes and the source MD5 and byte size in `results/objects/12_gse181919_oral_myeloid_counts.rds`. Reuse the cache only when those checks agree. [GEO GSE181919](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181919) · [File ledger](../results/tables/12_inputs.csv) · [Sample inventory](../results/tables/12_cohort_inventory.csv)

Apply the fixed DC gates to the extracted matrix and tabulate exclusive subtype calls, ambiguity and nonassignment by sample and author label. The classification table shows both one- and two-marker detection cutoffs. The paired fits use the latter. CD163-negative status means zero observed *CD163* UMIs. The hybrid pool uses only author-labelled DCs for migDC and the CD163-negative myeloid pool for cDC2. Count filtering, TMM, fitting and program inference follow Stage 11. [Classification](../results/tables/12_classification.csv) · [Eligibility](../results/tables/12_donor_eligibility.csv)

| External fit IDs | States and cell support |
|---|---|
| `primary`, `min_5` | Myeloid APC migDC–cDC2, at least 10 or 5 cells per state. |
| `cd163_neg_min10`, `cd163_neg_min5` | CD163-negative myeloid migDC–cDC2, at least 10 or 5 per state. |
| `hybrid_author_migDC_vs_cd163neg_cDC2` | Author DC migDC and CD163-negative myeloid cDC2, at least 10 per state. |
| `author_dc_only_min10`, `author_dc_only_min5` | Author DC migDC–cDC2, at least 10 or 5 per state. |
| `author_dc_migDC_vs_cDC1_min3` | Author DC migDC–cDC1, at least 3 per state. |
| `migDC_vs_cDC1`, `cDC1_vs_cDC2` | Myeloid APC alternate lineage contrasts, at least 5 per state. |

[External methods script](../R/12_external_validation.R) · [Fitting summary](../results/tables/12_sensitivity.csv)

Each successful external fit is matched to its discovery lineage contrast. Calculate Spearman ρ over shared, tested non-assignment genes. The hit concordance denominator is shared genes with discovery FDR < 0.05. The numerator counts equal signs. External FDR support counts concordant discovery hits also meeting FDR < 0.05 externally. Insufficiently supported contrasts retain eligibility counts and missing fit statistics. [Cross-cohort genes](../results/tables/12_cross_cohort_genes.csv) · [Concordance summary](../results/tables/12_cross_cohort_summary.csv)

For every fitted migDC–cDC2 external model, recalculate the usual four program scores on the genes that passed the filters in that model and in the discovery primary contrast. Pool independent patients with equal weight and apply the mean, interval, sign-flip and paired *t* tests to the shared-gene differences. Separately combine the discovery and corresponding external directional fry P values with Fisher's chi-square method. Each fry run uses its own eligible genes, so the shared-gene pooled estimate and Fisher P value represent complementary summaries. Record the cohort-specific directions and use BH over the four combined P values in each model. [Shared-gene donor scores](../results/tables/12_shared_program_donor_effects.csv) · [External fry](../results/tables/12_program_fry.csv) · [Combined summary](../results/tables/12_combined_program_summary.csv)

## Stage 13 tumour microenvironment

Calculate composition across all discovery libraries by joining atlas annotations, retained DC subtypes and sample metadata. The immune denominator comprises T cells, Tregs, NK cells, B cells, plasma cells, macrophages, DCs, pDCs, mast cells and neutrophils. The T/NK denominator includes T cells, Tregs and NK cells. Use the complete retained DC pool as denominator for its subtype proportions. The profiles also report migDC, Treg, T-cell and macrophage fractions of the immune pool and Treg fractions of the T/NK pool. [Composition](../results/tables/13_tme_cell_proportions.csv)

In primary tumours, aggregate integer counts separately by library for retained DC, migDC, Treg, T cells and macrophages. Require at least three migDC cells and a nonempty aggregate for the other populations. Estimate TMM factors on each population's entire transcriptome across its available libraries, then extract log₂ CPM with prior count 1. The T-cell inhibitory module averages *PDCD1*, *CTLA4*, *TIGIT*, *LAG3* and *HAVCR2*. Combine these features with paired Maturation, Regulatory expression and NF-κB feedback scores for the primary and three-cell sensitivity thresholds. Stage 13 recomputes the three-cell paired scores with Stage 11 settings when that run's donor scores were not exported. [Patient profiles](../results/tables/13_tme_donor_profiles.csv) · [TME analysis script](../R/13_tme_covariation.R)

Test each specified composition, cross-population ligand–receptor and paired-program axis on finite patient pairs. Require at least three patients and both features to vary. Use two-sided Spearman `cor.test(..., exact = FALSE)` and BH adjustment over all tested axes together. In the audited output, 42 of 46 axes were testable across the tumour profiles. The four unavailable axes used the low-coverage Regulatory expression panel at the three-cell threshold. Retained-DC *CCL17*–Treg *CCR4* reached ρ = 0.93 and FDR = 0.0362. [All TME axes](../results/tables/13_tme_ligand_receptor_correlations.csv)

## Stage 14 composite figures and audit

The five manuscript figures use the reviewed tables, Stage 11 pseudobulks and the lightweight Stage 05 DC object. The plotting script uses ggplot2 and patchwork and saves 300 dpi PNG and Cairo PDF. [Figure script](../R/14_manuscript_figures.R) · [Fig1](../results/figures/Fig1_discovery_atlas.png) · [Fig2](../results/figures/Fig2_pseudobulk_de_and_robustness.png) · [Fig3](../results/figures/Fig3_programs_and_inflammation_paradox.png) · [Fig4](../results/figures/Fig4_lineage_and_tme_context.png) · [Fig5](../results/figures/Fig5_external_validation.png)

The audited primary fit covers five patients, 4,151 genes and 933 hits at FDR < 0.05. The external primary APC fit covers six patients and yields ρ = 0.671 and 93.4% concordance among shared discovery hits. The four-patient hybrid fit yields ρ = 0.751 and 96.8% hit concordance. In the primary and five-cell combined analyses, the three positive programs are concordant in 11/11 and 12/12 patients, respectively. [Discovery model counts](../results/tables/11_sensitivity.csv) · [External replication](../results/tables/12_cross_cohort_summary.csv) · [Combined analysis](../results/tables/12_combined_program_summary.csv)

Run commands from the repository root with the prepared inputs and reviewed atlas objects. The full audit refits the primary model, verifies the Stage 09 and 11–14 evidence, and checks the current manuscript and repository documentation. [Data setup](../data/README.md) · [Audit script](../R/99_verify_all.R)

~~~bash
Rscript R/verify_primary.R
./run_all.sh 09 11 12 13 14 verify_all
Rscript R/99_verify_all.R
~~~

Software versions for the reviewed analyses are recorded with the results. The run used R 4.6.2, Seurat 5.6.0, edgeR 4.10.7 and limma 3.68.5. [Discovery session](../results/tables/11_sessionInfo.txt) · [External session](../results/tables/12_sessionInfo.txt)
