# Methods

## Discovery atlas

The atlas uses six GSE172577 primary oral tumours plus two GSE193766 tumour–adjacent pairs. The 10 libraries represent eight donors. The sample manifest uses GEO accession as the join key and records the GEO title, cohort, tissue and project donor ID. Rereading cached atlas objects refreshes those fields. The atlas stages use the same cell selection and UMI counts as the preceding release.

1. Read integer RNA-layer UMI counts from the project’s harmonised input snapshots. Map Ensembl IDs to symbols with `org.Hs.eg.db`, sum counts for duplicate symbols and use the intersection across samples.
2. Apply the QC thresholds in `R/00_setup.R` per library. The thresholds are 200–6,000 genes, at least 500 UMIs, mitochondrial fraction up to 15%, haemoglobin fraction up to 3%, and a sample-adaptive upper feature limit. Call and remove doublets with scDblFinder.
3. Log-normalise, select 2,000 highly variable genes, compute PCA, adjust the embedding by library with Harmony, and use Leiden clustering and UMAP. The integer RNA counts are retained separately.
4. Annotate major lineages with cluster markers and SingleR profiles. The DC candidate pool includes the major labels DC and pDC. Recluster this pool and use `config/annotation_dc.csv` to assign six retained DC states and two excluded clusters.
5. Render composition with two explicit denominators. Candidate DCs and retained DCs each use total QC-passing atlas cells as denominator. Subtype fractions use retained DCs. Fill unobserved subtype counts with zero for each library and use paired displays for the two tissue pairs.

The legacy cell-level Wilcoxon tests quantify expression differences between clusters. Seurat’s adjusted P values use its built-in Bonferroni correction. The corresponding overrepresentation tests use BH adjustment and the assay feature universe. Legacy GSEA ranks genes by average log₂ fold change. Slingshot orders cDC2 and migDC by transcriptomic similarity conditional on a cDC2 root. The decoupleR scores are expression-based summaries of candidate regulatory activity. These legacy outputs describe the discovery atlas.

## Donor-paired expression

Stage 11 sums the integer RNA counts by library and subtype. Only primary tumours enter the focal comparison. The primary model requires both cDC2 and migDC in the same donor’s tumour, with at least 10 cells per subtype and at least four eligible donors. cDC2 is the reference. The model uses expression filtering by subtype with `edgeR::filterByExpr`, TMM normalisation, and a robust edgeR negative-binomial quasi-likelihood fit with design `~ patient + dc_subtype`. The subtype coefficient is tested with the QL F test. BH adjustment uses the genes tested in each contrast. The analysis is conditional on the defined subtype labels and retrospective selection of cohorts and analysis settings. [Squair et al.](https://doi.org/10.1038/s41467-021-25960-2) [edgeR](https://bioconductor.org/packages/edgeR)

Two secondary contrasts test migDC versus cDC1 and cDC1 versus cDC2, using a five-cell threshold and the same four-donor minimum. Every sensitivity fit records eligible donors, testable genes and BH-adjusted results. The reported effects use the same direction convention in all models.

| Check | Method |
|---|---|
| Donor leverage | Refit the primary model after removing each donor. Compare log₂ fold changes over the intersecting tested gene universe and record direction stability. |
| Subtype cell support | Fit thresholds of 3, 5, 15 and 20 cells per subtype where at least four donors qualify. |
| Normalisation | Refit the primary contrast with TMMwsp. |
| Cell-count imbalance | Randomly subsample both subtypes in each donor to the same cell count, 40 times. Use the primary tested gene universe, recalculate TMM, and summarise paired program and panel effects across runs. The spread reflects subsampling sensitivity. |
| Subtype assignment | Reclassify the original cDC1, cDC2 and migDC pool with the fixed three-gene gates defined in `config/external_dc_gates.csv`. Require detection of at least two genes in one gate and reject ambiguous cells. Repeat at 5 and 10 cells per subtype. |
| Cohort representation | Summarise the program effects within each of the two source series. These results describe cohort consistency in the integrated discovery atlas. |

The tracked pseudobulk count matrix includes the ten eligible primary aggregates. `R/verify_primary.R` checks the library mapping and counts, then verifies the fitted log₂ fold changes, P values and FDR against the tracked results.

## Expression programs and Hallmark pathways

The program definitions in `config/dc_programs.csv` list maturation, NF-κB feedback, MHC class II and regulatory expression genes separately. The lists use literature and curated pathway knowledge. The feedback panel addresses an observation from the original exploratory atlas. Every panel excludes genes explicitly used in the focal subtype assignment, including the external gates. The regulatory panel is evaluated only if at least 60% of its genes pass the expression filter. The same requirement applies to all program scores. [Maier et al.](https://pubmed.ncbi.nlm.nih.gov/32269339/)

Within each donor, take the mean TMM log₂ CPM difference (prior count 1) across retained program genes. Report donor values, an unweighted mean and descriptive 95% t intervals. Use an exact paired, two-sided sign-flip test and BH adjustment across the defined programs. Five donors give a minimum attainable two-sided P value of 0.0625. The t interval adds an approximate estimate of effect spread.

Test Hallmark enrichment using CAMERA with the paired design and the gene set cache `data/genesets/h.all.symbols.gmt`. Exclude explicit assignment genes from the tested universe, retain the TMM effective library sizes, and test sets with at least ten remaining genes. Use the default inter-gene correlation of 0.01 and BH adjustment. CAMERA direction denotes differential expression of a gene set relative to the rest of the tested universe. [Wu and Smyth](https://pubmed.ncbi.nlm.nih.gov/22638577/)

## External GSE181919 validation

Stage 12 uses primary tumours from the oral cavity (`tissue.type=CA`, `subsite=OC`) and the authors’ broad dendritic-cell annotation. Match metadata to UMI columns by unique barcode. Subtype classification uses only the fixed gates, with ambiguous and unassigned cells counted and excluded from paired testing. Sensitivities use the five-cell minimum and a one-marker rule. The model uses one observation per donor–subtype. A comparison runs once at least four donors meet its cell threshold. [GEO GSE181919](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181919)

The outputs report effect-size rank correlation on common testable genes, direction concordance among discovery FDR hits, FDR support in the external cohort, and donor scores for shared program genes. The cross-cohort comparisons of genes omit the assignment panel. For shared program scores both series use the same genes. Scope and estimates are reported at each fitted threshold.

The current real-data run generated the metadata inventory. The UMI supplement is required for the external expression analysis. The implementation was exercised end to end with synthetic, independently shuffled barcode rows and columns to test sample alignment, classification, fitting and shared-program output.

## Software

The discovery analysis used R 4.6.2, Seurat 5.6.0 and edgeR 4.10.7. Package and platform details are in [discovery session information](../results/tables/11_sessionInfo.txt) and [the original atlas session](sessionInfo.txt).
