# Data Preparation

## 1. Quick verification of the primary model

The primary donor-paired `edgeR` differential expression model (`migDC` vs `cDC2`) can be verified directly from the tracked aggregate counts (`results/tables/11_primary_pseudobulk_counts.tsv.gz`) without downloading raw cell-level matrices:

```bash
Rscript -e 'install.packages("BiocManager"); BiocManager::install(c("edgeR", "statmod"))'
Rscript R/verify_primary.R
```

## 2. Full discovery pipeline (`GSE172577` & `GSE193766`)

Install required R/Bioconductor dependencies (tested on R `4.6.2`, Seurat `5.6.0`, edgeR `4.10.7`, limma `3.68.5`):

```bash
Rscript -e 'install.packages("BiocManager"); BiocManager::install(c(
  "Seurat", "harmony", "SingleR", "celldex", "scDblFinder", "clusterProfiler",
  "org.Hs.eg.db", "fgsea", "slingshot", "decoupleR", "msigdbr", "edgeR", "limma",
  "statmod", "data.table", "future", "future.apply", "dplyr", "tidyr", "ggplot2",
  "ggrepel", "patchwork"
))'
```

Place the 10 `<GSM>_seurat.rds` and `<GSM>_metaData.rds` files under `data/raw/GSE172577/` and `data/raw/GSE193766/` according to [`config/samples.csv`](../config/samples.csv) and verify integrity from `data/raw/`:

```bash
md5sum -c ../checksums.csv
```

If reusing pre-built `.rds` objects stored in an external directory, set `OSCC_OBJECT_DIR=/path/to/objects` before running downstream scripts (`11_donor_pseudobulk.R`, `13_tme_covariation.R`).

## 3. External validation cohort (`GSE181919`)

Download `GSE181919_Barcode_metadata.txt.gz` (`339,614` bytes) and `GSE181919_UMI_counts.txt.gz` (`127,878,601` bytes) from [GEO GSE181919](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181919) into `data/raw/GSE181919/`, then run:

```bash
./run_all.sh 12
```

`R/12_external_validation.R` streams only the target oral-cavity myeloid APC columns from the compressed matrix and caches the sparse integer `dgCMatrix` in `results/objects/12_gse181919_oral_myeloid_counts.rds`.
