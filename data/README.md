# Data preparation

## Refit the primary model

The main test can be reproduced directly with the public [primary aggregate counts](../results/tables/11_primary_pseudobulk_counts.tsv.gz) and [library map](../results/tables/11_primary_pseudobulk_libraries.csv).

~~~bash
Rscript -e 'install.packages("BiocManager"); BiocManager::install(c("edgeR","statmod"))'
Rscript R/verify_primary.R
~~~

## Rebuild the discovery atlas

Use a current R and matching Bioconductor release. The committed run used R 4.6.2 and Seurat 5.6.0.

~~~bash
Rscript -e 'install.packages("BiocManager"); BiocManager::install(c(
  "Seurat","harmony","SingleR","celldex","scDblFinder","clusterProfiler",
  "org.Hs.eg.db","fgsea","slingshot","decoupleR","msigdbr","edgeR","limma",
  "statmod","data.table","future","future.apply","dplyr","tidyr","ggplot2",
  "ggrepel","patchwork"))'
~~~

Place the harmonised `<GSM>_seurat.rds` and `<GSM>_metaData.rds` snapshots under `data/raw/<GSE>/`. Each Seurat snapshot supplies the integer RNA count layer. The metadata snapshot supplies barcode-indexed QC columns, including `nGene`, `nUMI` and `percent.mito`. The [sample manifest](../config/samples.csv) specifies the ten accessions. The previously prepared files can be verified with the following command from `data/raw`.

~~~bash
md5sum -c ../checksums.csv
~~~

The [GEO accessions](../docs/data_provenance.md) provide MTX and TSV supplements. A reconstruction from these files needs a `ReadMtx()` import with Ensembl IDs, barcode-level QC metadata, and a review of the cluster annotations after running the atlas. The versioned count aggregates support immediate reproduction of the primary model with the published atlas definition.

To reuse a local `05_dc_annotated.rds` that sits outside this checkout, point the stage 11 command at its containing directory.

~~~bash
OSCC_OBJECT_DIR=/path/to/results/objects Rscript R/11_donor_pseudobulk.R
~~~

## External GSE181919

Download the two supplements from [the GEO series page](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181919) and keep the original compressed files in `data/raw/GSE181919/`.

| File | GEO-listed size |
|---|---:|
| `GSE181919_Barcode_metadata.txt.gz` | 339,614 bytes |
| `GSE181919_UMI_counts.txt.gz` | 127,878,601 bytes |

~~~bash
Rscript R/12_external_validation.R --inventory-only
./run_all.sh 12
~~~

The first command requires the metadata file alone. The second requires both files and the stage 11 output tables. It streams the genes-by-cells TSV, so the complete dense gene-by-cell table stays off memory.
