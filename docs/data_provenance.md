# Data provenance

## Discovery cohorts

| Series | Libraries | Source |
|---|---|---|
| [GSE172577](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE172577) | Six primary oral-tongue tumours, SYSMH1–SYSMH6 | Peng et al., *Oral Oncology*, 2021, PMID 34044317 |
| [GSE193766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193766) | Two primary oral tumours and their paired adjacent samples | Chen et al., *Science Advances*, 2022, PMID 35613277 |

The authoritative labels are in [the sample manifest](../config/samples.csv). The local project donor codes P1–P8 follow that mapping. The four GSE193766 libraries have these accession and tissue assignments.

| GEO accession | GEO title | Project donor | Tissue |
|---|---|---|---|
| GSM5820668 | N1 | P7 | Adjacent |
| GSM5820669 | T1 | P7 | Tumour |
| GSM5820670 | N2 | P8 | Adjacent |
| GSM5820671 | T2 | P8 | Tumour |

The preceding atlas used the opposite tissue assignments in both pairs. This revision refreshes the labels using accession and regenerates the plots and tables that encode tissue. The use of Harmony by library and shared gene labels means a tissue-only relabel can reuse the existing clustering and integer count snapshot.

## Count inputs

GEO supplies accession-level MTX and TSV files. The original atlas analysis used project-local derivatives named `<GSM>_seurat.rds` and `<GSM>_metaData.rds`. The snapshots were provided in the project workspace and their file hashes are in [data/checksums.csv](../data/checksums.csv). This naming convention refers to the local derivatives. For example, the GEO supplement for GSM5258385 enumerates the corresponding SYSMH1 MTX and TSV files. [GEO GSM5258385](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5258385)

The `counts` layer in each local Seurat snapshot contains integer UMIs. The similarly named `*_mat.rds` files contain previously transformed values and are excluded from the count model. The patient-level results use the curated stage 05 DC snapshot derived from the counts. Direct construction from the GEO matrices may change QC and clustering, since the two input forms have not been compared barcode by barcode. A rebuilt atlas requires fresh review of cluster-to-subtype mappings.

The tracked [primary counts](../results/tables/11_primary_pseudobulk_counts.tsv.gz) contain one UMI profile per eligible primary tumour and subtype. Each profile aggregates at least 10 cells. The [library map](../results/tables/11_primary_pseudobulk_libraries.csv) records sample accession, subtype, patient, cell count and UMI total. The [snapshot ledger](../results/tables/11_input_snapshot.csv) records the hash and dimensions of the discovery input. The original raw files and cell-level objects are excluded from Git.

## External validation series

[GSE181919](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181919) is the Choi et al. head-and-neck squamous-cell carcinoma study. For oral-cavity validation, use the ten oral-cavity primary tumours selected by the source annotation, omitting non-oral sites and other tissue categories. The deposited annotation identifies 10,397 cells, including 688 DCs, across these tumours. The [metadata inventory](../results/tables/12_cohort_inventory.csv) provides cohort membership by study pseudonym and library. The [input ledger](../results/tables/12_inputs.csv) records the file state and hash. [Choi et al.](https://pubmed.ncbi.nlm.nih.gov/36828832/)

## Use of human data

This project performs secondary analysis of publicly deposited, de-identified data. Repository outputs use library and donor study identifiers. The underlying data follow the source studies’ terms.
