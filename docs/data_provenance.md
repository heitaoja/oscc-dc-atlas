# Data Provenance and Cohort Accessions

## 1. Discovery cohorts

| GEO Series | Libraries | Tissue Scope | Source Publication |
|---|---|---|---|
| [GSE172577](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE172577) | 6 (`GSM5258385`–`GSM5258390`, `SYSMH1`–`SYSMH6`) | Primary oral-tongue squamous cell carcinoma | Peng *et al.*, *Oral Oncol* 2021 (PMID: 34044317) |
| [GSE193766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193766) | 4 (`GSM5820668`–`GSM5820671`, `N1`/`T1`, `N2`/`T2`) | 2 primary OSCC tumours and paired adjacent mucosa | Chen *et al.*, *Sci Adv* 2022 (PMID: 35613277) |

Authoritative sample-to-donor and tissue assignments are tracked in [`config/samples.csv`](../config/samples.csv):

| GEO Accession | GEO Title | Donor ID | Tissue |
|---|---|---|---|
| `GSM5258385`–`GSM5258390` | `SYSMH1`–`SYSMH6` | `P1`–`P6` | `tumour` |
| `GSM5820668` | `N1` | `P7` | `normal` (adjacent) |
| `GSM5820669` | `T1` | `P7` | `tumour` |
| `GSM5820670` | `N2` | `P8` | `normal` (adjacent) |
| `GSM5820671` | `T2` | `P8` | `tumour` |

All discovery libraries were generated using 10x Genomics 3' scRNA-seq aligned to GRCh38. Integer UMI counts are extracted exclusively from the `counts` layer of `<GSM>_seurat.rds`; normalized `<GSM>_mat.rds` matrices are excluded. MD5 checksums for all 20 discovery input files are recorded in [`data/checksums.csv`](../data/checksums.csv).

## 2. Independent validation cohort

| GEO Series | Selected Subset | Cell Counts | Source Publication |
|---|---|---|---|
| [GSE181919](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181919) | 10 primary oral-cavity tumours (`tissue.type == "CA"`, `subsite == "OC"`) | `10,397` total cells; `2,902` myeloid APCs (`688` `Dendritic.cells` + `2,214` `Macrophages`) | Choi *et al.*, *Nat Commun* 2023 (PMID: 36828832) |

File sizes and MD5 checksums for `GSE181919_Barcode_metadata.txt.gz` and `GSE181919_UMI_counts.txt.gz` are recorded in [`results/tables/12_inputs.csv`](../results/tables/12_inputs.csv), and sample-level cell inventories are in [`results/tables/12_cohort_inventory.csv`](../results/tables/12_cohort_inventory.csv).

## 3. Tracked source count data

The aggregated donor-by-subtype integer UMI matrix (`results/tables/11_primary_pseudobulk_counts.tsv.gz`) and library metadata (`results/tables/11_primary_pseudobulk_libraries.csv`) are tracked directly in the repository to enable immediate verification of the primary differential expression results via `Rscript R/verify_primary.R`.

## 4. Human data governance

This study relies solely on publicly available, de-identified human scRNA-seq datasets deposited in NCBI GEO. All outputs use anonymized study library and donor identifiers in accordance with the original submissions' data access terms.
