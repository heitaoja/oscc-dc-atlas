# 数据来源与可复现性（Data provenance）

## 1. 数据集

| GEO | 样本数 | 组织 | 原文 | 用途 |
|---|---|---|---|---|
| [GSE172577](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE172577) | 6 | 舌部 OSCC 原发灶 | Peng Y. *et al.* Oral Oncol 2021 (**PMID 34044317**), TCF7+ T 细胞与 TLS | 主队列，全部为肿瘤 |
| [GSE193766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193766) | 4（2 癌 + 2 癌旁） | 口腔癌及配对正常组织 | Chen X. *et al.* Sci Adv 2022 (**PMID 35613277**), γδT17 / c-Maf | 提供肿瘤-癌旁配对（仅 2 对） |

两套数据都是 10x Genomics 3' 建库、GRCh38 比对，`.counts` 矩阵以 Cell Ranger 输出为准。
**注意**：两篇原文的主题分别是 T 细胞和 γδ T 细胞，**都不是为 DC 设计的**；
DC 在全部细胞里只占约 2%，这决定了本项目的所有 DC 结论都是"小样本、探索性"的。

## 2. 样本表

| 样本 | 数据集 | 分组 | 患者 |
|---|---|---|---|
| GSM5258385–GSM5258390 | GSE172577 | tumour | P1–P6 |
| GSM5820668 | GSE193766 | tumour | P7 |
| GSM5820669 | GSE193766 | normal | P7（配对） |
| GSM5820670 | GSE193766 | tumour | P8 |
| GSM5820671 | GSE193766 | normal | P8（配对） |

> 原论文摘要写"10 例 OSCC 样本"，实际是 **8 例患者、10 个文库，其中 2 个是癌旁正常组织**。

## 3. 文件说明（重要）

每个 GSM 在 GEO 附件里有 5 个对象。**用哪个、为什么**：

| 文件 | 内容 | 本项目是否使用 |
|---|---|---|
| `<GSM>_seurat.rds` | Seurat 对象；`counts` 层是**整数 UMI** | ✅ **唯一计数来源** |
| `<GSM>_mat.rds` | genes × cells，**log1p(CPM/1e4)**，非整数 | ❌ 不可当日志计数用 |
| `<GSM>_metaData.rds` | 每个 barcode 的 nGene/nUMI/percent.mito 及 QC 标记 | ✅ 仅用于交叉核对 |
| `<GSM>_sce.rds` | SingleCellExperiment，内容同 seurat | ❌ 重复 |
| `<GSM>_diagnostics.rds` | 质控诊断 list | ❌ 仅参考 |

`_mat.rds` 的判定依据（`R/01` 注释）：
取 GSM5820668 第 1 个细胞，`cor(_mat[,1], log1p(CPM10k)) = 0.9994`，最大值 12.65、非整数；
而 `_seurat.rds` 的 counts 最大 36,453 且全为整数。原论文脚本若把 `_mat` 当 counts 用，
会把已经标准化过的值再标准化一次。

## 4. 如何获取（重建环境）

```bash
# 项目内约定的目录结构
data/raw/GSE172577/<GSM>_seurat.rds
data/raw/GSE193766/<GSM>_seurat.rds
```

GEO 附件可以直接下载：

```bash
mkdir -p data/raw/GSE172577 data/raw/GSE193766
# GSE172577 的 6 个 GSM，GSE193766 的 4 个 GSM，仅下载 *_seurat.rds
# 例：
#   wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE172nnn/GSE172577/suppl/GSE172577_RAW.tar
#   tar -xf GSE172577_RAW.tar
```

> 若从 raw 10x 矩阵（`barcodes/features/matrix`）重建，`R/01` 需要把读取函数换成
> `Read10X()`；本仓库使用作者提供的 Seurat 对象以保持与论文数据的同一性，
> 读取函数集中在 `read_counts()` 一处，替换成本很低。

## 5. 校验

`data/checksums.csv` 记录了本仓库使用过的 20 个文件的 MD5。
复现时先核对：

```bash
cd data/raw && md5sum -c ../checksums.csv
```

## 6. 许可与伦理

两套数据均为 GEO 公开的人类去标识化数据，遵循各自原文的数据使用声明；
本项目不再做个体再识别，也不发布任何可回溯到患者的原始序列。
