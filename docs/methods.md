# 方法（Methods）

> 本文件描述的每一个参数都写在 `R/00_setup.R` 的 `PARAMS` 里，脚本中不出现魔法数字。
> 复现命令：`make all`（或 `make 01_load_qc 02_integrate_cluster ...`）。

## 0. 总览

```
data/raw/<GSE>/<GSM>_seurat.rds  (raw UMI counts)
        │
        ├─ 01  基因 ID 映射 → 稀疏化 → 每样本 QC → scDblFinder → 合并
        │
        ├─ 02  LogNormalize → 2000 HVG(vst) → Scale → PCA(50)
        │      → Harmony(batch = sample) → Leiden(0.5) → UMAP
        │
        ├─ 03  SingleR ×3 参考集（Blueprint / Monaco / HPCA）+ cluster marker
        │      → config/annotation_major.csv → 主细胞类型
        │
        ├─ 04  取 DC → 重新标准化/降维 → Harmony → Leiden(0.4) → 亚群 marker
        │
        ├─ 05  config/annotation_dc.csv → DC 亚型 → FindAllMarkers
        │      → migDC vs cDC2 差异分析（RNA assay, Wilcoxon, BH）
        │
        ├─ 06  GO-BP / KEGG（上调、下调分开）+ fgsea Hallmark
        ├─ 07  slingshot 拟时序（根节点固定为 cDC2）
        ├─ 08  decoupleR + MSigDB C3(TFT) → TF 活性
        └─ 09  DC 浸润比例 / 亚型构成 / cDC1:cDC2 平衡
```

## 1. 数据读取与基因 ID

- 计数来源：`<GSM>_seurat.rds` 的 `counts` 层（整数 UMI）。见 `docs/data_provenance.md` §3。
- Ensembl ID → HGNC symbol：`org.Hs.eg.db`（`mapIds`, `multiVals="first"`）；
  去版本号后映射，映射失败的基因丢弃。
- 同名基因（多个 Ensembl 映射到同一 symbol）用**稀疏聚合矩阵求和**，不做"取最大"，
  避免丢失计数。
- 两套数据的特征集不同（GSE172577: 22,161；GSE193766: 33,514），取交集后合并。

## 2. 质量控制

| 规则 | 阈值 | 依据 |
|---|---|---|
| 最少基因数 | `nFeature_RNA ≥ 200` | 空液滴/破损细胞 |
| 最少 UMI | `nCount_RNA ≥ 500` | 同上 |
| 最多基因数 | `nFeature_RNA ≤ 6000` | 双细胞/多细胞 |
| 自适应上界 | `log10(nFeature) ≤ median + 3×MAD`（每样本） | 各样本测序深度差异大，固定阈值会切掉高深度样本的真实细胞 |
| 线粒体 | `percent.mt ≤ 15%` | 口腔黏膜组织解离应激，阈值比血样宽 |
| 血红蛋白 | `percent.hb ≤ 3%` | 去除红细胞 |
| 双细胞 | `scDblFinder`（每样本独立，默认预期比例） | 组织解离样本必做 |

**实际结果**（`results/tables/01_cell_numbers.csv`）：

| 阶段 | 细胞数 | 占比 |
|---|---|---|
| Cell Ranger barcode | 113,448 | — |
| 数据提供方 QC 后（GEO 对象） | 108,569 | 95.7% |
| 本流程 QC 后 | 95,328 | 84.0% |
| 去双细胞后 | **85,164** | **75.1%** |

各文库 QC 去除率 6.3%–24.0%（GSM5258390 最高，其中位线粒体 7.6%），
双细胞率 8.2%–11.6%。原论文没有做双细胞去除，也没有做样本内自适应阈值；
这是重建版与原版最大的差别之一。

## 3. 标准化、整合、聚类

1. `LogNormalize`，scale factor 1e4。
2. 高变基因 2000 个（`vst`），仅用于降维，不改变表达矩阵。
3. `ScaleData` → `RunPCA(50)`。
4. **Harmony** 以 `sample` 为批次变量（10 个文库），`theta = 2`，使用前 30 个 PC。
   - 选 Harmony 而非 CCA anchors：细胞数 >7 万时 CCA 的内存与时间成本高，
     且 Harmony 对"同一组织、不同文库"的批次校正足够。
   - 判别标准：校正后 UMAP 上各样本应混合，而生物学亚群仍分离（`02_umap_by_sample.png`）。
5. `FindNeighbors(dims = 1:30)` → `FindClusters(resolution = 0.5, algorithm = 4)`（Leiden）。
6. `RunUMAP(dims = 1:30, seed = 20240601)`。

## 4. 主细胞类型注释（三级证据）

- **参考集投票**：SingleR × {Blueprint+ENCODE, Monaco immune, Human Primary Cell Atlas}，
  三个参考集的 `label.main` 分别产出，按 cluster 汇总（`results/tables/03_singler_labels.csv`）。
- **marker 证据**：`FindAllMarkers`（Wilcoxon，`min.pct = 0.1`，`logFC > 0.25`），
  每个 cluster 的前 10 个 marker（`03_cluster_markers.csv`）。
- **人工裁决**：两者一致才落定；不一致的 cluster 在 `config/annotation_major.csv`
  的 `evidence` 列写明理由。
- 与原论文的差别：原版用"单个基因名"当细胞类型（如 `cluster 5 = CDC20`、`cluster 23 = H2AC14`），
  且 SingleR 参考集被 `sub("-.*","")` 破坏；本版不允许"基因名 = 细胞类型"。

## 5. DC 亚群

- 取 `cell_type == "DC"` 的细胞（约 2%），`DietSeurat` 后重新走
  `NormalizeData → HVG(1500) → Scale → PCA(30) → Harmony(sample) → Leiden(0.4) → UMAP`。
- 亚群命名依据表 1（与原论文同一套 marker，补上 pDC 的正确 marker）：

| 亚型 | 核心 marker | 参考 |
|---|---|---|
| cDC1 | `CLEC9A`, `XCR1`, `BATF3`, `IRF8` | Sánchez-Paulete 2021 |
| cDC2 | `CD1C`, `CLEC10A`, `SIRPA`, `FCER1A` | Heidkamp 2020 |
| migDC（成熟/迁移型） | `CCR7`, `LAMP3`, `CCL19`, `CCL17`, `CD274` | Maier 2020; Zhang 2019 (LAMP3+ DC) |
| pDC | `LILRA4`, `CLEC4C`, `IL3RA`, `TCF4`, `IRF7` | — |
| LC-like | `CD207`, `CD1A` | — |

- **pDC 不再用 `CXCR3`**（原论文表 2-1 的错误）。

## 6. 差异表达

- 全部在 `RNA` assay（log-normalized）上做，**不在 integrated/harmony 空间做**。
- `FindMarkers(test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0.25)`。
- 显著性：BH 校正 `p_val_adj < 0.05`；效应量：`|avg_log2FC| ≥ 0.5`。
- 主对比：`migDC` vs `cDC2`（双向），另附各亚群 `FindAllMarkers`（仅上调）。
- 不把"细胞"当独立重复去推断患者水平差异（见 §10）。

## 7. 富集分析

- **方向敏感**：上调、下调基因集**分别**做 GO-BP 与 KEGG，不再合并（原论文的问题）。
- 背景基因集：本次检测到的全部基因（`universe`），不用全基因组。
- `clusterProfiler::enrichGO`（BH，`p < 0.05`, `q < 0.25`）、`enrichKEGG`（BH）。
- **GSEA**：`fgsea` + MSigDB Hallmark，排序统计量
  `sign(log2FC) × -log10(p_adj)`，避免人为 fold-change 阈值。

## 8. 拟时序

- 输入：cDC2 + migDC（可再加 cDC1），仅在 DC 子集内建图。
- `slingshot`，`clusterLabels = dc_subtype`，`start.clus = "cDC2"`（**写死在代码里**，可复现）。
- 伪时间 = 各 lineage 的均值；不可达细胞记为 `NA`（不填 0，避免污染）。
- 基因随伪时间变化：HVG 的 Spearman 相关 + 分箱均值曲线（20 箱）。
- **解释边界**：伪时间只是转录排序，不等于分化/谱系关系；原论文的"C2→C3 分化"在本版
  表述为"沿成熟轴的转录连续谱"。

## 9. TF 活性

- MSigDB C3（TFT:GTRD / TFT_Legacy）作为 TF–靶基因网络。
- `decoupleR::run_ulm`（单变量线性模型）计算每个细胞的 TF 活性，按亚型取均值。
- 替代原计划但未实现的 SCENIC（GENIE3 + cisTarget 需要额外数据库与数小时计算）。

## 10. 统计与报告规范

- 所有图都同时输出 PNG(300 dpi) 与 PDF；文字结论必须在 `results/tables/` 有对应文件。
- 随机种子固定（`SEED = 20240601`），所有随机过程显式传参。
- 版本信息：`make session` 生成 `docs/sessionInfo.txt`。
- 样本水平结论（如肿瘤 vs 癌旁）只做描述性统计，n = 2 对，不做检验。

## 11. 与原论文的方法学差异（速查）

| 环节 | 原论文 | 本版 |
|---|---|---|
| QC | 三套参数互相矛盾，无 doublet | 单一定义 + 自适应上界 + scDblFinder |
| 整合 | CCA anchors（dims 30/50 混用） | Harmony（batch = sample） |
| 注释 | 单基因名命名 + 被破坏的 SingleR | 三参考集 + marker + 人工裁决表 |
| malignant 细胞 | 声称做 CNV，实际未做 | 明确列为限制（未做 CNV） |
| DE | integrated assay | RNA assay |
| 富集 | 上下调混合，单张图 | 分方向 + universe + GSEA |
| 拟时序 | 鼠标点根节点 | 代码固定根节点 |
| TF | 计划 SCENIC，未运行 | decoupleR + C3 |
| 随机性 | 无种子 | 全程固定 |
