# 局限性与"下一步"（Limitations & next steps）

> 写这一页的目的不是自谦，而是**在面试里主动划清边界**。
> 凡是这里写了的，都是我已经想到、并且能解释为什么没做的；凡是没写又被人问到的，才叫漏洞。

## 1. 数据层面（最硬的天花板）

| 限制 | 具体影响 | 缓解/替代 |
|---|---|---|
| 两套数据都是**为 T 细胞设计**的（GSE172577: TCF7+ T/TLS；GSE193766: γδT17/c-Maf） | DC 只占约 1.7%，稀有亚群（pDC、LC-like）细胞数可能 < 100 | 亚群结论一律标注细胞数；稀有亚群只做描述，不做统计推断 |
| 8 例患者，癌旁仅 2 例 | 任何"肿瘤 vs 癌旁"的差异都无法做统计检验 | 只报告描述性比例 + 配对连线，不下结论 |
| 无临床信息（分期、HPV、预后） | 无法做生存/分层分析 | 明确写"不做预后推断" |
| 无配对 TCR/BCR、无 CNV 数据 | 无法确认 DC 与 T 细胞的克隆互作、无法区分恶性上皮 | 列为 future work |

## 2. 方法层面（我选择不做的）

| 未做 | 为什么没做 | 如果要做 |
|---|---|---|
| **inferCNV / CopyKAT** 恶性细胞判定 | 需要基因位置文件 + 数十分钟计算，且本项目的结论只涉及 DC（免疫细胞），不依赖恶性细胞判定 | `infercnv` + `gene_order_file`，或 `CopyKAT`；用上皮 cluster 做 reference |
| **SoupX / CellBubble** ambient RNA 校正 | 没有 raw（未过滤）barcode 矩阵，无法估计污染率 | 从 Cell Ranger 的 `raw_feature_bc_matrix` 重跑 |
| **RNA velocity (scVelo)** | 需要 spliced/unspliced 计数，GEO 只给了 counts | 重新比对（velocyto/STARsolo） |
| **CellChat / CellPhoneDB** 细胞通讯 | 需要完整的配受体库与置换检验，DC 细胞数少时假阳性高；且原论文的 CCL17-CCR4 结论需要 T 细胞侧数据配合 | `liana` 或 `CellChat`，但必须报告置换 p 值与细胞数 |
| **scenic / GENIE3** | 需要 cisTarget 数据库（数 GB）与数小时计算 | 已用 `decoupleR + MSigDB C3(TFT)` 做轻量替代 |
| **空间转录组验证** | 公共数据中无配对的 OSCC Visium | 可找同癌种 Visium 做 deconvolution（cell2location/RCTD） |
| **湿实验验证**（流式/IHC/qPCR） | 纯干实验项目 | 这是所有"纯生信"项目的共同短板；面试时主动说明，并给出验证方案 |

## 3. 统计层面

1. **细胞当重复的陷阱**：`FindMarkers` 的 p 值把每个细胞当独立观测，几万个细胞会让
   极小的差异也显著。因此本项目同时报告**效应量**（`avg_log2FC`、`pct.1 - pct.2`），
   并在报告里只用"效应量 + 方向"下结论。
2. **没有 pseudobulk**：n = 8 患者时，正确的做法是把每个患者的亚群表达加和成 pseudobulk
   再做 `limma-voom`/`DESeq2`。本项目在 `docs/next_steps.md` 里给出实现方案，
   但因为正常组织只有 2 例，pseudobulk 的统计功效依然很低。
3. **亚型比例**：只有 2 对肿瘤/癌旁，不做检验；只画配对点线图。
4. **多重检验**：所有差异分析用 BH；富集分析用 BH + q < 0.25；GSEA 用 `padj`。

## 4. 生物学解释层面

- **拟时序 ≠ 分化**：slingshot 给出的是转录相似性排序。原论文说"C2 向 C3 分化"，
  本版只说"沿成熟轴的转录连续谱"，并引用 Maier 2020 / Zhang 2019 的 mregDC/LAMP3+ DC
  作为独立支持；若要证明分化，需要 velocity 或谱系追踪。
- **TF 活性 ≠ TF 蛋白活性**：decoupleR 的 ULM 打分只反映靶基因集的整体表达，
  不含磷酸化、入核、辅因子信息。
- **通路富集受数据库版本影响**：KEGG 是动态数据库，`enrichKEGG` 的结果随抓取日期变化；
  `docs/sessionInfo.txt` 与运行日志记录了日期，复现时可能有细微差异。
- **CLEC10A → IDO1 → 色氨酸代谢**：原论文的核心假说之一，但在本项目的 DC 数据里
  `IDO1` 的表达极低（见 `results/tables/05_dc_subtype_dotplot` 对应数据），
  没有证据支持。本版**不保留**该结论。

## 5. 如果要继续做（按性价比排序）

1. **pseudobulk 差异分析**（1 天）：`decoupleR`/`muscat` 或 `limma-voom`，把细胞水平
   的假显著性去掉，这是最容易被面试官追问的点。
2. **恶性细胞判定**（1 天）：`inferCNV`，让"上皮/肿瘤细胞"的注释有依据。
3. **细胞通讯**（1 天）：`liana` + 置换检验，回答"migDC 是否在招募 Treg"。
4. **独立队列验证**（2–3 天）：用 HNSCC 的公共 scRNA/Visium 数据验证 migDC 比例与 CCR7 轴。
5. **可交互报告**（1 天）：`shiny` 或 `cellxgene` 部署，作为作品集的"可点开"部分。

## 6. 补充：两个数据集的特征集不对称（实现细节，但会影响结论）

GSE172577 的 Seurat 对象只有 **22,161 个 Ensembl ID**，GSE193766 有 **33,514 个**，
且前者是后者的子集。合并时取**交集**（得到 20,814 个符号），代价是：

- GSE193766 独有的基因被丢弃；
- *RACK1*、*ATP5F1E*、*ATP5MG*、*ATP5MC2* 等基因**不在 GSE172577 的特征集内**，
  因此不在合并矩阵中，原论文的"RACK1 + ATP 合酶模块"无法在本项目中复核。

为什么不取并集？并集会把"GSE172577 没测到"当成"表达为 0"，
在跨数据集差异分析中制造假信号。取交集是保守但正确的选择。
若必须评估这些基因，正确做法是**只在 GSE193766 的 4 个文库内单独分析**，
并在图注中注明"仅基于 2 例患者"。

复现检查：

```r
g1 <- readRDS("data/raw/GSE172577/GSM5258385_seurat.rds")
g2 <- readRDS("data/raw/GSE193766/GSM5820668_seurat.rds")
length(intersect(rownames(g1), rownames(g2)))   # 22161
c("RACK1","ATP5F1E","ATP5MG","ATP5MC2") %in% rownames(g1)  # 全 FALSE
```
