# data/

这个目录**不进入 git**（见根目录 `.gitignore`）。目录约定：

```
data/
├── raw/                     # 原始数据（需自行下载，约 1.1 GB）
│   ├── GSE172577/
│   │   ├── GSM5258385_seurat.rds
│   │   └── ...
│   └── GSE193766/
│       ├── GSM5820668_seurat.rds
│       └── ...
├── genesets/                # 基因集缓存（可提交，约 3 MB）
│   ├── h.all.symbols.gmt              # MSigDB Hallmark v2024.1.Hs
│   ├── c3.tft.v2024.1.Hs.symbols.gmt  # MSigDB C3 TFT
│   └── kegg_hsa_pathway.gmt           # 由 rest.kegg.jp 构建（06 脚本自动生成）
└── checksums.csv            # 原始文件 MD5
```

## 如何准备 `raw/`

只需要每个 GSM 的 **`*_seurat.rds`**（内含整数 UMI 计数）。
`*_mat.rds` 是 log1p(CPM)，不能当计数用——原因见
[`docs/data_provenance.md`](../docs/data_provenance.md) §3。

从 GEO 下载：

```bash
mkdir -p data/raw/GSE172577 data/raw/GSE193766

# GSE172577（6 个 GSM）
wget -P /tmp https://ftp.ncbi.nlm.nih.gov/geo/series/GSE172nnn/GSE172577/suppl/GSE172577_RAW.tar
tar -xf /tmp/GSE172577_RAW.tar -C data/raw/GSE172577 --wildcards '*_seurat.rds'

# GSE193766（4 个 GSM）
wget -P /tmp https://ftp.ncbi.nlm.nih.gov/geo/series/GSE193nnn/GSE193766/suppl/GSE193766_RAW.tar
tar -xf /tmp/GSE193766_RAW.tar -C data/raw/GSE193766 --wildcards '*_seurat.rds'
```

校验：

```bash
cd data/raw && md5sum -c ../checksums.csv
```

> 如果 GEO 附件里没有 `_seurat.rds`（不同批次的补充文件可能不同），
> 也可以下载 10x 三件套后用 `Read10X()` 读入；
> 只需改 `R/01_load_qc.R` 里的 `read_counts()` 一个函数，其余步骤不变。
