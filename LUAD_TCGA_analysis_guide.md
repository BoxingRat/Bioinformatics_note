# TCGA-LUAD 转录组差异分析说明（R）

## 分析目标
- 数据来源：TCGA 数据库，项目 `TCGA-LUAD`。
- 数据类型：RNA-seq `STAR - Counts`。
- 任务：
  1. 预处理（缺失值处理、标准化）；
  2. 差异表达基因筛选（阈值 `|log2FC| > 2` 且 `P < 0.05`）；
  3. 核心差异基因功能富集分析（GO / KEGG）。

## 脚本文件
- 主脚本：`TCGA_LUAD_DEG_enrichment.R`

## 运行前准备
建议使用 R >= 4.2，并安装以下包：

```r
install.packages(c("dplyr", "stringr", "ggplot2", "ggrepel", "tibble"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "TCGAbiolinks", "SummarizedExperiment", "edgeR", "limma",
  "clusterProfiler", "org.Hs.eg.db", "enrichplot"
))
```

## 运行方式
```bash
Rscript TCGA_LUAD_DEG_enrichment.R
```

## 下载太慢/中断时的替代方案（推荐 `gdc-client`）
脚本已内置下载参数，默认：
- `download_method <- "client"`
- `files_per_chunk <- 20`

即默认使用 `gdc-client` 下载（更适合大文件与不稳定网络，支持断点续传）。  
若本机未安装 `gdc-client`，脚本会报错提示并退出。

### 安装 `gdc-client`（示例）
1. 访问 GDC 官方下载页并下载对应系统的二进制程序。  
2. 解压后将 `gdc-client` 加入系统 `PATH`。  
3. 终端验证：
```bash
gdc-client --version
```

### 失败重试建议
- 保持 `download_method <- "client"`；
- 调小 `files_per_chunk`（如 10）减少单次批量请求压力；
- 中断后重新运行脚本，`gdc-client` 会利用已有文件继续下载。

### 若仍想使用 API
把脚本中的：
```r
download_method <- "client"
```
改为：
```r
download_method <- "api"
```

## 关键流程说明
1. **下载与读取数据**：使用 `TCGAbiolinks` 获取 LUAD 的 STAR counts。
   - 若 `SummarizedExperiment` 中存在多个 assay，脚本优先使用 `unstranded`。
2. **分组信息提取**：从 TCGA barcode 中提取样本类型，分为 `Tumor` 与 `Normal`。
3. **缺失值处理**：如存在 NA，按基因维度用中位数填补。
4. **标准化**：
   - 先过滤低表达基因（`filterByExpr`）；
   - 再进行 TMM 标准化（`calcNormFactors`）；
   - 导出 logCPM 标准化矩阵。
5. **差异分析**：`limma-voom` 建模，按 `Tumor vs Normal` 计算。
6. **DEG筛选阈值**：`abs(logFC) > 2` 且 `P.Value < 0.05`。
   - 额外输出 `gene symbol` 命名的差异基因表。
7. **核心差异基因定义**：显著 DEG 中按 `|logFC|` 取前 200（可按需求修改）。
8. **功能富集**：
   - GO（BP/CC/MF）；
   - KEGG 通路；
   - 输出表格与 dotplot 图。
9. **火山图标注**：按 `P.Value` 选取 top 20 显著差异基因并标注基因名。

## 输出结果（`results_luad/`）
- `TCGA_LUAD_SE.rds`：原始 SummarizedExperiment 对象。
- `LUAD_logCPM_normalized.csv`：标准化表达矩阵。
- `LUAD_DEG_all.csv`：全部基因差异分析结果。
- `LUAD_DEG_sig_logFC2_P0.05.csv`：显著差异基因结果。
- `LUAD_DEG_sig_logFC2_P0.05_symbol.csv`：`gene symbol` 命名的显著差异基因结果。
- `LUAD_core_DEG_top200.csv`：核心差异基因（默认 top200）。
- `LUAD_core_DEG_GO_enrichment.csv`：GO 富集结果。
- `LUAD_core_DEG_KEGG_enrichment.csv`：KEGG 富集结果。
- `LUAD_DEG_volcano.png`：差异火山图。
- `LUAD_core_DEG_GO_dotplot.png`：GO 富集图。
- `LUAD_core_DEG_KEGG_dotplot.png`：KEGG 富集图。

## 可选优化建议
- 若希望更稳健，建议在筛选时把 `P.Value` 改为 `adj.P.Val`（FDR 校正）。
- 可补充 GSEA（按全基因排序）以减少阈值依赖。
- 核心基因可结合 PPI 网络（STRING/Cytoscape）进一步筛选 hub genes。
