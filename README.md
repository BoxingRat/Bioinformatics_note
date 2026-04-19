# Bioinformatics_note

## GSE212192 差异分析（R / DESeq2）

脚本：`scripts/gse212192_deseq2_analysis.R`

### 功能
- 下载 GEO 数据集 `GSE212192`。
- 过滤基因：仅保留在超过一半样本中表达（count > 0）的基因。
- 使用 `DESeq2` 进行差异分析。
- 输出差异表达结果时提供 `gene_symbol` 字段（genesymbol 格式）。
- 绘制火山图并标注 top 20 差异基因。
- 进行 GO 与 KEGG 富集分析。
- 额外输出并可视化 `RIPK1-RIPK3-MLKL` 轴。

### 依赖包
`GEOquery`, `DESeq2`, `dplyr`, `tibble`, `stringr`, `ggplot2`, `ggrepel`, `clusterProfiler`, `org.Hs.eg.db`, `pheatmap`

### 运行方式（推荐先全自动）
```bash
Rscript scripts/gse212192_deseq2_analysis.R
```

脚本会自动：
- 猜测分组列 `group_col`；
- 在仅两组或关键词可识别的情况下自动识别 `control/treat`。
- 优先从 `GSEMatrix` 读取表达矩阵；若矩阵为空，会自动尝试从 GEO supplementary files 中解析 count matrix。
- 会自动清理样本名中的隐藏换行/空白（如 `WT1\\n`），减少样本对齐时意外丢样。
- 会自动清理分组标签中的误输入换行（例如手动输入 `WT1` 时误回车）。
- 会自动从 `characteristics_ch1*` 中提取形如 `group: WT` / `condition: disease` 的键值作为候选分组列，提高自动分组识别成功率。
- 读取 count matrix 时会识别“首列为空表头、基因ID在 rownames”的格式，避免把第一列样本（如 `WT1`）误当成基因列而丢失。

### 手动指定参数
```bash
Rscript scripts/gse212192_deseq2_analysis.R [group_col] [control_label] [treat_label] [outdir]
```

示例：
```bash
Rscript scripts/gse212192_deseq2_analysis.R source_name_ch1 control disease results/GSE212192
```

如自动识别失败，脚本会提示可用列名和分组水平，请按提示手动补充参数。

### 主要输出文件
- `run_config.txt`
- `DEG_results_genesymbol.csv`
- `volcano_top20.png`
- `GO_enrichment.csv`
- `GO_dotplot.png`
- `KEGG_enrichment.csv`
- `KEGG_dotplot.png`
- `RIPK1_RIPK3_MLKL_DEG.csv`
- `RIPK1_RIPK3_MLKL_heatmap.png`
