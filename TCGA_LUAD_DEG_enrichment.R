#!/usr/bin/env Rscript

# TCGA-LUAD 转录组分析流程
# 目标：
# 1) 下载并整理 TCGA LUAD RNA-seq 数据
# 2) 预处理（缺失值处理 + TMM标准化）
# 3) 差异表达分析（|log2FC| > 2, P < 0.05）
# 4) 核心差异基因功能富集分析（GO/KEGG）

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(dplyr)
  library(stringr)
  library(edgeR)
  library(limma)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(ggrepel)
  library(tibble)
})

set.seed(123)

# ========== 0. 参数设置 ==========
out_dir <- "results_luad"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 下载参数（网络不稳定时推荐使用 gdc-client）
# 可选: "api" 或 "client"
download_method <- "client"
files_per_chunk <- 20

# ========== 1. 下载 TCGA-LUAD STAR Counts ==========
query <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# ========== 1.1 下载数据（支持断点续传/失败重试） ==========
if (download_method == "client" && !nzchar(Sys.which("gdc-client"))) {
  stop(
    paste0(
      "当前设置 download_method='client'，但系统未检测到 gdc-client。\n",
      "请先安装并加入 PATH，或将 download_method 改为 'api'。"
    )
  )
}

message("开始下载，方法: ", download_method)
GDCdownload(
  query,
  method = download_method,
  files.per.chunk = files_per_chunk
)
se <- GDCprepare(query)

# 保存原始SummarizedExperiment对象
saveRDS(se, file.path(out_dir, "TCGA_LUAD_SE.rds"))

# ========== 2. 构建表达矩阵与分组 ==========
# STAR - Counts 数据可能包含多个assay（unstranded/stranded等），优先使用unstranded
assay_names <- assayNames(se)
if ("unstranded" %in% assay_names) {
  expr_raw <- assay(se, "unstranded")
} else {
  expr_raw <- assay(se)
  message("未找到 'unstranded' assay，已使用默认 assay: ", assay_names[1])
}
meta <- as.data.frame(colData(se))

# 使用样本条形码第14-15位判定样本类型：
# 01-09 为肿瘤样本；10-19 为正常样本
sample_type_code <- substr(colnames(expr_raw), 14, 15)
group <- ifelse(as.integer(sample_type_code) < 10, "Tumor", "Normal")

# 只保留原发肿瘤与实体正常样本
keep_samples <- group %in% c("Tumor", "Normal")
expr_raw <- expr_raw[, keep_samples, drop = FALSE]
group <- factor(group[keep_samples], levels = c("Normal", "Tumor"))
meta <- meta[keep_samples, , drop = FALSE]

# ========== 3. 缺失值处理 ==========
# RNA-seq count理论上应为整数，若出现NA，这里用基因在样本中的中位数填补
if (anyNA(expr_raw)) {
  message("检测到缺失值，进行中位数填补...")
  expr_raw <- apply(expr_raw, 1, function(x) {
    x[is.na(x)] <- median(x, na.rm = TRUE)
    x
  })
  expr_raw <- t(expr_raw)
}

# 确保行为基因、列为样本
expr_raw <- as.matrix(expr_raw)

# ========== 4. 过滤低表达基因 + TMM标准化 ==========
dge <- DGEList(counts = expr_raw, group = group)
keep_gene <- filterByExpr(dge, group = group)
dge <- dge[keep_gene, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")

# 保存标准化后的logCPM
logcpm <- cpm(dge, log = TRUE, prior.count = 1)
write.csv(logcpm, file.path(out_dir, "LUAD_logCPM_normalized.csv"), quote = FALSE)

# ========== 5. 差异表达分析（limma-voom） ==========
design <- model.matrix(~ group)
voom_obj <- voom(dge, design, plot = FALSE)
fit <- lmFit(voom_obj, design)
fit <- eBayes(fit)

# groupTumor系数代表 Tumor vs Normal
deg <- topTable(
  fit,
  coef = "groupTumor",
  number = Inf,
  sort.by = "P"
)

deg <- deg %>%
  rownames_to_column(var = "Gene") %>%
  mutate(
    regulation = case_when(
      logFC > 2 & P.Value < 0.05 ~ "Up",
      logFC < -2 & P.Value < 0.05 ~ "Down",
      TRUE ~ "NS"
    )
  )

# 差异结果补充 gene symbol（若可映射）
deg$Gene_clean <- str_replace(deg$Gene, "\\..*$", "")
deg_symbol_map <- bitr(
  unique(deg$Gene_clean),
  fromType = "ENSEMBL",
  toType = c("ENSEMBL", "SYMBOL"),
  OrgDb = org.Hs.eg.db
)

deg <- deg %>%
  left_join(deg_symbol_map, by = c("Gene_clean" = "ENSEMBL")) %>%
  mutate(GeneSymbol = ifelse(is.na(SYMBOL) | SYMBOL == "", Gene_clean, SYMBOL)) %>%
  select(-SYMBOL)

write.csv(deg, file.path(out_dir, "LUAD_DEG_all.csv"), row.names = FALSE)

deg_sig <- deg %>%
  filter(abs(logFC) > 2, P.Value < 0.05)
write.csv(deg_sig, file.path(out_dir, "LUAD_DEG_sig_logFC2_P0.05.csv"), row.names = FALSE)
write.csv(
  deg_sig %>% select(GeneSymbol, everything()),
  file.path(out_dir, "LUAD_DEG_sig_logFC2_P0.05_symbol.csv"),
  row.names = FALSE
)

message("显著差异基因数：", nrow(deg_sig))

# ========== 6. 核心差异基因筛选 ==========
# 这里定义“核心差异基因”为：显著DEG中按|logFC|排序前200（不足200则全部）
core_deg <- deg_sig %>%
  arrange(desc(abs(logFC))) %>%
  slice_head(n = 200)

write.csv(core_deg, file.path(out_dir, "LUAD_core_DEG_top200.csv"), row.names = FALSE)

# ========== 7. ID转换（ENSEMBL -> ENTREZID） ==========
id_map <- bitr(
  core_deg$Gene_clean,
  fromType = "ENSEMBL",
  toType = c("ENTREZID", "SYMBOL"),
  OrgDb = org.Hs.eg.db
)

entrez_ids <- unique(id_map$ENTREZID)

if (length(entrez_ids) < 10) {
  warning("可用于富集分析的ENTREZID数量过少，结果可能不稳定。")
}

# ========== 8. GO富集分析 ==========
ego <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
  go_res <- as.data.frame(ego)
  write.csv(go_res, file.path(out_dir, "LUAD_core_DEG_GO_enrichment.csv"), row.names = FALSE)

  p_go <- dotplot(ego, showCategory = 20, split = "ONTOLOGY") +
    facet_grid(ONTOLOGY ~ ., scales = "free") +
    ggtitle("GO enrichment of core DEGs (LUAD)")
  ggsave(file.path(out_dir, "LUAD_core_DEG_GO_dotplot.png"), p_go, width = 10, height = 8)
}

# ========== 9. KEGG富集分析 ==========
ekegg <- enrichKEGG(
  gene = entrez_ids,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2
)

if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
  kegg_res <- as.data.frame(ekegg)
  write.csv(kegg_res, file.path(out_dir, "LUAD_core_DEG_KEGG_enrichment.csv"), row.names = FALSE)

  p_kegg <- dotplot(ekegg, showCategory = 20) +
    ggtitle("KEGG enrichment of core DEGs (LUAD)")
  ggsave(file.path(out_dir, "LUAD_core_DEG_KEGG_dotplot.png"), p_kegg, width = 10, height = 6)
}

# ========== 10. 火山图 ==========
volcano_df <- deg %>%
  mutate(
    negLog10P = -log10(P.Value),
    group_plot = case_when(
      logFC > 2 & P.Value < 0.05 ~ "Up",
      logFC < -2 & P.Value < 0.05 ~ "Down",
      TRUE ~ "NS"
    )
  )

p_vol <- ggplot(volcano_df, aes(x = logFC, y = negLog10P, color = group_plot)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("Up" = "#D62728", "Down" = "#1F77B4", "NS" = "grey70")) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_bw() +
  labs(title = "TCGA-LUAD Differential Expression Volcano Plot",
       x = "log2 Fold Change",
       y = "-log10(P value)")

# 标注 top20 差异基因（按 P 值从小到大）
top20_label <- volcano_df %>%
  filter(abs(logFC) > 2, P.Value < 0.05) %>%
  arrange(P.Value) %>%
  slice_head(n = 20)

if (nrow(top20_label) > 0) {
  p_vol <- p_vol +
    geom_text_repel(
      data = top20_label,
      aes(label = GeneSymbol),
      size = 3,
      box.padding = 0.35,
      point.padding = 0.2,
      max.overlaps = Inf,
      show.legend = FALSE
    )
}

ggsave(file.path(out_dir, "LUAD_DEG_volcano.png"), p_vol, width = 8, height = 6)

message("分析完成。结果目录：", normalizePath(out_dir))
