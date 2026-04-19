#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GEOquery)
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(pheatmap)
})

print_usage <- function() {
  cat(
    "\n用法（参数可选，推荐先自动识别）：\n",
    "  Rscript scripts/gse212192_deseq2_analysis.R [group_col] [control_label] [treat_label] [outdir]\n\n",
    "示例1（全自动）：\n",
    "  Rscript scripts/gse212192_deseq2_analysis.R\n\n",
    "示例2（手动指定分组列）：\n",
    "  Rscript scripts/gse212192_deseq2_analysis.R source_name_ch1\n\n",
    "示例3（全手动）：\n",
    "  Rscript scripts/gse212192_deseq2_analysis.R source_name_ch1 control disease results/GSE212192\n\n",
    sep = ""
  )
}

pick_group_col <- function(pheno) {
  candidates <- colnames(pheno)
  score_col <- function(colname) {
    v <- as.character(pheno[[colname]])
    v <- str_trim(v)
    v <- v[!is.na(v) & v != ""]
    if (length(v) == 0) return(-Inf)

    n_unique <- length(unique(v))
    # 排除几乎每个样本都不同的列
    if (n_unique > length(v) * 0.8) return(-Inf)
    if (n_unique < 2 || n_unique > 6) return(-Inf)

    name_bonus <- ifelse(str_detect(tolower(colname), "group|condition|source|title|characteristics|char_"), 2, 0)
    binary_bonus <- ifelse(n_unique == 2, 3, 0)
    return(name_bonus + binary_bonus)
  }

  scores <- sapply(candidates, score_col)
  best <- names(which.max(scores))
  if (length(best) == 0 || is.infinite(scores[[best]])) {
    return(NA_character_)
  }
  best
}

extract_group_from_characteristics <- function(pheno) {
  char_cols <- grep("^characteristics_ch1", colnames(pheno), value = TRUE)
  if (length(char_cols) == 0) return(pheno)

  for (cc in char_cols) {
    raw <- as.character(pheno[[cc]])
    raw <- normalize_group_label(raw)
    has_kv <- grepl(":", raw)
    if (!any(has_kv)) next

    key <- ifelse(has_kv, tolower(str_trim(sub(":.*$", "", raw))), NA_character_)
    val <- ifelse(has_kv, str_trim(sub("^[^:]+:\\s*", "", raw)), raw)
    key_levels <- unique(na.omit(key))
    if (length(key_levels) != 1) next

    new_col <- paste0("char_", gsub("[^a-z0-9]+", "_", key_levels[1]))
    pheno[[new_col]] <- val
  }
  pheno
}

pick_labels <- function(levels_vec) {
  lv <- as.character(levels_vec)
  low <- tolower(lv)

  ctrl_idx <- which(str_detect(low, "control|ctrl|normal|healthy|wt|mock|vehicle|untreated"))
  trt_idx <- which(str_detect(low, "disease|case|model|treated|stim|infection|ko|kd|mut|injury|sepsis"))

  if (length(ctrl_idx) >= 1 && length(trt_idx) >= 1 && ctrl_idx[1] != trt_idx[1]) {
    return(list(control = lv[ctrl_idx[1]], treat = lv[trt_idx[1]], method = "keyword"))
  }

  if (length(lv) == 2) {
    if (length(ctrl_idx) == 1) {
      control <- lv[ctrl_idx]
      treat <- setdiff(lv, control)[1]
      return(list(control = control, treat = treat, method = "2level+ctrl_keyword"))
    }
    # 两组但无关键词，按字母序稳定选择
    lv_sorted <- sort(lv)
    return(list(control = lv_sorted[1], treat = lv_sorted[2], method = "2level_alphabetical"))
  }

  return(NULL)
}

pick_non_empty_eset <- function(gse_list) {
  if (length(gse_list) == 0) return(NULL)
  dims <- sapply(gse_list, function(es) nrow(exprs(es)) * ncol(exprs(es)))
  if (all(dims == 0)) return(NULL)
  gse_list[[which.max(dims)]]
}

read_count_matrix_from_file <- function(fp) {
  ext <- tools::file_ext(fp)
  sep <- ifelse(ext %in% c("csv"), ",", "\t")
  df <- tryCatch(
    utils::read.table(fp, header = TRUE, sep = sep, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = ""),
    error = function(e) NULL
  )
  if (is.null(df) || nrow(df) < 10 || ncol(df) < 2) return(NULL)

  default_rn <- as.character(seq_len(nrow(df)))
  has_nondefault_rownames <- !is.null(rownames(df)) && !all(rownames(df) == default_rn)

  # 情况A：文件第一列已被 read.table 当作 rownames（常见于首列表头为空的 count matrix）
  if (has_nondefault_rownames) {
    mat <- as.matrix(df)
    suppressWarnings(storage.mode(mat) <- "numeric")
    if (sum(!is.na(mat)) == 0) return(NULL)
    rownames(mat) <- make.unique(as.character(rownames(df)))
    return(mat)
  }

  # 情况B：第一列是明确的 gene id 列
  first_col <- df[[1]]
  first_col_num_ratio <- mean(!is.na(suppressWarnings(as.numeric(first_col))))
  if (first_col_num_ratio < 0.5) {
    mat <- as.matrix(df[, -1, drop = FALSE])
    suppressWarnings(storage.mode(mat) <- "numeric")
    if (sum(!is.na(mat)) == 0) return(NULL)
    rownames(mat) <- make.unique(as.character(first_col))
    return(mat)
  }

  # 情况C：第一列也是数值，无法判定 gene id，保留全部列并生成占位基因名（避免丢失第一样本列）
  mat <- as.matrix(df)
  suppressWarnings(storage.mode(mat) <- "numeric")
  if (sum(!is.na(mat)) == 0) return(NULL)
  rownames(mat) <- paste0("gene_", seq_len(nrow(mat)))
  mat
}

normalize_sample_id <- function(x) {
  x <- as.character(x)
  x <- gsub("\ufeff", "", x, fixed = TRUE)   # 去 BOM
  x <- gsub("[\r\n\t]", "", x)               # 去隐藏换行/制表
  x <- str_trim(x)                           # 去首尾空格
  x
}

normalize_group_label <- function(x) {
  x <- as.character(x)
  x <- gsub("[\r\n]", "", x)
  x <- str_trim(x)
  x
}

load_expr_with_fallback <- function(gse_id, outdir) {
  gse_list <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = TRUE)
  eset <- pick_non_empty_eset(gse_list)

  if (!is.null(eset)) {
    expr <- exprs(eset)
    pheno <- pData(eset)
    expr_ids <- normalize_sample_id(colnames(expr))
    pheno_ids <- normalize_sample_id(rownames(pheno))
    # 尽量按同一套 sample id 统一，避免 make.unique 分别处理导致错位丢样
    if (length(expr_ids) == length(pheno_ids)) {
      unified <- ifelse(expr_ids != "", expr_ids, pheno_ids)
      unified[unified == ""] <- paste0("Sample", seq_along(unified))[unified == ""]
      unified <- make.unique(unified)
      colnames(expr) <- unified
      rownames(pheno) <- unified
    } else {
      colnames(expr) <- make.unique(expr_ids)
      rownames(pheno) <- make.unique(pheno_ids)
    }
    return(list(expr = expr, pheno = pheno, source = "GSEMatrix"))
  }

  message("GSEMatrix 中未发现可用表达矩阵，尝试从 supplementary files 读取 count matrix...")
  supp_dir <- file.path(outdir, "supp_files")
  dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
  GEOquery::getGEOSuppFiles(gse_id, baseDir = supp_dir, makeDirectory = TRUE)
  gse_supp_dir <- file.path(supp_dir, gse_id)

  tar_files <- list.files(gse_supp_dir, pattern = "\\.tar$", full.names = TRUE)
  for (tf in tar_files) {
    untar(tf, exdir = gse_supp_dir)
  }

  files <- list.files(gse_supp_dir, recursive = TRUE, full.names = TRUE)
  files <- files[grepl("\\.(txt|tsv|csv|gz)$", files, ignore.case = TRUE)]
  files <- files[!grepl("series_matrix", basename(files), ignore.case = TRUE)]
  if (length(files) == 0) {
    stop("未找到可读取的 supplementary 表达矩阵文件。")
  }

  # 优先尝试名称更像 count/matrix 的文件
  priority <- grepl("count|raw|matrix|read", basename(files), ignore.case = TRUE)
  files <- c(files[priority], files[!priority])

  for (fp in files) {
    mat <- read_count_matrix_from_file(fp)
    if (!is.null(mat)) {
      sample_names <- colnames(mat)
      if (is.null(sample_names) || length(sample_names) == 0) {
        sample_names <- paste0("Sample", seq_len(ncol(mat)))
      }
      sample_names <- make.unique(normalize_sample_id(sample_names))
      colnames(mat) <- sample_names

      pheno <- if (length(gse_list) >= 1) pData(gse_list[[1]]) else data.frame(row.names = colnames(mat))
      if (nrow(pheno) > 0) {
        rownames(pheno) <- make.unique(normalize_sample_id(rownames(pheno)))
      }
      if (!all(colnames(mat) %in% rownames(pheno))) {
        inter <- intersect(colnames(mat), rownames(pheno))
        if (length(inter) >= 2) {
          mat <- mat[, inter, drop = FALSE]
          pheno <- pheno[inter, , drop = FALSE]
        } else {
          pheno <- data.frame(group = rep("unknown", length(sample_names)),
                              row.names = sample_names,
                              stringsAsFactors = FALSE)
        }
      } else {
        pheno <- pheno[colnames(mat), , drop = FALSE]
      }
      return(list(expr = mat, pheno = pheno, source = paste0("supp:", basename(fp))))
    }
  }

  stop("supplementary files 中未能识别有效 count matrix。")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && args[[1]] %in% c("-h", "--help")) {
  print_usage()
  quit(save = "no", status = 0)
}

group_col <- ifelse(length(args) >= 1, args[[1]], NA_character_)
control_label <- ifelse(length(args) >= 2, args[[2]], NA_character_)
treat_label <- ifelse(length(args) >= 3, args[[3]], NA_character_)
outdir <- ifelse(length(args) >= 4, args[[4]], "results/GSE212192")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

message("[1/8] 下载并读取 GSE212192...")
dat <- load_expr_with_fallback("GSE212192", outdir)
expr <- dat$expr
pheno <- dat$pheno
pheno <- extract_group_from_characteristics(pheno)
message(sprintf("表达矩阵来源: %s; 维度: %d genes x %d samples", dat$source, nrow(expr), ncol(expr)))

if (is.na(group_col) || group_col == "") {
  group_col <- pick_group_col(pheno)
  if (is.na(group_col)) {
    print_usage()
    stop(sprintf(
      "无法自动识别分组列。请手动指定 group_col。可选列：%s",
      paste(colnames(pheno), collapse = ", ")
    ))
  }
  message(sprintf("自动识别 group_col = '%s'", group_col))
}

if (!group_col %in% colnames(pheno)) {
  print_usage()
  stop(sprintf("group_col '%s' 不在样本信息中。可选列：%s", group_col, paste(colnames(pheno), collapse = ", ")))
}

message("[2/8] 构建分组信息...")
pheno$group <- normalize_group_label(pheno[[group_col]])
pheno <- pheno[!is.na(pheno$group) & pheno$group != "", , drop = FALSE]

all_levels <- unique(pheno$group)
if (length(all_levels) < 2) {
  stop(sprintf("分组列 '%s' 只有 %d 个水平，无法做差异分析。", group_col, length(all_levels)))
}

if (is.na(control_label) || control_label == "" || is.na(treat_label) || treat_label == "") {
  auto <- pick_labels(all_levels)
  if (is.null(auto)) {
    print_usage()
    stop(sprintf(
      "检测到分组水平为：%s。请手动传入 control_label 与 treat_label。",
      paste(sort(all_levels), collapse = ", ")
    ))
  }
  control_label <- auto$control
  treat_label <- auto$treat
  message(sprintf("自动识别分组: control='%s', treat='%s' (%s)", control_label, treat_label, auto$method))
}

control_label <- normalize_group_label(control_label)
treat_label <- normalize_group_label(treat_label)

if (!(control_label %in% all_levels) || !(treat_label %in% all_levels)) {
  stop(sprintf("control_label (%s) 或 treat_label (%s) 不在 group 列水平中：%s",
               control_label, treat_label, paste(sort(all_levels), collapse = ", ")))
}

keep_samples <- pheno$group %in% c(control_label, treat_label)
pheno <- pheno[keep_samples, , drop = FALSE]
expr <- expr[, rownames(pheno), drop = FALSE]
pheno$group <- factor(pheno$group, levels = c(control_label, treat_label))

message(sprintf("使用样本数: %d（%s=%d, %s=%d）",
                nrow(pheno),
                control_label, sum(pheno$group == control_label),
                treat_label, sum(pheno$group == treat_label)))

writeLines(
  c(
    sprintf("group_col=%s", group_col),
    sprintf("control_label=%s", control_label),
    sprintf("treat_label=%s", treat_label),
    sprintf("n_samples=%d", nrow(pheno))
  ),
  con = file.path(outdir, "run_config.txt")
)

message("[3/8] 过滤基因：保留在一半以上样本里表达（count > 0）的基因...")
expr <- round(expr)
expr[expr < 0] <- 0
keep_gene <- rowSums(expr > 0) > (ncol(expr) / 2)
expr_filt <- expr[keep_gene, , drop = FALSE]
message(sprintf("原始基因数: %d; 过滤后基因数: %d", nrow(expr), nrow(expr_filt)))

message("[4/8] DESeq2 差异分析...")
col_data <- pheno %>%
  rownames_to_column("sample") %>%
  dplyr::select(sample, group) %>%
  column_to_rownames("sample")

dds <- DESeqDataSetFromMatrix(countData = expr_filt, colData = col_data, design = ~ group)
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", treat_label, control_label))
res_df <- as.data.frame(res) %>% rownames_to_column("gene_id") %>% arrange(padj)

message("[5/8] 将差异基因名称统一为 genesymbol...")
id_type <- ifelse(all(str_detect(res_df$gene_id, "^ENSG")), "ENSEMBL",
                  ifelse(all(str_detect(res_df$gene_id, "^[0-9]+$")), "ENTREZID", "SYMBOL"))

map_df <- tryCatch({
  bitr(unique(res_df$gene_id), fromType = id_type, toType = c("SYMBOL", "ENTREZID"), OrgDb = org.Hs.eg.db)
}, error = function(e) {
  message("ID 映射失败，默认使用原始 gene_id 作为 SYMBOL。")
  data.frame()
})

if (nrow(map_df) > 0) {
  colnames(map_df)[1] <- "gene_id"
  res_df <- res_df %>%
    left_join(map_df, by = "gene_id") %>%
    mutate(gene_symbol = ifelse(is.na(SYMBOL), gene_id, SYMBOL))
} else {
  res_df <- res_df %>% mutate(gene_symbol = gene_id, ENTREZID = NA_character_)
}

res_df <- res_df %>% mutate(significance = case_when(
  !is.na(padj) & padj < 0.05 & log2FoldChange >= 1 ~ "Up",
  !is.na(padj) & padj < 0.05 & log2FoldChange <= -1 ~ "Down",
  TRUE ~ "NS"
))

write.csv(res_df, file.path(outdir, "DEG_results_genesymbol.csv"), row.names = FALSE)

message("[6/8] 绘制火山图（标注 top 20 差异基因）...")
volcano_df <- res_df %>% mutate(neg_log10_padj = -log10(ifelse(is.na(padj), 1, padj)))
top20 <- volcano_df %>% filter(!is.na(padj)) %>% arrange(padj) %>% slice_head(n = 20)

p_vol <- ggplot(volcano_df, aes(x = log2FoldChange, y = neg_log10_padj, color = significance)) +
  geom_point(alpha = 0.7, size = 1.8) +
  scale_color_manual(values = c(Up = "#D73027", Down = "#4575B4", NS = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = 2, color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey40") +
  ggrepel::geom_text_repel(data = top20, aes(label = gene_symbol), size = 3, max.overlaps = Inf, box.padding = 0.5) +
  theme_bw(base_size = 12) +
  labs(title = "GSE212192 Volcano Plot", x = "log2 Fold Change", y = "-log10 adjusted P", color = "Category")

ggsave(file.path(outdir, "volcano_top20.png"), p_vol, width = 8, height = 6, dpi = 300)

message("[7/8] GO / KEGG 富集分析...")
sig <- res_df %>% filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1)
entrez_sig <- unique(na.omit(sig$ENTREZID))

if (length(entrez_sig) >= 10) {
  ego <- enrichGO(gene = entrez_sig, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "ALL", pAdjustMethod = "BH", qvalueCutoff = 0.2, readable = TRUE)
  write.csv(as.data.frame(ego), file.path(outdir, "GO_enrichment.csv"), row.names = FALSE)
  p_go <- dotplot(ego, showCategory = 20) + ggtitle("GO Enrichment")
  ggsave(file.path(outdir, "GO_dotplot.png"), p_go, width = 10, height = 8, dpi = 300)

  ekegg <- enrichKEGG(gene = entrez_sig, organism = "hsa", pAdjustMethod = "BH", qvalueCutoff = 0.2)
  ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  write.csv(as.data.frame(ekegg), file.path(outdir, "KEGG_enrichment.csv"), row.names = FALSE)
  p_kegg <- dotplot(ekegg, showCategory = 20) + ggtitle("KEGG Enrichment")
  ggsave(file.path(outdir, "KEGG_dotplot.png"), p_kegg, width = 10, height = 8, dpi = 300)
} else {
  message("显著差异基因可映射的 ENTREZID 少于 10，跳过 GO/KEGG 富集。")
}

message("[8/8] 查看 RIPK1-RIPK3-MLKL 轴...")
target_axis <- c("RIPK1", "RIPK3", "MLKL")
axis_df <- res_df %>%
  filter(gene_symbol %in% target_axis) %>%
  dplyr::select(gene_id, gene_symbol, log2FoldChange, pvalue, padj, significance)
write.csv(axis_df, file.path(outdir, "RIPK1_RIPK3_MLKL_DEG.csv"), row.names = FALSE)

vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

if (nrow(map_df) > 0) {
  symbol_map <- res_df %>% dplyr::select(gene_id, gene_symbol) %>% distinct()
  vsd_df <- as.data.frame(vsd_mat) %>% rownames_to_column("gene_id") %>% left_join(symbol_map, by = "gene_id")
  axis_expr <- vsd_df %>% filter(gene_symbol %in% target_axis)
} else {
  axis_expr <- as.data.frame(vsd_mat) %>% rownames_to_column("gene_symbol") %>% filter(gene_symbol %in% target_axis)
}

if (nrow(axis_expr) > 0) {
  axis_mat <- axis_expr %>%
    dplyr::select(any_of(c("gene_symbol", "gene_id")), everything()) %>%
    dplyr::select(-any_of("gene_id")) %>%
    column_to_rownames("gene_symbol") %>%
    as.matrix()

  anno <- data.frame(group = col_data$group)
  rownames(anno) <- rownames(col_data)

  png(file.path(outdir, "RIPK1_RIPK3_MLKL_heatmap.png"), width = 2400, height = 1200, res = 300)
  pheatmap(axis_mat, annotation_col = anno, scale = "row", cluster_rows = FALSE, cluster_cols = TRUE, fontsize_row = 12,
           main = "RIPK1-RIPK3-MLKL Axis (VST)")
  dev.off()
} else {
  message("未在结果中找到 RIPK1/RIPK3/MLKL，已输出空结果表用于核对。")
}

message("分析完成。输出目录: ", outdir)
