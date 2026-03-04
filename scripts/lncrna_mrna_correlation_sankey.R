#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Hmisc)
  library(ggplot2)
  library(ggalluvial)
})

# -----------------------------
# 预设参数（可直接在这里修改）
# -----------------------------
default_lnc_expr_file <- "lncrna_expression.tsv"
default_mrna_expr_file <- "mrna_expression.tsv"
default_out_prefix <- "result/tcga"
default_cor_cutoff <- 0.3
default_fdr_cutoff <- 0.05

# PDF 输出参数
pdf_width <- 14
pdf_height <- 10
max_links_for_plot <- 200

# -----------------------------
# 命令行参数解析
# 用法：
# Rscript scripts/lncrna_mrna_correlation_sankey.R \
#   <lncrna_expression.tsv> <mrna_expression.tsv> <output_prefix> [cor_cutoff] [fdr_cutoff]
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

lnc_expr_file <- ifelse(length(args) >= 1, args[1], default_lnc_expr_file)
mrna_expr_file <- ifelse(length(args) >= 2, args[2], default_mrna_expr_file)
out_prefix <- ifelse(length(args) >= 3, args[3], default_out_prefix)
cor_cutoff <- ifelse(length(args) >= 4, as.numeric(args[4]), default_cor_cutoff)
fdr_cutoff <- ifelse(length(args) >= 5, as.numeric(args[5]), default_fdr_cutoff)

if (is.na(cor_cutoff) || is.na(fdr_cutoff)) {
  stop("cor_cutoff 和 fdr_cutoff 必须是数值。")
}

if (length(args) > 0 && length(args) < 3) {
  warning(
    paste0(
      "检测到仅传入了部分参数。推荐完整用法：Rscript scripts/lncrna_mrna_correlation_sankey.R ",
      "<lncrna_expression.tsv> <mrna_expression.tsv> <output_prefix> [cor_cutoff=0.3] [fdr_cutoff=0.05]"
    )
  )
}

# 创建输出目录（如果不存在）
out_dir <- dirname(out_prefix)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# -----------------------------
# 读取两个表达矩阵
# 要求：第一列列名为 gene，后续列为样本
# -----------------------------
lnc_expr <- fread(lnc_expr_file)
mrna_expr <- fread(mrna_expr_file)

if (!"gene" %in% colnames(lnc_expr)) {
  stop("lncRNA 表达矩阵必须包含名为 gene 的第一列。")
}
if (!"gene" %in% colnames(mrna_expr)) {
  stop("mRNA 表达矩阵必须包含名为 gene 的第一列。")
}

# 获取两个矩阵的共同样本
lnc_samples <- setdiff(colnames(lnc_expr), "gene")
mrna_samples <- setdiff(colnames(mrna_expr), "gene")
common_samples <- intersect(lnc_samples, mrna_samples)

if (length(common_samples) < 3) {
  stop("lncRNA 与 mRNA 表达矩阵共同样本数不足（至少需要 3 个）。")
}

# 去重基因，避免重复行导致相关性结果歧义
lnc_expr <- unique(lnc_expr, by = "gene")
mrna_expr <- unique(mrna_expr, by = "gene")

# 构建数值矩阵：行为样本，列为基因
lnc_mat <- as.matrix(lnc_expr[, ..common_samples])
colnames(lnc_mat) <- lnc_expr$gene
mrna_mat <- as.matrix(mrna_expr[, ..common_samples])
colnames(mrna_mat) <- mrna_expr$gene

# 强制转数值，非数值会变成 NA
mode(lnc_mat) <- "numeric"
mode(mrna_mat) <- "numeric"

if (ncol(lnc_mat) < 2 || ncol(mrna_mat) < 2) {
  stop("至少需要 2 个 lncRNA 基因和 2 个 mRNA 基因。")
}

all_mat <- cbind(lnc_mat, mrna_mat)

# -----------------------------
# 相关性分析（Spearman）
# -----------------------------
corr <- rcorr(all_mat, type = "spearman")
r_mat <- corr$r[colnames(lnc_mat), colnames(mrna_mat), drop = FALSE]
p_mat <- corr$P[colnames(lnc_mat), colnames(mrna_mat), drop = FALSE]

edges <- as.data.table(as.table(r_mat))
setnames(edges, c("lncRNA", "mRNA", "rho"))
edges[, pvalue := as.vector(p_mat)]
edges[, fdr := p.adjust(pvalue, method = "BH")]

# 去掉无法计算相关性的组合（通常由全 NA 或常量列引起）
edges <- edges[!is.na(rho) & !is.na(pvalue) & !is.na(fdr)]

edges_sig <- edges[abs(rho) >= cor_cutoff & fdr <= fdr_cutoff]
setorder(edges_sig, -abs(rho))

fwrite(edges, paste0(out_prefix, "_all_correlations.tsv"), sep = "\t")
fwrite(edges_sig, paste0(out_prefix, "_significant_edges.tsv"), sep = "\t")

if (nrow(edges_sig) == 0) {
  message("在当前阈值下未发现显著的 lncRNA-mRNA 配对，仅输出结果表格。")
  quit(save = "no")
}

# -----------------------------
# 构建静态桑基图（PDF）
# -----------------------------
plot_edges <- copy(edges_sig)
if (is.finite(max_links_for_plot) && nrow(plot_edges) > max_links_for_plot) {
  plot_edges <- plot_edges[1:max_links_for_plot]
}

alluv_long <- rbind(
  data.table(
    alluvium = seq_len(nrow(plot_edges)),
    axis = "lncRNA",
    stratum = plot_edges$lncRNA,
    weight = abs(plot_edges$rho),
    source_lnc = plot_edges$lncRNA
  ),
  data.table(
    alluvium = seq_len(nrow(plot_edges)),
    axis = "mRNA",
    stratum = plot_edges$mRNA,
    weight = abs(plot_edges$rho),
    source_lnc = plot_edges$lncRNA
  )
)

alluv_long[, axis := factor(axis, levels = c("lncRNA", "mRNA"))]

lnc_levels <- unique(plot_edges$lncRNA)
lnc_colors <- grDevices::hcl.colors(length(lnc_levels), palette = "Dynamic")
names(lnc_colors) <- lnc_levels

p <- ggplot(
  alluv_long,
  aes(x = axis, stratum = stratum, alluvium = alluvium, y = weight, fill = source_lnc)
) +
  geom_flow(alpha = 0.78, color = "gray35", width = 0.2) +
  geom_stratum(width = 0.35, color = "gray20", fill = "gray92") +
  scale_fill_manual(values = lnc_colors, guide = "none") +
  labs(
    title = "lncRNA-mRNA Correlation Sankey (Spearman)",
    subtitle = paste0(
      "|rho| >= ", cor_cutoff,
      ", FDR <= ", fdr_cutoff,
      "; plotted links: ", nrow(plot_edges),
      "/", nrow(edges_sig),
      "; common samples: ", length(common_samples)
    ),
    x = NULL,
    y = "Flow weight (|rho|)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(face = "bold")
  )

pdf_file <- paste0(out_prefix, "_sankey.pdf")
ggsave(pdf_file, p, width = pdf_width, height = pdf_height, device = "pdf")

message("分析完成，输出文件如下：")
message("- ", paste0(out_prefix, "_all_correlations.tsv"))
message("- ", paste0(out_prefix, "_significant_edges.tsv"))
message("- ", pdf_file)
