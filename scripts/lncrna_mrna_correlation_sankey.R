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
default_expr_file <- "expression.tsv"
default_lnc_file <- "lncrna_list.txt"
default_mrna_file <- "mrna_list.txt"
default_out_prefix <- "result/tcga"
default_cor_cutoff <- 0.3
default_fdr_cutoff <- 0.05

# PDF 输出参数
# max_links_for_plot: 为避免连线过多导致图像不可读，默认仅绘制相关性绝对值最高的前 N 条
# 若想绘制全部显著边，可设为 Inf
pdf_width <- 14
pdf_height <- 10
max_links_for_plot <- 200

# -----------------------------
# 命令行参数解析
# 若不传参，则使用上方预设参数
# 若传参，则按顺序覆盖预设参数
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

expr_file <- ifelse(length(args) >= 1, args[1], default_expr_file)
lnc_file <- ifelse(length(args) >= 2, args[2], default_lnc_file)
mrna_file <- ifelse(length(args) >= 3, args[3], default_mrna_file)
out_prefix <- ifelse(length(args) >= 4, args[4], default_out_prefix)
cor_cutoff <- ifelse(length(args) >= 5, as.numeric(args[5]), default_cor_cutoff)
fdr_cutoff <- ifelse(length(args) >= 6, as.numeric(args[6]), default_fdr_cutoff)

if (is.na(cor_cutoff) || is.na(fdr_cutoff)) {
  stop("cor_cutoff and fdr_cutoff must be numeric.")
}

if (length(args) > 0 && length(args) < 4) {
  warning(
    paste0(
      "Detected partial arguments. Recommended full usage: Rscript scripts/lncrna_mrna_correlation_sankey.R ",
      "<expression.tsv> <lncrna_list.txt> <mrna_list.txt> <output_prefix> [cor_cutoff=0.3] [fdr_cutoff=0.05]"
    )
  )
}

# 创建输出目录（如果不存在）
out_dir <- dirname(out_prefix)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# -----------------------------
# 读取输入数据
# -----------------------------
expr <- fread(expr_file)
if (!"gene" %in% colnames(expr)) {
  stop("Expression file must contain a first column named 'gene'.")
}

lnc_genes <- fread(lnc_file, header = FALSE)$V1
mrna_genes <- fread(mrna_file, header = FALSE)$V1

# 仅保留两组基因的并集，减少后续计算量
expr_sub <- expr[gene %in% unique(c(lnc_genes, mrna_genes))]
if (nrow(expr_sub) == 0) {
  stop("No overlapping genes found between expression matrix and input gene lists.")
}

# 转为矩阵，行为基因，列为样本
expr_mat <- as.matrix(expr_sub[, -1, with = FALSE])
rownames(expr_mat) <- expr_sub$gene

# 分别筛选在表达矩阵中真实存在的 lncRNA 和 mRNA
lnc_keep <- intersect(lnc_genes, rownames(expr_mat))
mrna_keep <- intersect(mrna_genes, rownames(expr_mat))

if (length(lnc_keep) < 2 || length(mrna_keep) < 2) {
  stop("Need at least 2 lncRNAs and 2 mRNAs after overlap filtering.")
}

# -----------------------------
# 相关性分析（Spearman）
# -----------------------------
lnc_mat <- t(expr_mat[lnc_keep, , drop = FALSE])
mrna_mat <- t(expr_mat[mrna_keep, , drop = FALSE])
all_mat <- cbind(lnc_mat, mrna_mat)

corr <- rcorr(all_mat, type = "spearman")
r_mat <- corr$r[colnames(lnc_mat), colnames(mrna_mat), drop = FALSE]
p_mat <- corr$P[colnames(lnc_mat), colnames(mrna_mat), drop = FALSE]

# 整理为边表，并做 BH 多重校正
edges <- as.data.table(as.table(r_mat))
setnames(edges, c("lncRNA", "mRNA", "rho"))
edges[, pvalue := as.vector(p_mat)]
edges[, fdr := p.adjust(pvalue, method = "BH")]

# 按阈值筛选显著边
edges_sig <- edges[abs(rho) >= cor_cutoff & fdr <= fdr_cutoff]
setorder(edges_sig, -abs(rho))

# 输出表格结果
fwrite(edges, paste0(out_prefix, "_all_correlations.tsv"), sep = "\t")
fwrite(edges_sig, paste0(out_prefix, "_significant_edges.tsv"), sep = "\t")

if (nrow(edges_sig) == 0) {
  message("No significant lncRNA-mRNA pairs under current thresholds. Only table outputs were generated.")
  quit(save = "no")
}

# -----------------------------
# 构建静态桑基图（PDF）
# -----------------------------
plot_edges <- copy(edges_sig)
if (is.finite(max_links_for_plot) && nrow(plot_edges) > max_links_for_plot) {
  plot_edges <- plot_edges[1:max_links_for_plot]
}

# 将每条边转换为 alluvial 的长表结构
alluv_long <- rbind(
  data.table(
    alluvium = seq_len(nrow(plot_edges)),
    axis = "lncRNA",
    stratum = plot_edges$lncRNA,
    weight = abs(plot_edges$rho),
    sign = ifelse(plot_edges$rho >= 0, "Positive", "Negative"),
    source_lnc = plot_edges$lncRNA
  ),
  data.table(
    alluvium = seq_len(nrow(plot_edges)),
    axis = "mRNA",
    stratum = plot_edges$mRNA,
    weight = abs(plot_edges$rho),
    sign = ifelse(plot_edges$rho >= 0, "Positive", "Negative"),
    source_lnc = plot_edges$lncRNA
  )
)

alluv_long[, axis := factor(axis, levels = c("lncRNA", "mRNA"))]

# 丰富配色：按 lncRNA 来源为每条流分配颜色（颜色数量随 lncRNA 数量扩展）
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
      "/", nrow(edges_sig)
    ),
    x = NULL,
    y = "Flow weight (|rho|)",
    fill = "Source lncRNA"
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

message("Done. Outputs:")
message("- ", paste0(out_prefix, "_all_correlations.tsv"))
message("- ", paste0(out_prefix, "_significant_edges.tsv"))
message("- ", pdf_file)
