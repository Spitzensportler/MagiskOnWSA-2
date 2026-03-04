# TCGA lncRNA-mRNA 相关性分析与桑基图流程

这个流程适合你现在的场景：
- 一组 lncRNA（数量较多）
- 一组 mRNA（筛选后数量较少）
- 想做两组之间的相关性分析并绘制桑基图

> 当前版本输出的是**静态 PDF 桑基图**（非交互式）。

## 1. 输入文件准备

需要 3 个文件：

1) 表达矩阵：`expression.tsv`
- 第一列列名必须是 `gene`
- 后续每一列是样本
- 值建议是标准化后的表达值（如 log2(TPM+1)）

示例：

```text
gene	Sample1	Sample2	Sample3
LINC00001	2.1	1.8	2.5
TP53	6.2	5.7	6.8
```

2) lncRNA 列表：`lncrna_list.txt`
- 每行一个基因符号

3) mRNA 列表：`mrna_list.txt`
- 每行一个基因符号

## 2. 参数提前赋值（推荐）

脚本顶部提供了可直接修改的默认参数：

- `default_expr_file`
- `default_lnc_file`
- `default_mrna_file`
- `default_out_prefix`
- `default_cor_cutoff`
- `default_fdr_cutoff`
- `pdf_width` / `pdf_height`（PDF 尺寸）
- `max_links_for_plot`（最多绘制多少条显著连线，避免图太乱）

你可以先在脚本里改好这些参数，然后直接运行：

```bash
Rscript scripts/lncrna_mrna_correlation_sankey.R
```

## 3. 命令行运行（可覆盖默认参数）

```bash
Rscript scripts/lncrna_mrna_correlation_sankey.R \
  expression.tsv \
  lncrna_list.txt \
  mrna_list.txt \
  result/tcga \
  0.3 \
  0.05
```

参数含义：
- 第 1~4 个参数：表达矩阵、lncRNA 列表、mRNA 列表、输出前缀
- 第 5 个参数：相关系数阈值（默认 0.3）
- 第 6 个参数：FDR 阈值（默认 0.05）

## 4. 输出结果

脚本会输出 3 个文件：

- `result/tcga_all_correlations.tsv`：全部 lncRNA-mRNA 配对的相关系数、p 值、FDR
- `result/tcga_significant_edges.tsv`：满足阈值过滤后的显著配对
- `result/tcga_sankey.pdf`：静态桑基图（颜色丰富，按 lncRNA 来源着色）

## 5. 统计建议

- 相关性方法：Spearman（脚本中已固定）
- 多重检验：Benjamini-Hochberg（BH）校正
- 推荐在同一癌种、同一批样本中进行分析，避免批次效应和混杂因素影响
- 若显著边过多导致图拥挤，可调小 `max_links_for_plot`

## 6. 依赖包

脚本依赖以下 R 包：

- `data.table`
- `Hmisc`
- `ggplot2`
- `ggalluvial`

如果缺包，请先安装：

```r
install.packages(c("data.table", "Hmisc", "ggplot2", "ggalluvial"))
```
