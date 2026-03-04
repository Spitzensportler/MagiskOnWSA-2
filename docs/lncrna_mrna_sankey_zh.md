# TCGA lncRNA-mRNA 相关性分析与桑基图流程（仅需两个表达矩阵）

你现在只需要提供 **2 个文件**：
- `lncRNA 表达矩阵`
- `mRNA 表达矩阵`

不再需要单独提供基因列表文件。

## 1. 输入文件格式

两个文件格式一致：
- 第一列列名必须是 `gene`
- 后续每一列是样本
- 值建议是标准化后的表达值（如 `log2(TPM+1)`）

例如：

`lncrna_expression.tsv`

```text
gene	Sample1	Sample2	Sample3
LINC00001	2.1	1.8	2.5
LINC00002	3.0	2.7	3.2
```

`mrna_expression.tsv`

```text
gene	Sample1	Sample2	Sample3
TP53	6.2	5.7	6.8
EGFR	4.1	4.4	4.0
```

> 注意：两个矩阵需要有共同样本名（至少 3 个）。

## 2. 参数提前赋值（推荐）

脚本顶部可直接修改默认参数：
- `default_lnc_expr_file`
- `default_mrna_expr_file`
- `default_out_prefix`
- `default_cor_cutoff`
- `default_fdr_cutoff`
- `pdf_width` / `pdf_height`
- `max_links_for_plot`

改好后可直接运行：

```bash
Rscript scripts/lncrna_mrna_correlation_sankey.R
```

## 3. 命令行运行（覆盖默认参数）

```bash
Rscript scripts/lncrna_mrna_correlation_sankey.R \
  lncrna_expression.tsv \
  mrna_expression.tsv \
  result/tcga \
  0.3 \
  0.05
```

参数含义：
- 第 1 个参数：lncRNA 表达矩阵
- 第 2 个参数：mRNA 表达矩阵
- 第 3 个参数：输出前缀
- 第 4 个参数：相关系数阈值（默认 0.3）
- 第 5 个参数：FDR 阈值（默认 0.05）

## 4. 输出结果

脚本会输出 3 个文件：
- `result/tcga_all_correlations.tsv`：全部 lncRNA-mRNA 配对的相关系数、p 值、FDR
- `result/tcga_significant_edges.tsv`：满足阈值过滤后的显著配对
- `result/tcga_sankey.pdf`：静态桑基图（按 lncRNA 来源进行丰富配色）

## 5. 方法说明

- 相关性方法：Spearman
- 多重检验：Benjamini-Hochberg（BH）校正
- 仅使用两类矩阵的共同样本进行计算

## 6. 依赖包

```r
install.packages(c("data.table", "Hmisc", "ggplot2", "ggalluvial"))
```
