# ========================= 0. 安装和加载所需包 =================================
# 如果尚未安装，请取消注释并运行（首次使用时安装）
# install.packages(c("tidyverse", "data.table", "ggplot2"))
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("MCPcounter")
# pak::pkg_install("GfellerLab/EPIC")   # 使用 pak 安装 EPIC

library(MCPcounter)
library(EPIC)
library(data.table)
library(tidyverse)
library(ggplot2)

# 设置工作目录（请修改为你的实际路径）
setwd("E:/NETS/")

# ========================= 1. 读取Bulk表达矩阵 =================================
expr_file <- "bulk rna数据/原数据/1_CGGA_118_Filtered/CGGA.mRNAseq_693_filtered.txt"
expr_raw <- fread(expr_file, data.table = FALSE)
rownames(expr_raw) <- expr_raw[, 1]   # 第一列为基因名
expr_raw <- expr_raw[, -1]            # 移除基因名列

# 转换为数值矩阵（注意：fread 可能将部分列读为字符，需强制转换）
expr_mat <- as.matrix(expr_raw)
mode(expr_mat) <- "numeric"
# 转换时产生的 NA 警告是正常的（非数值列已被移除），这里用 0 填充 NA
expr_mat[is.na(expr_mat)] <- 0

# 去除全零基因（可选）
expr_mat <- expr_mat[rowSums(expr_mat) > 0, ]
cat("表达矩阵维度 (基因 × 样本):", dim(expr_mat), "\n")

# ========================= 2. 基因名标准化 ====================================
# 去除基因名中的点号及后缀（例如 "TP53.1" -> "TP53"）
rownames(expr_mat) <- gsub("\\..*", "", rownames(expr_mat))

# 重要：MCPcounter 和 EPIC 需要标准基因符号（HGNC）。如果你的数据是 Ensembl ID，
# 请使用 biomaRt 转换为基因符号，否则后续会匹配失败。这里假设已经是基因符号。

# ========================= 3. 运行 MCP-counter ================================
cat("运行 MCP-counter...\n")
# MCPcounter.estimate 会自动匹配内置基因列表，不需要手动过滤
# 注意：表达矩阵需要是 基因 × 样本，且为原始计数或 TPM/FPKM（不要 log 转换）
mcp_res <- tryCatch(
  MCPcounter.estimate(
    expression = expr_mat,
    featuresType = "HUGO_symbols"   # 指定基因名类型
  ),
  error = function(e) {
    message("MCP-counter 运行出错: ", e$message)
    return(NULL)
  }
)

if (!is.null(mcp_res)) {
  mcp_scores <- as.data.frame(t(mcp_res))
  colnames(mcp_scores) <- paste0("MCP_", colnames(mcp_scores))
  cat("MCP-counter 完成，结果维度:", dim(mcp_scores), "\n")
} else {
  cat("MCP-counter 运行失败，请检查基因名格式。\n")
  mcp_scores <- NULL
}

# ========================= 4. 运行 EPIC =======================================
cat("运行 EPIC...\n")
# EPIC 要求表达矩阵为 基因 × 样本，且值应为 TPM/FPKM（非负实数）
epic_res <- tryCatch(
  EPIC::EPIC(
    expr_mat,
    scale = TRUE   # 自动缩放，提高稳定性
  ),
  error = function(e) {
    message("EPIC 运行出错: ", e$message)
    return(NULL)
  }
)

if (!is.null(epic_res)) {
  epic_scores <- as.data.frame(epic_res$cellFractions)
  colnames(epic_scores) <- paste0("EPIC_", colnames(epic_scores))
  cat("EPIC 完成，结果维度:", dim(epic_scores), "\n")
} else {
  cat("EPIC 运行失败，请检查表达矩阵格式。\n")
  epic_scores <- NULL
}

# ========================= 5. 合并结果并保存 ==================================
if (!is.null(mcp_scores) & !is.null(epic_scores)) {
  all_scores <- cbind(mcp_scores, epic_scores)
  all_scores <- rownames_to_column(all_scores, var = "Sample")
  fwrite(all_scores, "Immune_Infiltration_Scores_MCP_EPIC.csv")
  cat("结果已保存至: Immune_Infiltration_Scores_MCP_EPIC.csv\n")
} else if (!is.null(mcp_scores)) {
  fwrite(rownames_to_column(mcp_scores, "Sample"), "MCP_scores_only.csv")
  cat("仅保存 MCP-counter 结果: MCP_scores_only.csv\n")
} else if (!is.null(epic_scores)) {
  fwrite(rownames_to_column(epic_scores, "Sample"), "EPIC_scores_only.csv")
  cat("仅保存 EPIC 结果: EPIC_scores_only.csv\n")
} else {
  stop("两个工具均运行失败，请检查数据格式。")
}

# ========================= 6. 可选：简单可视化 ================================
if (interactive() && !is.null(epic_scores)) {
  library(pheatmap)
  pheatmap(epic_scores, main = "EPIC Immune Cell Fractions", scale = "column")
}
cat("脚本执行完毕。\n")






# ========================= 0. 加载包 ========================================
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(data.table)
library(patchwork)    # 用于组合图形

# ========================= 1. 读取 NETS 评分和临床数据 ========================
nets_file <- "E:/NETS/bulk rna数据/结果/bulk rna数据 4-14/NETS_signature_with_clinical.csv"
nets_data <- fread(nets_file, data.table = FALSE)

# 检查必要的列是否存在
if (!"NETS_score_raw" %in% colnames(nets_data)) {
  stop("文件中没有找到 NETS_score_raw 列，请检查列名。")
}
if (!"sample" %in% colnames(nets_data)) {
  # 如果列名是 Sample 而不是 sample，尝试转换
  if ("Sample" %in% colnames(nets_data)) {
    colnames(nets_data)[colnames(nets_data) == "Sample"] <- "sample"
  } else {
    stop("文件中没有找到样本标识列（sample 或 Sample）。")
  }
}

# 按 NETS_score_raw 中位数分为高/低两组
median_score <- median(nets_data$NETS_score_raw, na.rm = TRUE)
nets_data$NETS_group <- ifelse(nets_data$NETS_score_raw > median_score, "High", "Low")
nets_data$NETS_group <- factor(nets_data$NETS_group, levels = c("Low", "High"))

cat("NETS_score_raw 中位数:", median_score, "\n")
cat("分组样本数:\n")
print(table(nets_data$NETS_group))

# ========================= 2. 读取免疫浸润结果 ================================
immune_file <- "E:/NETS/Immune_Infiltration_Scores_MCP_EPIC.csv"
immune_scores <- fread(immune_file, data.table = FALSE)

# 检查样本名列（假设为 Sample）
if (!"Sample" %in% colnames(immune_scores)) {
  stop("免疫浸润文件中没有找到 Sample 列。")
}

# 统一将样本名列改为 sample，方便合并
colnames(immune_scores)[colnames(immune_scores) == "Sample"] <- "sample"

# ========================= 3. 合并数据 ========================================
merged_data <- merge(nets_data[, c("sample", "NETS_group", "NETS_score_raw")],
                     immune_scores,
                     by = "sample", all = FALSE)

cat("合并后的样本数:", nrow(merged_data), "\n")
if (nrow(merged_data) == 0) {
  stop("合并后无数据，请检查样本名是否匹配。")
}

# ========================= 4. 提取 EPIC 和 MCP 相关列 =========================
epic_cols <- grep("^EPIC_", colnames(merged_data), value = TRUE)
mcp_cols  <- grep("^MCP_", colnames(merged_data), value = TRUE)

cat("EPIC 细胞类型:", epic_cols, "\n")
cat("MCP 细胞类型:", mcp_cols, "\n")

# ========================= 5. 定义绘图函数（箱线图 + 统计） ====================
# 确保加载必要的包
library(dplyr)
library(tidyr)

plot_immune_comparison <- function(data, feature_cols, group_col = "NETS_group", 
                                   title_prefix = "", y_label = "Score / Fraction") {
  # 将数据转换为长格式
  long_data <- data %>%
    dplyr::select(all_of(c(group_col, feature_cols))) %>%   # 明确使用 dplyr::select
    tidyr::pivot_longer(cols = all_of(feature_cols),
                        names_to = "CellType",
                        values_to = "Value")
  
  # 计算统计显著性（Wilcoxon检验）
  stat_results <- long_data %>%
    group_by(CellType) %>%
    summarise(p_value = wilcox.test(Value ~ .data[[group_col]])$p.value,
              .groups = "drop") %>%
    mutate(signif = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    ))
  
  # 绘图
  p <- ggplot(long_data, aes(x = .data[[group_col]], y = Value, fill = .data[[group_col]])) +
    geom_boxplot(outlier.shape = 21, outlier.size = 1.5, alpha = 0.7) +
    geom_jitter(width = 0.2, size = 0.8, alpha = 0.5) +
    facet_wrap(~ CellType, scales = "free_y", ncol = 4) +
    labs(title = title_prefix,
         x = "NETS Score Group", y = y_label) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_rect(fill = "lightgrey"),
          strip.text = element_text(face = "bold")) +
    scale_fill_manual(values = c("Low" = "#2E9FDF", "High" = "#E67E22"))
  
  return(list(plot = p, statistics = stat_results))
}

# ========================= 6. 绘制 EPIC 结果 ==================================
epic_res <- plot_immune_comparison(merged_data, epic_cols, 
                                   title_prefix = "EPIC - Immune Cell Fractions by NETS Group",
                                   y_label = "Estimated Fraction")
print(epic_res$plot)
ggsave("EPIC_boxplot_by_NETS_group.pdf", epic_res$plot, width = 12, height = 8)
cat("EPIC 箱线图已保存为 EPIC_boxplot_by_NETS_group.pdf\n")

# 输出 EPIC 统计结果
epic_stats <- epic_res$statistics
print(epic_stats)
fwrite(epic_stats, "EPIC_statistical_results.csv")

# ========================= 7. 绘制 MCP-counter 结果 ===========================
mcp_res <- plot_immune_comparison(merged_data, mcp_cols,
                                  title_prefix = "MCP-counter - Immune/Stromal Scores by NETS Group",
                                  y_label = "MCP-counter Score")
print(mcp_res$plot)
ggsave("MCP_boxplot_by_NETS_group.pdf", mcp_res$plot, width = 12, height = 8)
cat("MCP-counter 箱线图已保存为 MCP_boxplot_by_NETS_group.pdf\n")

fwrite(mcp_res$statistics, "MCP_statistical_results.csv")

# ========================= 8. （可选）汇总热图 =================================
# 如果需要展示高低两组各细胞类型的平均差异，可以绘制差异热图
mean_diff <- merged_data %>%
  group_by(NETS_group) %>%
  summarise(across(all_of(c(epic_cols, mcp_cols)), mean, na.rm = TRUE)) %>%
  column_to_rownames("NETS_group") %>%
  t() %>%
  as.data.frame() %>%
  mutate(Diff = High - Low) %>%
  rownames_to_column("CellType")

library(pheatmap)
# 准备差异矩阵（仅展示 High - Low 的差值）
diff_matrix <- mean_diff %>%
  select(CellType, Diff) %>%
  column_to_rownames("CellType") %>%
  as.matrix()

pheatmap(diff_matrix,
         cluster_rows = TRUE, cluster_cols = FALSE,
         main = "Difference in Immune Infiltration (High - Low NETS)",
         color = colorRampPalette(c("#2E9FDF", "white", "#E67E22"))(50),
         fontsize_row = 8, width = 6, height = 10,
         filename = "Heatmap_NETS_immune_diff.pdf")

cat("差异热图已保存为 Heatmap_NETS_immune_diff.pdf\n")

cat("\n===== 分析完成 =====")





# ========================= 完整脚本（带p值标注） =========================
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(data.table)
library(patchwork)

# 设置工作目录（请根据实际修改）
setwd("E:/NETS/")

# 1. 读取 NETS 评分和临床数据
nets_file <- "bulk rna数据/结果/bulk rna数据 4-14/NETS_signature_with_clinical.csv"
nets_data <- fread(nets_file, data.table = FALSE)
colnames(nets_data)[tolower(colnames(nets_data)) == "sample"] <- "sample"
nets_data$NETS_group <- ifelse(nets_data$NETS_score_raw > median(nets_data$NETS_score_raw, na.rm=TRUE), "High", "Low")
nets_data$NETS_group <- factor(nets_data$NETS_group, levels = c("Low", "High"))

# 2. 读取免疫浸润结果
immune_file <- "Immune_Infiltration_Scores_MCP_EPIC.csv"
immune_scores <- fread(immune_file, data.table = FALSE)
colnames(immune_scores)[colnames(immune_scores) == "Sample"] <- "sample"

# 3. 合并
merged_data <- merge(nets_data[, c("sample", "NETS_group")], immune_scores, by = "sample")

# 4. 提取细胞类型列
epic_cols <- grep("^EPIC_", colnames(merged_data), value = TRUE)
mcp_cols  <- grep("^MCP_", colnames(merged_data), value = TRUE)

# 5. 绘图函数（带 p 值）
plot_immune_comparison <- function(data, feature_cols, group_col = "NETS_group", 
                                   title_prefix = "", y_label = "Score / Fraction",
                                   p_label = "p.format") {
  long_data <- data %>%
    dplyr::select(all_of(c(group_col, feature_cols))) %>%
    tidyr::pivot_longer(cols = all_of(feature_cols), names_to = "CellType", values_to = "Value")
  
  p <- ggplot(long_data, aes(x = .data[[group_col]], y = Value, fill = .data[[group_col]])) +
    geom_boxplot(outlier.shape = 21, outlier.size = 1.5, alpha = 0.7) +
    geom_jitter(width = 0.2, size = 0.8, alpha = 0.5) +
    facet_wrap(~ CellType, scales = "free_y", ncol = 4) +
    labs(title = title_prefix, x = "NETS Score Group", y = y_label) +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_rect(fill = "lightgrey"), strip.text = element_text(face = "bold")) +
    scale_fill_manual(values = c("Low" = "#2E9FDF", "High" = "#E67E22")) +
    stat_compare_means(comparisons = list(c("Low", "High")),
                       label = p_label, method = "wilcox.test", hide.ns = FALSE,
                       label.y.npc = 0.95)
  
  stat_results <- long_data %>%
    group_by(CellType) %>%
    summarise(p_value = wilcox.test(Value ~ .data[[group_col]])$p.value, .groups = "drop") %>%
    mutate(signif = case_when(p_value < 0.001 ~ "***", p_value < 0.01 ~ "**",
                              p_value < 0.05 ~ "*", TRUE ~ "ns"))
  return(list(plot = p, statistics = stat_results))
}

# 6. 绘制并保存
epic_res <- plot_immune_comparison(merged_data, epic_cols, 
                                   title_prefix = "EPIC - Immune Cell Fractions by NETS Group",
                                   y_label = "Estimated Fraction", p_label = "p.format")
ggsave("EPIC_boxplot_with_pvalue.pdf", epic_res$plot, width = 12, height = 8)

mcp_res <- plot_immune_comparison(merged_data, mcp_cols,
                                  title_prefix = "MCP-counter - Immune/Stromal Scores by NETS Group",
                                  y_label = "MCP-counter Score", p_label = "p.format")
ggsave("MCP_boxplot_with_pvalue.pdf", mcp_res$plot, width = 12, height = 8)

# 保存统计表格
fwrite(epic_res$statistics, "EPIC_stats.csv")
fwrite(mcp_res$statistics, "MCP_stats.csv")

cat("分析完成，带 p 值的箱线图已保存。\n")
