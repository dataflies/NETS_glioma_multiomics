# 重新读取54样本原始FPKM（未log）
file_54_raw <- "E:/NETS/bulk rna数据/原数据/1_CGGA_54_Filtered/CGGA.mRNAseq_325_filtered_genes.txt"
expr_54_raw <- fread(file_54_raw, data.table = FALSE)
rownames(expr_54_raw) <- expr_54_raw[[1]]   # 第一列为基因名
expr_54_raw <- expr_54_raw[, -1]            # 移除基因名列
expr_mat_54 <- as.matrix(expr_54_raw)
mode(expr_mat_54) <- "numeric"
expr_mat_54[is.na(expr_mat_54)] <- 0
expr_mat_54 <- expr_mat_54[rowSums(expr_mat_54) > 0, ]  # 去除全零基因
cat("54样本原始表达矩阵（基因×样本）:", dim(expr_mat_54), "\n")

# 基因名标准化（去除点号后缀，例如 "TP53.1" -> "TP53"）
rownames(expr_mat_54) <- gsub("\\..*", "", rownames(expr_mat_54))

# ==================== MCP ====================
library(MCPcounter)
mcp_res_54 <- tryCatch(
  MCPcounter.estimate(expression = expr_mat_54, featuresType = "HUGO_symbols"),
  error = function(e) { message("MCP-counter 出错: ", e$message); return(NULL) }
)

if (!is.null(mcp_res_54)) {
  mcp_scores_54 <- as.data.frame(t(mcp_res_54))
  colnames(mcp_scores_54) <- paste0("MCP_", colnames(mcp_scores_54))
  mcp_scores_54$sample <- rownames(mcp_scores_54)
  cat("MCP-counter 完成，维度:", dim(mcp_scores_54), "\n")
} else {
  stop("MCP-counter 运行失败")
}

# ==================== EPIC ====================
library(EPIC)
epic_res_54 <- tryCatch(
  EPIC(expr_mat_54, scale = TRUE),
  error = function(e) { message("EPIC 出错: ", e$message); return(NULL) }
)

if (!is.null(epic_res_54)) {
  epic_scores_54 <- as.data.frame(epic_res_54$cellFractions)
  colnames(epic_scores_54) <- paste0("EPIC_", colnames(epic_scores_54))
  epic_scores_54$sample <- rownames(epic_scores_54)
  cat("EPIC 完成，维度:", dim(epic_scores_54), "\n")
} else {
  warning("EPIC 运行失败，将只分析 MCP-counter 结果")
  epic_scores_54 <- NULL
}

# ====================分组可视化 ====================


# 假设 merged_54 已经包含 RiskScore（见之前的生存分析）
# 如果 merged_54 不在当前环境中，可重新读取风险评分文件
risk_file <- "E:/NETS/bulk rna数据/结果/54_validation/54_samples_risk_scores.csv"
if (!exists("merged_54")) {
  risk_data <- fread(risk_file, data.table = FALSE)
  # 读取临床并合并（参考之前生存分析代码）
  clinical_54 <- fread("E:/NETS/bulk rna数据/原数据/1_CGGA_54_Filtered/CGGA.mRNAseq_325_clinical_filtered.txt", data.table = FALSE)
  clinical_54$sample <- clinical_54$CGGA_ID
  clinical_54$OS.status <- as.numeric(clinical_54$`Censor (alive=0; dead=1)`)
  colnames(clinical_54)[colnames(clinical_54) == "OS"] <- "OS.time"
  merged_54 <- merge(risk_data, clinical_54, by.x = "Sample", by.y = "sample")
}

# 按风险评分中位数分组
median_risk <- median(merged_54$RiskScore, na.rm = TRUE)
merged_54$NETS_group <- ifelse(merged_54$RiskScore > median_risk, "High", "Low")
merged_54$NETS_group <- factor(merged_54$NETS_group, levels = c("Low", "High"))
cat("NETS_group 样本数:\n")
print(table(merged_54$NETS_group))

# MCP 数据合并
mcp_plot <- merge(merged_54[, c("Sample", "NETS_group")], 
                  mcp_scores_54, by.x = "Sample", by.y = "sample")
# EPIC 数据合并（如果存在）
if (!is.null(epic_scores_54)) {
  epic_plot <- merge(merged_54[, c("Sample", "NETS_group")], 
                     epic_scores_54, by.x = "Sample", by.y = "sample")
}

library(ggplot2)
library(ggpubr)
library(tidyr)

# 通用绘图函数
plot_immune <- function(data, cell_cols, group_col = "NETS_group", 
                        title = "", ylab = "Score/Fraction") {
  long_data <- data %>%
    dplyr::select(all_of(c(group_col, cell_cols))) %>%   # 明确指定 dplyr::select
    tidyr::pivot_longer(cols = all_of(cell_cols), names_to = "CellType", values_to = "Value")
  
  p <- ggplot(long_data, aes(x = .data[[group_col]], y = Value, fill = .data[[group_col]])) +
    geom_boxplot(outlier.shape = 21, outlier.size = 1.2, alpha = 0.7) +
    geom_jitter(width = 0.2, size = 0.8, alpha = 0.5) +
    facet_wrap(~ CellType, scales = "free_y", ncol = 4) +
    labs(title = title, x = "NETS Risk Group", y = ylab) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_rect(fill = "lightgrey"),
          strip.text = element_text(face = "bold")) +
    scale_fill_manual(values = c("Low" = "#2E9FDF", "High" = "#E67E22")) +
    stat_compare_means(comparisons = list(c("Low", "High")),
                       label = "p.format", method = "wilcox.test", hide.ns = FALSE,
                       label.y.npc = 0.95)
  return(p)
}

# 绘制 MCP 结果
mcp_cols <- grep("^MCP_", colnames(mcp_plot), value = TRUE)
p_mcp <- plot_immune(mcp_plot, mcp_cols, title = "MCP-counter Scores by NETS Group (54 samples)")
ggsave("54_MCP_boxplot.pdf", p_mcp, width = 12, height = 8)
cat("MCP 箱线图已保存: 54_MCP_boxplot.pdf\n")

# 绘制 EPIC 结果（若有）
if (exists("epic_plot")) {
  epic_cols <- grep("^EPIC_", colnames(epic_plot), value = TRUE)
  p_epic <- plot_immune(epic_plot, epic_cols, title = "EPIC Cell Fractions by NETS Group (54 samples)", ylab = "Estimated Fraction")
  ggsave("54_EPIC_boxplot.pdf", p_epic, width = 12, height = 8)
  cat("EPIC 箱线图已保存: 54_EPIC_boxplot.pdf\n")
}


# MCP 统计（修正 select 冲突）
mcp_stats <- mcp_plot %>%
  dplyr::select(all_of(c("NETS_group", mcp_cols))) %>%
  tidyr::pivot_longer(cols = all_of(mcp_cols), names_to = "CellType", values_to = "Value") %>%
  group_by(CellType) %>%
  summarise(p_value = wilcox.test(Value ~ NETS_group)$p.value, .groups = "drop") %>%
  mutate(Significance = case_when(p_value < 0.001 ~ "***", p_value < 0.01 ~ "**",
                                  p_value < 0.05 ~ "*", TRUE ~ "ns"))

fwrite(mcp_stats, "54_MCP_statistics.csv")
cat("MCP 统计表已保存: 54_MCP_statistics.csv\n")

# EPIC 统计（如果存在）
if (exists("epic_plot")) {
  epic_stats <- epic_plot %>%
    dplyr::select(all_of(c("NETS_group", epic_cols))) %>%
    tidyr::pivot_longer(cols = all_of(epic_cols), names_to = "CellType", values_to = "Value") %>%
    group_by(CellType) %>%
    summarise(p_value = wilcox.test(Value ~ NETS_group)$p.value, .groups = "drop") %>%
    mutate(Significance = case_when(p_value < 0.001 ~ "***", p_value < 0.01 ~ "**",
                                    p_value < 0.05 ~ "*", TRUE ~ "ns"))
  fwrite(epic_stats, "54_EPIC_statistics.csv")
  cat("EPIC 统计表已保存: 54_EPIC_statistics.csv\n")
}

library(pheatmap)

# 计算MCP细胞丰度高低组均值差异
# 计算 MCP 细胞丰度高低组均值差异
mcp_mean <- mcp_plot %>%
  group_by(NETS_group) %>%
  summarise(across(all_of(mcp_cols), mean, na.rm = TRUE)) %>%
  column_to_rownames("NETS_group") %>%
  t() %>%
  as.data.frame() %>%
  mutate(Diff = High - Low) %>%
  rownames_to_column("CellType")

# 构建差异矩阵（不使用 select，避免冲突）
diff_vec <- mcp_mean$Diff
names(diff_vec) <- mcp_mean$CellType
diff_mat <- as.matrix(diff_vec)
colnames(diff_mat) <- "High - Low"

# 可选：按差异值排序（从大到小）
diff_mat <- diff_mat[order(diff_mat[,1], decreasing = TRUE), , drop = FALSE]

# 绘制热图
pheatmap(diff_mat, 
         cluster_rows = TRUE,   # 虽然行列数>1才能聚类，但这里只有一列，会自动忽略列聚类
         cluster_cols = FALSE,
         main = "MCP-score Difference (High - Low NETS group, 54 samples)",
         color = colorRampPalette(c("#2E9FDF", "white", "#E67E22"))(50),
         filename = "54_MCP_heatmap.pdf", width = 6, height = 8)

if (exists("epic_plot")) {
  epic_mean <- epic_plot %>%
    group_by(NETS_group) %>%
    summarise(across(all_of(epic_cols), mean, na.rm = TRUE)) %>%
    column_to_rownames("NETS_group") %>%
    t() %>%
    as.data.frame() %>%
    mutate(Diff = High - Low) %>%
    rownames_to_column("CellType")
  
  diff_vec_epic <- epic_mean$Diff
  names(diff_vec_epic) <- epic_mean$CellType
  diff_mat_epic <- as.matrix(diff_vec_epic)
  colnames(diff_mat_epic) <- "High - Low"
  diff_mat_epic <- diff_mat_epic[order(diff_mat_epic[,1], decreasing = TRUE), , drop = FALSE]
  
  pheatmap(diff_mat_epic,
           cluster_rows = TRUE, cluster_cols = FALSE,
           main = "EPIC Cell Fraction Difference (High - Low NETS group, 54 samples)",
           color = colorRampPalette(c("#2E9FDF", "white", "#E67E22"))(50),
           filename = "54_EPIC_heatmap.pdf", width = 6, height = 8)
}

# ==================== 修正：确保需要的变量和包 ====================
library(pheatmap)
library(tidyverse)   # 包含 tibble 的 column_to_rownames

# 如果 mcp_cols 或 epic_cols 未定义，重新提取
if (!exists("mcp_cols")) {
  mcp_cols <- grep("^MCP_", colnames(mcp_plot), value = TRUE)
}
if (exists("epic_plot") && !exists("epic_cols")) {
  epic_cols <- grep("^EPIC_", colnames(epic_plot), value = TRUE)
}

# 构建 MCP 矩阵
mcp_mat <- mcp_plot %>%
  select(Sample, all_of(mcp_cols)) %>%
  column_to_rownames("Sample")

# 构建 EPIC 矩阵（如果存在）并合并
if (exists("epic_plot") && !is.null(epic_plot)) {
  epic_mat <- epic_plot %>%
    select(Sample, all_of(epic_cols)) %>%
    column_to_rownames("Sample")
  # 合并两个矩阵（按行名，即样本名）
  immune_mat <- merge(mcp_mat, epic_mat, by = "row.names")
  rownames(immune_mat) <- immune_mat$Row.names
  immune_mat <- immune_mat[, -1]
} else {
  immune_mat <- mcp_mat
}

# 计算 Spearman 相关系数矩阵
cor_matrix <- cor(immune_mat, method = "spearman", use = "pairwise.complete.obs")

# 定义颜色（-1 到 1）
my_colors <- colorRampPalette(c("#2E9FDF", "white", "#E67E22"))(100)

# 绘制热图并保存
pdf("Immune_cell_correlation_heatmap_54.pdf", width = 10, height = 8)
pheatmap(cor_matrix,
         color = my_colors,
         breaks = seq(-1, 1, length.out = 101),
         main = "Spearman correlation among immune cell fractions (n=54)",
         display_numbers = TRUE,
         number_format = "%.2f",
         fontsize_number = 6,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         treeheight_row = 30,
         treeheight_col = 30)
dev.off()

cat("热图已保存为: Immune_cell_correlation_heatmap_54.pdf\n")


