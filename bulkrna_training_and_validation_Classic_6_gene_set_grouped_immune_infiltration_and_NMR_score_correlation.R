# ========================= 0. 安装和加载所需包 =================================
# 首次使用请取消注释并运行
# install.packages(c("tidyverse", "data.table", "ggplot2", "ggpubr", "patchwork"))
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("MCPcounter")
# pak::pkg_install("GfellerLab/EPIC")   # EPIC 从 GitHub 安装
# install.packages("UCell")             # UCell 包

library(MCPcounter)
library(EPIC)
library(UCell)          # 基因集评分
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(patchwork)

# 设置工作目录（请修改为您的实际路径）
setwd("E:/NETS/")

# ========================= 1. 定义基因集 ====================================
geneset <- c("PTGS2", "MME", "SLC2A3", "MPO", "ELANE", "PADI4")
# 若 ELANE/PADI4 表示两个基因，按上述写法；若是一个融合基因，请改为 "ELANE/PADI4"

# ========================= 2. 读取两个批次的表达矩阵 =========================
# 假设两个表达矩阵文件分别为：
#   - 118 批次：bulk rna数据/原数据/1_CGGA_118_Filtered/CGGA.mRNAseq_693_filtered.txt
#   - 54 批次：bulk rna数据/原数据/2_CGGA_54_Filtered/CGGA.mRNAseq_325_filtered.txt
# 如果只有单一表达矩阵，可只读取一个，但按照第二组代码，应该是两个独立数据集。

expr_file_118 <- "bulk rna数据/原数据/1_CGGA_118_Filtered/CGGA.mRNAseq_693_filtered.txt"
expr_file_54  <- "bulk rna数据/原数据/1_CGGA_54_Filtered/CGGA.mRNAseq_325_filtered_genes.txt"

# 读取并整理表达矩阵的函数
read_expr_matrix <- function(file_path) {
  expr <- fread(file = file_path, data.table = FALSE)   # 关键修改：加上 file=
  rownames(expr) <- expr[, 1]
  expr <- expr[, -1]
  expr_mat <- as.matrix(expr)
  mode(expr_mat) <- "numeric"
  expr_mat[is.na(expr_mat)] <- 0
  expr_mat <- expr_mat[rowSums(expr_mat) > 0, ]
  rownames(expr_mat) <- gsub("\\..*", "", rownames(expr_mat))
  return(expr_mat)
}

# 重新读取两个文件
expr_118 <- read_expr_matrix(expr_file_118)
expr_54  <- read_expr_matrix(expr_file_54)
cat("118 批次表达矩阵维度:", dim(expr_118), "\n")
cat("54 批次表达矩阵维度:", dim(expr_54), "\n")

## ==================== 使用平均表达量作为评分 ====================
# 对每个样本，计算 geneset 中基因表达量的平均值（或中位数）

# 计算 118 批次
avg_expr_118 <- colMeans(expr_118[geneset, , drop = FALSE], na.rm = TRUE)
# 如果担心极端值，也可以使用中位数
# avg_expr_118 <- apply(expr_118[geneset, , drop = FALSE], 2, median, na.rm = TRUE)

# 计算 54 批次
avg_expr_54 <- colMeans(expr_54[geneset, , drop = FALSE], na.rm = TRUE)

# 构建数据框（变量名保持 UCell_score 以兼容后续代码）
ucell_df_118 <- data.frame(
  Sample = colnames(expr_118),
  UCell_score = avg_expr_118,
  stringsAsFactors = FALSE
)

ucell_df_54 <- data.frame(
  Sample = colnames(expr_54),
  UCell_score = avg_expr_54,
  stringsAsFactors = FALSE
)

# 检查分布
summary(ucell_df_118$UCell_score)
summary(ucell_df_54$UCell_score)

median_118 <- median(ucell_df_118$UCell_score, na.rm = TRUE)
median_54  <- median(ucell_df_54$UCell_score, na.rm = TRUE)
cat("118 中位数:", median_118, "\n")
cat("54 中位数:", median_54, "\n")
# 查看分组数量
table(ifelse(ucell_df_118$UCell_score > median_118, "High", "Low"))
table(ifelse(ucell_df_54$UCell_score > median_54, "High", "Low"))



# ========================= 4. 运行 MCP-counter 和 EPIC ========================
# 注意：MCPcounter 和 EPIC 需要原始的（非 log 转换）TPM/FPKM 值，且基因名应为 HUGO 符号。

run_immune_deconvolution <- function(expr_mat, batch_name) {
  cat("运行", batch_name, "的 MCP-counter...\n")
  mcp_res <- tryCatch(
    MCPcounter.estimate(expression = expr_mat, featuresType = "HUGO_symbols"),
    error = function(e) { message("MCP 失败: ", e$message); return(NULL) }
  )
  if (!is.null(mcp_res)) {
    mcp_scores <- as.data.frame(t(mcp_res))
    colnames(mcp_scores) <- paste0("MCP_", colnames(mcp_scores))
    mcp_scores <- rownames_to_column(mcp_scores, var = "Sample")
  } else {
    mcp_scores <- NULL
  }
  
  cat("运行", batch_name, "的 EPIC...\n")
  epic_res <- tryCatch(
    EPIC::EPIC(expr_mat, scale = TRUE),
    error = function(e) { message("EPIC 失败: ", e$message); return(NULL) }
  )
  if (!is.null(epic_res)) {
    epic_scores <- as.data.frame(epic_res$cellFractions)
    colnames(epic_scores) <- paste0("EPIC_", colnames(epic_scores))
    epic_scores <- rownames_to_column(epic_scores, var = "Sample")
  } else {
    epic_scores <- NULL
  }
  
  # 合并免疫结果
  if (!is.null(mcp_scores) & !is.null(epic_scores)) {
    immune_all <- merge(mcp_scores, epic_scores, by = "Sample", all = TRUE)
  } else if (!is.null(mcp_scores)) {
    immune_all <- mcp_scores
  } else if (!is.null(epic_scores)) {
    immune_all <- epic_scores
  } else {
    stop("两个工具均运行失败")
  }
  return(immune_all)
}

immune_118 <- run_immune_deconvolution(expr_118, "118")
immune_54  <- run_immune_deconvolution(expr_54, "54")

# ========================= 5. 读取 NETS 风险评分（来自第二组代码） ============
# 假设风险评分文件为：
#   - 1_CGGA_118_risk_scores.csv
#   - 1_CGGA_54_risk_scores.csv
# 它们包含 Sample 和 risk_raw 列

risk_118 <- read.csv("1_CGGA_118_risk_scores.csv", stringsAsFactors = FALSE)
risk_118 <- risk_118[, c("sample", "risk_raw")]
risk_54  <- read.csv("1_CGGA_54_risk_scores.csv", stringsAsFactors = FALSE)
risk_54  <- risk_54[, c("sample", "risk_raw")]

# ========================= 6. 合并所有数据（按批次） ==========================
# 修正风险数据的列名
risk_118 <- read.csv("1_CGGA_118_risk_scores.csv", stringsAsFactors = FALSE)
colnames(risk_118)[1] <- "Sample"   # 假设第一列为样本ID
risk_118 <- risk_118[, c("Sample", "risk_raw")]

risk_54 <- read.csv("1_CGGA_54_risk_scores.csv", stringsAsFactors = FALSE)
colnames(risk_54)[1] <- "Sample"
risk_54 <- risk_54[, c("Sample", "risk_raw")]

# 分别合并两个批次（不合并在一起）
data_118 <- merge_batch_data(ucell_df_118, immune_118, risk_118, "118")
data_54  <- merge_batch_data(ucell_df_54, immune_54, risk_54, "54")

# ========================= 7. 分别对每个批次进行相关性分析 =====================
# 定义一个函数，对单个批次数据执行全套分析（相关性 + 可视化）
analyze_batch <- function(data, batch_name) {
  # 提取免疫列
  mcp_cols  <- grep("^MCP_", colnames(data), value = TRUE)
  epic_cols <- grep("^EPIC_", colnames(data), value = TRUE)
  
  # 计算 UCell vs 免疫细胞相关性
  cor_mcp <- data.frame(CellType = character(), R = numeric(), P = numeric(), stringsAsFactors = FALSE)
  for (col in mcp_cols) {
    if (sum(complete.cases(data$UCell_score, data[[col]])) >= 3) {
      test <- cor.test(data$UCell_score, data[[col]], method = "spearman", use = "complete.obs")
      cor_mcp <- rbind(cor_mcp, data.frame(CellType = col, R = test$estimate, P = test$p.value))
    }
  }
  
  cor_epic <- data.frame(CellType = character(), R = numeric(), P = numeric(), stringsAsFactors = FALSE)
  for (col in epic_cols) {
    if (sum(complete.cases(data$UCell_score, data[[col]])) >= 3) {
      test <- cor.test(data$UCell_score, data[[col]], method = "spearman", use = "complete.obs")
      cor_epic <- rbind(cor_epic, data.frame(CellType = col, R = test$estimate, P = test$p.value))
    }
  }
  
  # UCell vs NETS 风险
  if (sum(complete.cases(data$UCell_score, data$risk_raw)) >= 3) {
    cor_net <- cor.test(data$UCell_score, data$risk_raw, method = "spearman", use = "complete.obs")
    cor_net_df <- data.frame(Variable1 = "UCell_score", Variable2 = "risk_raw",
                             R = cor_net$estimate, P = cor_net$p.value)
  } else {
    cor_net_df <- data.frame(Variable1 = "UCell_score", Variable2 = "risk_raw", R = NA, P = NA)
  }
  
  # 保存结果
  fwrite(cor_mcp,    paste0("Correlation_UCell_vs_MCP_", batch_name, ".csv"))
  fwrite(cor_epic,   paste0("Correlation_UCell_vs_EPIC_", batch_name, ".csv"))
  fwrite(cor_net_df, paste0("Correlation_UCell_vs_NETS_", batch_name, ".csv"))
  
  # 返回结果列表（方便后续可视化）
  return(list(mcp = cor_mcp, epic = cor_epic, net = cor_net_df))
}

# 执行分析
res_118 <- analyze_batch(data_118, "118")
res_54  <- analyze_batch(data_54, "54")
# ---- 可视化 ----
# 散点图：UCell vs NETS
# 118 批次
p1 <- ggplot(data_118, aes(x = UCell_score, y = risk_raw)) +
  geom_point(color = "#2E9FDF", size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  stat_cor(method = "spearman", label.x.npc = 0.05, label.y.npc = 0.95) +
  labs(title = "118 Cohort: UCell vs NETS Risk", x = "UCell Score", y = "NETS Risk Score") +
  theme_minimal()
ggsave("Scatter_UCell_vs_NETS_118.pdf", p1, width = 5, height = 4)

# 54 批次
p2 <- ggplot(data_54, aes(x = UCell_score, y = risk_raw)) +
  geom_point(color = "#E67E22", size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  stat_cor(method = "spearman", label.x.npc = 0.05, label.y.npc = 0.95) +
  labs(title = "54 Cohort: UCell vs NETS Risk", x = "UCell Score", y = "NETS Risk Score") +
  theme_minimal()
ggsave("Scatter_UCell_vs_NETS_54.pdf", p2, width = 5, height = 4)

# 合并显示（分面）
data_combined <- bind_rows(
  data_118 %>% mutate(Cohort = "118"),
  data_54  %>% mutate(Cohort = "54")
)
p_combined <- ggplot(data_combined, aes(x = UCell_score, y = risk_raw, color = Cohort)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
  stat_cor(method = "spearman", label.x.npc = 0.05, label.y.npc = 0.95, aes(color = Cohort)) +
  facet_wrap(~ Cohort, scales = "free") +
  labs(title = "UCell vs NETS Risk by Cohort", x = "UCell Score", y = "NETS Risk Score") +
  theme_minimal() +
  theme(legend.position = "none")
ggsave("Scatter_UCell_vs_NETS_split.pdf", p_combined, width = 10, height = 4)

# 条形图：免疫相关性
# 合并两个批次的相关性结果
cor_compare <- bind_rows(
  res_118$mcp %>% mutate(Cohort = "118", Type = "MCP"),
  res_54$mcp  %>% mutate(Cohort = "54",  Type = "MCP"),
  res_118$epic %>% mutate(Cohort = "118", Type = "EPIC"),
  res_54$epic  %>% mutate(Cohort = "54",  Type = "EPIC")
)

# 绘制分面条形图
p_bar <- ggplot(cor_compare, aes(x = reorder(CellType, R), y = R, fill = Cohort)) +
  geom_col(position = "dodge") +
  facet_wrap(~ Type, scales = "free_x") +
  coord_flip() +
  labs(title = "Spearman Correlation (UCell vs Immune Cells)", 
       x = "Cell Type", y = "R") +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("118" = "#2E9FDF", "54" = "#E67E22"))
ggsave("Barplot_Correlation_UCell_Immune_split.pdf", p_bar, width = 10, height = 8)
  # 可选：散点图矩阵（选取部分细胞）—— 如果列存在
  selected_cells <- c("MCP_T_cells", "MCP_B_cells", "MCP_Monocytes", 
                      "EPIC_CD8_T_cells", "EPIC_B_cells", "EPIC_Macrophages")
  selected_cols <- intersect(selected_cells, colnames(data))
  if (length(selected_cols) > 0) {
    library(GGally)
    p_pairs <- ggpairs(data, columns = c("UCell_score", selected_cols), 
                       upper = list(continuous = wrap("cor", method = "spearman"))) +
      theme_minimal()
    ggsave(paste0("Pairs_UCell_Immune_", batch_name, ".pdf"), p_pairs, width = 12, height = 12)
  }
}

# 对两个批次分别执行分析
analyze_batch(data_118, "118")
analyze_batch(data_54, "54")

# ========================= 加载必要的包 =====================================
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(data.table)
library(patchwork)   # 如需组合图形

# ========================= 定义绘图函数（基于 UCell 分组） ===================
# 参数：
#   data: 包含 UCell_score 和免疫列的数据框
#   feature_cols: 免疫细胞列名（MCP_* 或 EPIC_*）
#   group_col: 分组列名（这里固定为 "UCell_group"）
#   title_prefix: 图标题前缀
#   y_label: y轴标签
#   p_label: p值显示格式（"p.format" 或 "p.signif"）

plot_immune_by_ucell <- function(data, feature_cols, 
                                 title_prefix = "", y_label = "Score / Fraction",
                                 p_label = "p.format") {
  # 1. 根据 UCell_score 中位数分组
  median_ucell <- median(data$UCell_score, na.rm = TRUE)
  data <- data %>%
    mutate(UCell_group = ifelse(UCell_score > median_ucell, "High", "Low"),
           UCell_group = factor(UCell_group, levels = c("Low", "High")))
  
  cat("UCell 中位数:", median_ucell, "\n")
  cat("分组样本数:\n")
  print(table(data$UCell_group))
  
  # 2. 转换为长格式
  long_data <- data %>%
    dplyr::select(all_of(c("UCell_group", feature_cols))) %>%
    tidyr::pivot_longer(cols = all_of(feature_cols),
                        names_to = "CellType",
                        values_to = "Value")
  
  # 3. 计算统计显著性（Wilcoxon检验）
  stat_results <- long_data %>%
    group_by(CellType) %>%
    summarise(p_value = wilcox.test(Value ~ UCell_group)$p.value,
              .groups = "drop") %>%
    mutate(signif = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    ))
  
  # 4. 绘图（使用 ggpubr::stat_compare_means 添加 p 值）
  p <- ggplot(long_data, aes(x = UCell_group, y = Value, fill = UCell_group)) +
    geom_boxplot(outlier.shape = 21, outlier.size = 1.5, alpha = 0.7) +
    geom_jitter(width = 0.2, size = 0.8, alpha = 0.5) +
    facet_wrap(~ CellType, scales = "free_y", ncol = 4) +
    labs(title = title_prefix,
         x = "UCell Score Group", y = y_label) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_rect(fill = "lightgrey"),
          strip.text = element_text(face = "bold")) +
    scale_fill_manual(values = c("Low" = "#2E9FDF", "High" = "#E67E22")) +
    stat_compare_means(comparisons = list(c("Low", "High")),
                       label = p_label, method = "wilcox.test", hide.ns = FALSE,
                       label.y.npc = 0.95)   # 将p值放在每个分面的顶部
  
  return(list(plot = p, statistics = stat_results, group_data = data))
}

# ========================= 为两个批次分别绘制 ================================

# ---------- 118 批次 ----------
# 确保 data_118 已存在（包含 UCell_score, MCP_*, EPIC_*）
mcp_cols_118 <- grep("^MCP_", colnames(data_118), value = TRUE)
epic_cols_118 <- grep("^EPIC_", colnames(data_118), value = TRUE)

# MCP 箱线图
mcp_res_118 <- plot_immune_by_ucell(
  data_118, 
  feature_cols = mcp_cols_118,
  title_prefix = "118 Cohort: MCP-counter Scores by UCell Group",
  y_label = "MCP Score",
  p_label = "p.format"
)
ggsave("MCP_boxplot_by_UCell_118.pdf", mcp_res_118$plot, width = 12, height = 8)
fwrite(mcp_res_118$statistics, "MCP_stats_UCell_118.csv")

# EPIC 箱线图
epic_res_118 <- plot_immune_by_ucell(
  data_118,
  feature_cols = epic_cols_118,
  title_prefix = "118 Cohort: EPIC Fractions by UCell Group",
  y_label = "EPIC Fraction",
  p_label = "p.format"
)
ggsave("EPIC_boxplot_by_UCell_118.pdf", epic_res_118$plot, width = 12, height = 8)
fwrite(epic_res_118$statistics, "EPIC_stats_UCell_118.csv")

# ---------- 54 批次 ----------
mcp_cols_54 <- grep("^MCP_", colnames(data_54), value = TRUE)
epic_cols_54 <- grep("^EPIC_", colnames(data_54), value = TRUE)

mcp_res_54 <- plot_immune_by_ucell(
  data_54,
  feature_cols = mcp_cols_54,
  title_prefix = "54 Cohort: MCP-counter Scores by UCell Group",
  y_label = "MCP Score",
  p_label = "p.format"
)
ggsave("MCP_boxplot_by_UCell_54.pdf", mcp_res_54$plot, width = 12, height = 8)
fwrite(mcp_res_54$statistics, "MCP_stats_UCell_54.csv")

epic_res_54 <- plot_immune_by_ucell(
  data_54,
  feature_cols = epic_cols_54,
  title_prefix = "54 Cohort: EPIC Fractions by UCell Group",
  y_label = "EPIC Fraction",
  p_label = "p.format"
)
ggsave("EPIC_boxplot_by_UCell_54.pdf", epic_res_54$plot, width = 12, height = 8)
fwrite(epic_res_54$statistics, "EPIC_stats_UCell_54.csv")

cat("\n===== 所有箱线图已生成并保存 =====\n")

# ========================= 导出基因集评分为 Excel 文件 =========================
# 安装并加载必要的包（如果尚未安装）
# install.packages("openxlsx")   # 或使用 writexl
library(openxlsx)

# 1. 分别保存两个批次
write.xlsx(ucell_df_118, "GeneSignature_Score_118.xlsx", rowNames = FALSE)
write.xlsx(ucell_df_54,  "GeneSignature_Score_54.xlsx",  rowNames = FALSE)

# 2. （可选）合并为一个文件，用两个 sheet 分别存放
wb <- createWorkbook()
addWorksheet(wb, "Batch_118")
writeData(wb, "Batch_118", ucell_df_118)
addWorksheet(wb, "Batch_54")
writeData(wb, "Batch_54", ucell_df_54)
saveWorkbook(wb, "GeneSignature_Scores_AllBatches.xlsx", overwrite = TRUE)

cat("基因集评分已导出为 Excel 文件。\n")