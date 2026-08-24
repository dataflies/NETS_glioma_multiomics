# ========================= 0. 加载必需的 R 包 =========================
library(data.table)
library(tidyverse)
library(sva)      # ComBat 批次校正
library(limma)    # 可选：数据规范化

# ========================= 1. 读取原始数据 ============================
# 54 样本数据
file_54 <- "E:/NETS/bulk rna数据/原数据/1_CGGA_54_Filtered/CGGA.mRNAseq_325_filtered_genes.txt"
expr_54_raw <- fread(file_54, data.table = FALSE)
rownames(expr_54_raw) <- expr_54_raw[, 1]
expr_54_raw <- expr_54_raw[, -1]
expr_54_raw <- as.matrix(expr_54_raw)
mode(expr_54_raw) <- "numeric"
cat("54 样本数据维度（基因 × 样本）:", dim(expr_54_raw), "\n")

# 删除 CGGA_256 样本（如果存在）
if ("CGGA_256" %in% colnames(expr_54_raw)) {
  expr_54_raw <- expr_54_raw[, colnames(expr_54_raw) != "CGGA_256", drop = FALSE]
  cat("已从 54 样本数据中删除 CGGA_256 列。\n")
} else {
  cat("54 样本数据中未发现 CGGA_256 列，无需删除。\n")
}

# 118 样本数据
file_118 <- "E:/NETS/bulk rna数据/原数据/1_CGGA_118_Filtered/CGGA.mRNAseq_693_filtered.txt"
expr_118_raw <- fread(file_118, data.table = FALSE)
rownames(expr_118_raw) <- expr_118_raw[, 1]
expr_118_raw <- expr_118_raw[, -1]
expr_118_raw <- as.matrix(expr_118_raw)
mode(expr_118_raw) <- "numeric"
cat("118 样本数据维度（基因 × 样本）:", dim(expr_118_raw), "\n")

# 删除 CGGA_256 样本（如果存在）
if ("CGGA_256" %in% colnames(expr_118_raw)) {
  expr_118_raw <- expr_118_raw[, colnames(expr_118_raw) != "CGGA_256", drop = FALSE]
  cat("已从 118 样本数据中删除 CGGA_256 列。\n")
} else {
  cat("118 样本数据中未发现 CGGA_256 列，无需删除。\n")
}

# ========================= 2. 取基因交集 ============================
common_genes <- intersect(rownames(expr_54_raw), rownames(expr_118_raw))
cat("共同基因数量:", length(common_genes), "\n")

expr_54 <- expr_54_raw[common_genes, , drop = FALSE]
expr_118 <- expr_118_raw[common_genes, , drop = FALSE]

# ========================= 3. log2 转换 ==============================
# 避免 log2(0) -> 使用 log2(FPKM + 1)
log_expr_54 <- log2(expr_54 + 1)
log_expr_118 <- log2(expr_118 + 1)

# ========================= 4. ComBat 批次校正 =======================
combined_log <- cbind(log_expr_54, log_expr_118)
batch <- c(rep(1, ncol(expr_54)), rep(2, ncol(expr_118)))   # 批次标签
combined_combat <- ComBat(dat = combined_log, batch = batch, par.prior = TRUE)

# 拆分回各自数据集
expr_54_combat <- combined_combat[, 1:ncol(expr_54)]
expr_118_combat <- combined_combat[, (ncol(expr_54)+1):ncol(combined_combat)]

cat("ComBat 批次校正完成。\n")

# ========================= 5. 检查负值和 NA =========================
# 理论上 log2(FPKM+1) 后数据正态分布，ComBat 可能产生少量负值
cat("54 样本中负值数量:", sum(expr_54_combat < 0), "\n")
cat("118 样本中负值数量:", sum(expr_118_combat < 0), "\n")

# NA 检查
if (any(is.na(expr_54_combat))) expr_54_combat[is.na(expr_54_combat)] <- 0
if (any(is.na(expr_118_combat))) expr_118_combat[is.na(expr_118_combat)] <- 0

# ========================= 6. 提取风险评分基因 ======================
# 假设弹性网络分析结果保存在 elasticnet_cox_results.RData
load("E:/NETS/elasticnet_cox_results.RData")   # 包含 beta_final 和 final_genes

my_genes <- c("IFI30", "DCN", "DUSP1", "MS4A6A", "CRTAP", "TAGLN")
my_genes <- intersect(my_genes, final_genes)  # 保证基因存在
cat("用于风险评分的基因:", my_genes, "\n")

# 提取系数并确保顺序一致
beta_vec <- beta_final[my_genes]
beta_vec <- beta_vec[my_genes]
cat("系数:\n"); print(beta_vec)

# 提取基因表达矩阵（样本 × 基因）
expr_sig_54 <- t(expr_54_combat[my_genes, , drop = FALSE])
expr_sig_118 <- t(expr_118_combat[my_genes, , drop = FALSE])

# ========================= 7. 计算风险评分 ==========================
risk_raw_54 <- as.numeric(as.matrix(expr_sig_54) %*% beta_vec)
risk_raw_118 <- as.numeric(as.matrix(expr_sig_118) %*% beta_vec)

# 0-1 归一化函数
normalize_01 <- function(x) (x - min(x, na.rm=TRUE)) / (max(x, na.rm=TRUE) - min(x, na.rm=TRUE))

risk_01_54 <- normalize_01(risk_raw_54)
risk_01_118 <- normalize_01(risk_raw_118)

# scale 标准化
risk_scale_54 <- scale(risk_raw_54)[,1]
risk_scale_118 <- scale(risk_raw_118)[,1]

# ========================= 8. 保存结果 ==============================
result_54 <- data.frame(
  sample = colnames(expr_54),
  risk_raw = risk_raw_54,
  risk_01 = risk_01_54,
  risk_scale = risk_scale_54
)
result_118 <- data.frame(
  sample = colnames(expr_118),
  risk_raw = risk_raw_118,
  risk_01 = risk_01_118,
  risk_scale = risk_scale_118
)

fwrite(result_54, "CGGA_54_risk_scores.csv")
fwrite(result_118, "CGGA_118_risk_scores.csv")

cat("风险评分已保存完成。\n")
