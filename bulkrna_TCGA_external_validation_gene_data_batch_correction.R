# ========================= 0. 加载必需 R 包 =========================
library(data.table)
library(tidyverse)
library(sva)      # ComBat 批次校正

# ========================= 1. 读取参考数据（CGGA_118） ===================
file_ref <- "E:/NETS/bulk rna数据/原数据/1_CGGA_118_Filtered/CGGA.mRNAseq_693_filtered.txt"
expr_ref_raw <- fread(file_ref, data.table = FALSE)
rownames(expr_ref_raw) <- expr_ref_raw[, 1]   # 第一列为基因名
expr_ref_raw <- expr_ref_raw[, -1]            # 去掉基因名列
expr_ref_raw <- as.matrix(expr_ref_raw)
mode(expr_ref_raw) <- "numeric"
cat("参考数据（118样本）维度（基因 × 样本）:", dim(expr_ref_raw), "\n")

# 删除 CGGA_256 样本（若存在）
if ("CGGA_256" %in% colnames(expr_ref_raw)) {
  expr_ref_raw <- expr_ref_raw[, colnames(expr_ref_raw) != "CGGA_256", drop = FALSE]
  cat("已从参考数据中删除 CGGA_256 列。\n")
}

# ========================= 2. 读取待校正数据（TCGA） =====================
file_tcga <- "E:/NETS/bulk&影像 TCGA外部验证/原数据/1_TCGA_Filtered_Data/gene_filtered.csv"
expr_tcga_raw <- fread(file_tcga, data.table = FALSE)
rownames(expr_tcga_raw) <- expr_tcga_raw[, 1]   # 第一列为基因名
expr_tcga_raw <- expr_tcga_raw[, -1]            # 去掉基因名列
expr_tcga_raw <- as.matrix(expr_tcga_raw)
mode(expr_tcga_raw) <- "numeric"
cat("TCGA数据维度（基因 × 样本）:", dim(expr_tcga_raw), "\n")

# ========================= 3. 取基因交集 ==============================
common_genes <- intersect(rownames(expr_ref_raw), rownames(expr_tcga_raw))
cat("共同基因数量:", length(common_genes), "\n")

# 仅保留共同基因，并按相同顺序排列（参考数据顺序为准）
expr_ref <- expr_ref_raw[common_genes, , drop = FALSE]
expr_tcga <- expr_tcga_raw[common_genes, , drop = FALSE]

# ========================= 4. log2(expr + 1) 转换 ====================
log_ref <- log2(expr_ref + 1)
log_tcga <- log2(expr_tcga + 1)

# ========================= 5. 合并数据并构建批次向量 ==================
combined_log <- cbind(log_ref, log_tcga)
batch <- c(rep(1, ncol(log_ref)), rep(2, ncol(log_tcga)))   # 1=参考，2=TCGA

# ========================= 6. ComBat 批次校正（以第1批次为参考） ======
# 使用 ref.batch 参数指定参考批次（1 为参考），其他批次将校正至参考批次分布
combined_combat <- ComBat(dat = combined_log, 
                          batch = batch, 
                          ref.batch = 1,        # 以 CGGA_118 为参考
                          par.prior = TRUE)

# 拆分校正后的数据
expr_ref_combat <- combined_combat[, 1:ncol(log_ref), drop = FALSE]
expr_tcga_combat <- combined_combat[, (ncol(log_ref) + 1):ncol(combined_combat), drop = FALSE]

cat("ComBat 批次校正完成（参考批次 = CGGA_118）。\n")

# ========================= 7. 检查校正后的异常值 ======================
cat("校正后的参考数据中负值数量:", sum(expr_ref_combat < 0), "\n")
cat("校正后的TCGA数据中负值数量:", sum(expr_tcga_combat < 0), "\n")
# 处理 NA（若存在）
if (any(is.na(expr_ref_combat))) expr_ref_combat[is.na(expr_ref_combat)] <- 0
if (any(is.na(expr_tcga_combat))) expr_tcga_combat[is.na(expr_tcga_combat)] <- 0

# ========================= 8. 保存校正后的 TCGA 表达矩阵 ==============
# 将校正后 TCGA 数据转换为数据框，并添加基因名作为第一列
tcga_out <- as.data.frame(expr_tcga_combat)
tcga_out <- cbind(gene = rownames(tcga_out), tcga_out)
fwrite(tcga_out, "E:/NETS/bulk&影像 TCGA外部验证/原数据/1_TCGA_Filtered_Data/TCGA_combat_adjusted.csv")
cat("校正后的 TCGA 表达矩阵已保存至：\nE:/NETS/bulk&影像 TCGA外部验证/原数据/1_TCGA_Filtered_Data/TCGA_combat_adjusted.csv\n")

# 可选：同时保存校正后的参考数据（用于后续验证）
ref_out <- as.data.frame(expr_ref_combat)
ref_out <- cbind(gene = rownames(ref_out), ref_out)
fwrite(ref_out, "E:/NETS/bulk rna数据/原数据/1_CGGA_118_Filtered/CGGA_118_combat_adjusted.csv")
cat("校正后的参考数据已保存至同目录。\n")

# ========================= 9. 加载弹性网络结果，提取基因系数 ==============
load("E:/NETS/elasticnet_cox_results.RData")   # 加载RData，内部应包含系数信息

# 查看加载了哪些对象（调试用）
cat("RData中加载的对象：", ls(), "\n")

# ---- 智能识别基因系数对象 ----
# 常见变量名：beta_final, coef, coefficients, final_genes, genes
if (exists("beta_final") && exists("final_genes")) {
  # 若beta_final为命名向量或数据框，直接提取
  if (is.vector(beta_final) && !is.null(names(beta_final))) {
    coef_vec <- beta_final
    gene_list <- names(coef_vec)
  } else if (is.data.frame(beta_final)) {
    # 假设第一列为基因，第二列为系数
    coef_vec <- setNames(beta_final[, 2], beta_final[, 1])
    gene_list <- rownames(beta_final)  # 或第一列
  } else {
    stop("无法解析 beta_final 对象，请手动指定基因和系数。")
  }
} else if (exists("coef_df")) {
  # 若存在数据框 coef_df
  coef_vec <- setNames(coef_df[, 2], coef_df[, 1])
  gene_list <- coef_df[, 1]
} else {
  # 若无法自动识别，则提示用户手动输入（这里使用已知的6个基因）
  warning("未找到标准变量名，将使用默认6个基因（IFI30, DCN, DUSP1, MS4A6A, CRTAP, TAGLN）及固定系数。")
  gene_list <- c("IFI30", "DCN", "DUSP1", "MS4A6A", "CRTAP", "TAGLN")
  # 手动指定系数（请根据实际结果修改！）
  coef_vec <- c(0.12, -0.08, 0.15, 0.10, -0.05, 0.20)  
  names(coef_vec) <- gene_list
}

# 确保只保留基因列表中的基因，并按照基因顺序排列系数
my_genes <- intersect(gene_list, rownames(expr_tcga_combat))   # 与校正后矩阵交集
if (length(my_genes) == 0) {
  stop("校正后的表达矩阵中未找到任何指定基因，请检查基因名是否匹配。")
}
# 提取对应系数（按my_genes顺序）
beta_vec <- coef_vec[my_genes]
cat("用于风险评分的基因及系数：\n")
print(beta_vec)

# ========================= 10. 计算TCGA样本的风险评分 ====================
# 从校正后的表达矩阵提取基因（行为基因，列为样本）
expr_sig <- expr_tcga_combat[my_genes, , drop = FALSE]   # 基因 × 样本

# 转置为样本 × 基因，便于矩阵乘法
expr_sig_t <- t(expr_sig)   # 样本 × 基因

# 计算原始风险评分：每个样本 = sum(表达值 * 系数)
risk_raw <- as.numeric(as.matrix(expr_sig_t) %*% beta_vec)
names(risk_raw) <- colnames(expr_tcga_combat)   # 样本ID

cat("风险评分计算完成，样本数：", length(risk_raw), "\n")

# ========================= 11. 读取临床数据并合并评分 ====================
clin_file <- "E:/NETS/bulk&影像 TCGA外部验证/原数据/1_TCGA_Filtered_Data/clinical_filtered_final.csv"
clinical <- fread(clin_file, data.table = FALSE)
cat("临床数据维度（行×列）：", dim(clinical), "\n")

# 查看临床数据列名，找到样本ID列（通常为 "sample" 或 "bcr_patient_barcode"）
colnames_clin <- colnames(clinical)
cat("临床数据列名：", colnames_clin, "\n")

# 尝试多种常见样本ID列名
id_col <- intersect(colnames_clin, c("sample", "Sample", "SAMPLE", "bcr_patient_barcode", "patient", "id"))
if (length(id_col) == 0) {
  # 如果找不到，假设第一列为样本ID
  id_col <- colnames_clin[1]
  cat("未找到标准样本ID列，将使用第一列作为样本ID：", id_col, "\n")
} else {
  id_col <- id_col[1]
  cat("使用样本ID列：", id_col, "\n")
}

# 将风险评分转换为数据框，并确保样本ID与临床数据匹配
risk_df <- data.frame(sample_id = names(risk_raw), risk_score = risk_raw, stringsAsFactors = FALSE)

# 合并：左连接，保留所有临床样本（若某些样本无评分，则NA）
clinical_merged <- merge(clinical, risk_df, by.x = id_col, by.y = "sample_id", all.x = TRUE)

# 检查合并后的维度
cat("合并后临床数据维度：", dim(clinical_merged), "\n")

# ========================= 12. 保存结果 ================================
out_file <- "E:/NETS/bulk&影像 TCGA外部验证/原数据/1_TCGA_Filtered_Data/TCGA_clinical_with_risk.csv"
fwrite(clinical_merged, out_file)
cat("包含风险评分的临床数据已保存至：\n", out_file, "\n")

# 可选：单独保存风险评分表（样本ID + 评分）
risk_out_file <- "E:/NETS/bulk&影像 TCGA外部验证/原数据/1_TCGA_Filtered_Data/TCGA_risk_scores.csv"
fwrite(risk_df, risk_out_file)
cat("风险评分单独保存至：\n", risk_out_file, "\n")
