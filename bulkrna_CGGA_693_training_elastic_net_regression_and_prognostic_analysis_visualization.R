# ========================= 0. 加载包 ========================================
library(WGCNA)
library(tidyverse)
library(data.table)
library(randomForest)
library(caret)
library(glmnet)
library(ggplot2)
library(survminer)   # 可选，用于生存分析

setwd("E:/NETS/")    # 请修改为实际路径
out_dir <- "E:/NETS/bulk rna数据/结果/bulk rna数据 4-14"

# ========================= 1. 定义标准 ssGSEA 函数 ==========================
ssgsea_standard <- function(expr, gene_set) {
  # expr: 样本 × 基因矩阵（已log转换，无NA）
  # gene_set: 基因名向量
  gene_set <- intersect(gene_set, colnames(expr))
  if (length(gene_set) < 3) {
    stop("基因集大小 < 3，无法计算 ssGSEA。请检查基因名匹配情况。")
  }
  n_samples <- nrow(expr)
  scores <- numeric(n_samples)
  for (i in 1:n_samples) {
    x <- expr[i, ]
    names(x) <- colnames(expr)
    x <- x[!is.na(x) & !is.nan(x)]
    if (length(x) == 0) {
      scores[i] <- NA
      next
    }
    x_sorted <- sort(x, decreasing = TRUE)
    is_in_set <- names(x_sorted) %in% gene_set
    N_hit <- sum(is_in_set)
    N_miss <- length(is_in_set) - N_hit
    if (N_hit == 0 || N_miss == 0) {
      scores[i] <- NA
      next
    }
    P_hit <- cumsum(is_in_set) / N_hit
    P_miss <- cumsum(!is_in_set) / N_miss
    diff <- P_hit - P_miss
    ES <- max(diff)   # 富集分数取最大偏离
    scores[i] <- ES
  }
  if (sum(!is.na(scores)) > 1) {
    scores <- scale(scores)[, 1]
  }
  return(scores)
}

# ========================= 2. 读取数据 ======================================
# 2.1 差异基因列表（单细胞上调）
deg <- fread("单细胞数据/结果/DE_up_genes_filtered.csv", data.table = FALSE)
target_genes <- deg$gene
cat("原始目标基因数量:", length(target_genes), "\n")

# 2.2 Bulk 表达矩阵（FPKM）
fpkm_raw <- fread("bulk rna数据/原数据/1_CGGA_118_Filtered/CGGA.mRNAseq_693_filtered.txt",
                  data.table = FALSE)
rownames(fpkm_raw) <- fpkm_raw$Gene_Name
fpkm_raw <- fpkm_raw[, -1]                         # 移除 Gene_Name 列
colnames(fpkm_raw) <- as.character(colnames(fpkm_raw))

# 转换为数值矩阵
fpkm <- as.matrix(fpkm_raw)
mode(fpkm) <- "numeric"
fpkm[is.na(fpkm)] <- 0
cat("全基因表达矩阵维度（基因 × 样本）:", dim(fpkm), "\n")

# ========================= 3. 基因名匹配检查 ================================
common_genes <- intersect(target_genes, rownames(fpkm))
cat("实际匹配到的目标基因数量:", length(common_genes), "\n")
if (length(common_genes) == 0) {
  stop("错误：没有匹配到任何目标基因！请检查基因名格式（大小写、版本号、ENSG vs SYMBOL）。")
}

# ========================= 4. 准备用于 ssGSEA 的全基因矩阵 ==================
# 注意：ssGSEA 需要全基因背景（包含所有基因），且为样本×基因，已 log2 转换
expr_all_for_gsva <- t(fpkm)                        # 样本 × 基因
expr_all_for_gsva <- as.matrix(expr_all_for_gsva)
mode(expr_all_for_gsva) <- "numeric"
# 缺失值按列（基因）中位数填充（一般fpkm已处理为0，但保留稳健操作）
expr_all_for_gsva <- apply(expr_all_for_gsva, 2, function(x) {
  if (any(is.na(x))) {
    med <- median(x, na.rm = TRUE)
    x[is.na(x)] <- med
  }
  x
})
expr_all_for_gsva <- log2(expr_all_for_gsva + 1)    # log2 转换
colnames(expr_all_for_gsva) <- rownames(fpkm)       # 基因名
rownames(expr_all_for_gsva) <- colnames(fpkm)       # 样本名

# 计算原始基因集的 ssGSEA 评分（目标表型）
score_original <- ssgsea_standard(expr_all_for_gsva, common_genes)
names(score_original) <- rownames(expr_all_for_gsva)
cat("有效原始评分数量:", sum(!is.na(score_original)), "\n")
summary(score_original)

if (sum(!is.na(score_original)) == 0) {
  stop("原始评分为全 NA，请检查 ssGSEA 函数或基因集是否为空。")
}

# 输出原始ssGSEA评分表格
ssgsea_score_df <- data.frame(
  Sample = names(score_original),
  SSGSEA_Score = score_original
)
write.csv(ssgsea_score_df, "ssgsea_score_original.csv", row.names = FALSE)
cat("原始ssGSEA评分表格已保存至 ssgsea_score_original.csv\n")

# ==================== 在计算高变基因之前，去除噪音基因（模仿单细胞过滤） ====================
# 定义需要过滤的基因名称模式（与单细胞分析保持一致）
bad_prefix <- c("^RPS", "^RPL", "^MT-", "^KRT", "^HSP", "^EEF", "^TUB", "^ACT",
                "^IFIT", "^ISG", "^OAS", "^HLA", "^S100")

# 找出所有匹配上述模式的基因（基于 rownames(fpkm)）
bad_genes <- unique(grep(paste(bad_prefix, collapse="|"), rownames(fpkm), value=TRUE, ignore.case=FALSE))
cat("待过滤的噪音基因数量:", length(bad_genes), "\n")

# 从表达矩阵中剔除这些基因
fpkm_filtered <- fpkm[!rownames(fpkm) %in% bad_genes, ]
cat("过滤后剩余基因数量:", nrow(fpkm_filtered), "\n")

# 提示：如果您后续还需要保留原始 fpkm 用于其他分析（如 ssGSEA 的全基因背景），
# 可以继续使用原始 fpkm 计算 ssGSEA 评分（因为 ssGSEA 需要全基因）。
# 但高变基因提取、弹性网络、风险模型构建全部使用 fpkm_filtered。

# 将后续所有使用 fpkm 的地方替换为 fpkm_filtered
# 例如：
#   - 计算高变基因时：gene_mad <- apply(fpkm_filtered, 1, mad, na.rm=TRUE)
#   - 构建 expr_high 时：expr_high <- t(fpkm_filtered[top_genes, , drop=FALSE])

# ========================= 5. 提取高变基因（用于特征筛选）===================
# 使用 median absolute deviation (MAD) 筛选 top 5000
gene_mad <- apply(fpkm_filtered, 1, mad, na.rm = TRUE)
top_genes <- names(sort(gene_mad, decreasing = TRUE)[1:min(5000, nrow(fpkm))])
cat("高变基因数量:", length(top_genes), "\n")

# 直接取子集并转置（fpkm 已无 NA，且为基因×样本）
expr_high <- t(fpkm_filtered[top_genes, , drop = FALSE])   # 样本 × 基因
expr_high <- as.matrix(expr_high)
mode(expr_high) <- "numeric"
expr_high <- log2(expr_high + 1)                  # log2 转换

# 此时行列名已自动继承：行名 = colnames(fpkm)，列名 = top_genes
cat("高变基因表达矩阵维度（样本 × 基因）:", dim(expr_high), "\n")

# 匹配样本（确保评分和表达矩阵样本一致）
common_samples <- intersect(rownames(expr_high), names(score_original))
X <- expr_high[common_samples, , drop = FALSE]
y <- score_original[common_samples]
cat("最终用于弹性网络的样本数:", length(y), "，特征数:", ncol(X), "\n")

# ========================= 1. 弹性网络筛选正系数基因 ========================
set.seed(123)
cv_enet <- cv.glmnet(x = X, y = y, alpha = 0.5, family = "gaussian",
                     nfolds = 5, type.measure = "mse")
lambda_min <- cv_enet$lambda.min
coef_min <- coef(cv_enet, s = lambda_min)
coef_vec <- as.matrix(coef_min)[-1, 1]
names(coef_vec) <- rownames(coef_min)[-1]

# 正系数基因（与 NETs 评分正相关）
pos_genes <- names(coef_vec[coef_vec > 0])
cat("弹性网络正系数基因数量:", length(pos_genes), "\n")

# 按系数绝对值排序
pos_genes_sorted <- names(sort(abs(coef_vec[pos_genes]), decreasing = TRUE))

# 保存弹性网络正系数基因（按系数绝对值排序）
writeLines(pos_genes_sorted, file.path(out_dir, "elasticnet_positive_genes.txt"))
cat("正系数基因已保存至:", file.path(out_dir, "elasticnet_positive_genes.txt"), "\n")

# ========================= 2. Bootstrap 稳定性选择 ==========================
set.seed(123)
n_bootstrap <- 100
freq <- rep(0, ncol(X))
names(freq) <- colnames(X)

for (i in 1:n_bootstrap) {
  idx <- sample(1:nrow(X), replace = TRUE)
  cv_i <- cv.glmnet(X[idx, ], y[idx], alpha = 0.5, nfolds = 5)
  coef_i <- coef(cv_i, s = cv_i$lambda.min)
  sel <- rownames(coef_i)[which(coef_i != 0)][-1]
  freq[sel] <- freq[sel] + 1
}

stable_genes <- names(freq[freq >= 60])   # 60% 稳定性阈值
cat("Bootstrap 稳定性基因（频率≥60%）: ", length(stable_genes), "\n")

# 保存 Bootstrap 稳定性基因（频率≥60%）
writeLines(stable_genes, file.path(out_dir, "bootstrap_stable_genes.txt"))
cat("稳定性基因已保存至:", file.path(out_dir, "bootstrap_stable_genes.txt"), "\n")

# 取正系数基因与稳定基因的交集（最终签名）
final_genes <- intersect(pos_genes_sorted, stable_genes)
cat("最终签名基因数量:", length(final_genes), "\n")
writeLines(final_genes, "elasticnet_stable_signature.txt")

# ========================= 3. 计算风险评分 ==================================
# 使用最终签名基因在原始表达矩阵（log2 FPKM+1）中的表达量加权求和
# 注意：这里使用全基因矩阵 expr_all_for_gsva（样本×基因），只取 final_genes
expr_sig <- expr_all_for_gsva[, final_genes, drop = FALSE]
beta_final <- coef_vec[final_genes]   # 从弹性网络系数中提取对应系数
risk_score <- as.matrix(expr_sig) %*% beta_final
risk_score <- scale(risk_score)[,1]   # 标准化，便于比较

# ==================== 风险评分与原始 ssGSEA 评分的相关性 ====================
# 确保输出文件夹存在（若之前已定义 out_dir 则直接使用）
if (!exists("out_dir")) {
  out_dir <- "E:/NETS/bulk rna数据/结果/bulk rna数据 4-14"
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
}

# 共同样本
common_cor <- intersect(names(score_original), names(risk_score))
if (length(common_cor) > 0) {
  df_cor <- data.frame(
    Sample = common_cor,
    ssGSEA = score_original[common_cor],
    RiskScore = risk_score[common_cor]
  )
  
  # Pearson 相关性检验
  cor_test <- cor.test(df_cor$ssGSEA, df_cor$RiskScore, method = "pearson")
  r_val <- round(cor_test$estimate, 3)
  p_val <- cor_test$p.value
  # 若 p 值极小，用科学计数法显示
  p_label <- ifelse(p_val < 0.001, format(p_val, scientific = TRUE, digits = 3),
                    round(p_val, 4))
  
  # 绘图：散点 + 线性回归线 + 95% 置信带（渐近线方差）
  library(ggplot2)
  p_cor <- ggplot(df_cor, aes(x = ssGSEA, y = RiskScore)) +
    geom_point(alpha = 0.6, size = 2, color = "#2c7bb6") +
    geom_smooth(method = "lm", se = TRUE, color = "#d7191c", fill = "grey70", 
                level = 0.95) +
    labs(x = "Original ssGSEA Score", 
         y = "Risk Score (scaled)",
         title = "Correlation between Risk Score and ssGSEA Score") +
    annotate("text", 
             x = min(df_cor$ssGSEA, na.rm = TRUE), 
             y = max(df_cor$RiskScore, na.rm = TRUE),
             label = paste0("Pearson's r = ", r_val, "\np = ", p_label),
             hjust = 0, vjust = 1, size = 4.5, color = "black") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # 显示图形
  print(p_cor)
  
  # 保存为 PDF 和 PNG
  ggsave(file.path(out_dir, "RiskScore_vs_ssGSEA_correlation.pdf"), 
         p_cor, width = 6, height = 5)
  ggsave(file.path(out_dir, "RiskScore_vs_ssGSEA_correlation.png"), 
         p_cor, width = 6, height = 5, dpi = 300)
  
  cat("相关性分析完成：r =", r_val, ", p =", p_val, "\n")
  cat("图形已保存至:", file.path(out_dir, "RiskScore_vs_ssGSEA_correlation.pdf/png"), "\n")
} else {
  cat("警告：风险评分与 ssGSEA 评分无重叠样本，无法计算相关性。\n")
}


# ==================== 风险评分高低组的差异分析 + GO BP 富集 ====================
# 1. 分组（中位数）
risk_group <- ifelse(risk_score > median(risk_score), "High", "Low")
table(risk_group)

# 2. 准备表达矩阵（基因 × 样本，log2(FPKM+1)）
expr_for_diff <- t(expr_all_for_gsva)  # 基因 × 样本
expr_for_diff <- expr_for_diff[, names(risk_score)]  # 对齐样本

# 3. 差异分析（limma 简洁版）
library(limma)
group <- factor(risk_group, levels = c("Low", "High"))
design <- model.matrix(~ group)
fit <- lmFit(expr_for_diff, design)
fit <- eBayes(fit)
res <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH")

# 筛选上调基因（High vs Low 表达升高，logFC > 0.5 且 adj.P.Val < 0.05）
up_genes <- res[res$logFC > 1 & res$adj.P.Val < 0.05, ]
cat("上调基因数量:", nrow(up_genes), "\n")
if (nrow(up_genes) == 0) {
  cat("没有显著上调基因，尝试降低阈值 logFC > 0.3\n")
  up_genes <- res[res$logFC > 0.3 & res$adj.P.Val < 0.05, ]
}

# 4. GO BP 富集（仅 Biological Process）
if (exists("up_genes") && nrow(up_genes) > 0) {
  library(clusterProfiler)
  library(org.Hs.eg.db)   # 人类注释包，如果是其他物种需修改
  gene_list <- rownames(up_genes)
  
  # 转换基因名为 Entrez ID（可选，但推荐）
  entrez_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  if (nrow(entrez_ids) > 0) {
    go_bp <- enrichGO(gene          = entrez_ids$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID",
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.2,
                      readable      = TRUE)
    
    if (!is.null(go_bp) && nrow(go_bp) > 0) {
      # 保存结果表格
      write.csv(as.data.frame(go_bp), file.path(out_dir, "GO_BP_high_vs_low.csv"), row.names = FALSE)
      
      # 绘制气泡图（前10个条目）
      pdf(file.path(out_dir, "GO_BP_bubble.pdf"), width = 8, height = 6)
      print(dotplot(go_bp, showCategory = 10, title = "GO BP Enrichment in High Risk Group"))
      dev.off()
      
      cat("GO BP 富集完成，结果保存至:", file.path(out_dir, "GO_BP_high_vs_low.csv"), "\n")
    } else {
      cat("没有显著富集的 GO BP 条目。\n")
    }
  } else {
    cat("基因名无法转换为 Entrez ID，跳过富集。\n")
  }
} else {
  cat("没有显著上调基因，跳过 GO BP 富集。\n")
}
# 预排序GSEA：基于所有基因的logFC排序，检验NETs基因集是否在高风险组富集
library(fgsea)
# gene_list: 所有基因的logFC（High vs Low），从高到低排序
gene_list <- res$logFC
names(gene_list) <- rownames(res)
gene_list <- sort(gene_list, decreasing = TRUE)

# 您的原始NETs基因集（common_genes）
net_set <- list(NETs = intersect(common_genes, names(gene_list)))

# 运行GSEA
fgsea_res <- fgsea(net_set, gene_list, minSize = 5, maxSize = 500)
print(fgsea_res)
# ========================= 4. Cox 回归分析 ==================================
library(survival)
# 读取临床数据
clinical <- fread("bulk rna数据/原数据/1_CGGA_118_Filtered/CGGA.mRNAseq_693_clinical_filtered.txt", 
                  data.table = FALSE)

# 重命名列，以匹配后续 Surv 函数的标准名称
colnames(clinical)[colnames(clinical) == "CGGA_ID"] <- "sample"
colnames(clinical)[colnames(clinical) == "OS"] <- "OS.time"
colnames(clinical)[colnames(clinical) == "Censor (alive=0; dead=1)"] <- "OS.status"

# 确保 OS.status 为数值型（0/1）
clinical$OS.status <- as.numeric(clinical$OS.status)
# 创建风险评分数据框（包含样本名）
risk_df <- data.frame(sample = names(risk_score),
                      risk = as.numeric(risk_score),
                      stringsAsFactors = FALSE)

# 合并（假设 clinical 中有 'sample' 列）
merged <- merge(clinical, risk_df, by = "sample", all = FALSE)  # 或 all.x = TRUE

# 检查合并后的数据
head(merged)
# 检查样本匹配：clinical 中的 sample 是否与 rownames(expr_sig) 完全一致
common_samples <- intersect(clinical$sample, rownames(expr_sig))
cat("临床数据中的样本数:", nrow(clinical), "\n")
cat("风险评分中的样本数:", nrow(expr_sig), "\n")
cat("共同样本数:", length(common_samples), "\n")

if (length(common_samples) == 0) {
  stop("临床样本与表达样本完全不一致，请检查 sample ID 格式（例如大小写、前缀等）。")
}

# 只保留共同样本
clinical <- clinical[clinical$sample %in% common_samples, ]
expr_sig <- expr_sig[common_samples, , drop = FALSE]
risk_score <- risk_score[common_samples]

# 合并数据（确保顺序一致）
merged <- data.frame(sample = rownames(expr_sig), 
                     risk = risk_score,
                     clinical[match(rownames(expr_sig), clinical$sample), 
                              c("OS.time", "OS.status", "Age", "Grade", "IDH_mutation_status")])
# 注意：临床文件中的列名可能是 "Age", "Grade", "IDH_mutation_status" 等，请根据实际列名修改。
# 若没有 Age 等变量，多因素 Cox 中可移除。

# 转换分级为数值（如果需要）
# merged$Grade <- as.numeric(factor(merged$Grade, levels = c("WHO II", "WHO III", "WHO IV")))
# 假设你已经加载了临床数据 clinical，包含样本名、OS.time、OS.status、age、grade、IDH等
# clinical <- fread("CGGA_clinical.csv", data.table = FALSE)
# 合并风险评分与临床数据
merged <- data.frame(sample = rownames(expr_sig), risk = risk_score)
merged <- merge(merged, clinical, by = "sample")

# 4.1 单因素 Cox
cox_uni <- coxph(Surv(OS.time, OS.status) ~ risk, data = merged)
summary(cox_uni)

# 4.2 多因素 Cox（根据实际临床变量调整）
# 确保分类变量格式正确
merged$Grade <- factor(merged$Grade, levels = c("WHO II", "WHO III", "WHO IV"))
merged$IDH_mutation_status <- factor(merged$IDH_mutation_status, levels = c("Wildtype", "Mutant"))

# 拟合模型（使用正确的列名）
cox_multi <- coxph(Surv(OS.time, OS.status) ~ risk + Age + Grade + IDH_mutation_status, 
                   data = merged)
summary(cox_multi)

# 比例风险假设检验
ph_test <- cox.zph(cox_multi)
print(ph_test)

# 若发现某个变量违反 PH 假设（p < 0.05），可考虑分层或时变系数
if (any(ph_test$table[, "p"] < 0.05)) {
  warning("某些变量违反了比例风险假设，建议检查 Schoenfeld 残差图或考虑分层分析。")
}

# 4.3 比例风险假设检验
ph_test <- cox.zph(cox_multi)
print(ph_test)
# 若有变量违反 PH 假设（p < 0.05），可考虑分层或时变系数，这里仅提示
if (any(ph_test$table[, "p"] < 0.05)) {
  warning("某些变量违反了比例风险假设，请检查 cox.zph 结果。")
}

# ==================== Cox 回归可视化 ====================
# 确保输出目录存在
if (!exists("out_dir")) {
  out_dir <- "E:/NETS/bulk rna数据/结果/bulk rna数据 4-14"
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
}

# 1. 多因素 Cox 森林图 (Forest plot)
library(survminer)
# 绘制森林图（会自动显示HR和置信区间）
forest_multi <- ggforest(cox_multi, 
                         data = merged,
                         main = "Hazard ratios of risk score and clinical variables",
                         cpositions = c(0.02, 0.22, 0.4),  # 调整列位置
                         fontsize = 1.0,
                         refLabel = "reference",
                         noDigits = 3)

# 保存森林图（PDF和PNG）
pdf(file.path(out_dir, "Cox_multi_forest.pdf"), width = 8, height = 6)
print(forest_multi)
dev.off()

png(file.path(out_dir, "Cox_multi_forest.png"), width = 2400, height = 1800, res = 300)
print(forest_multi)
dev.off()

cat("多因素Cox森林图已保存至:", out_dir, "\n")

# 2. （可选）单因素Cox森林图 - 仅风险评分
forest_uni <- ggforest(cox_uni, 
                       data = merged,
                       main = "Univariate Cox: Risk Score",
                       cpositions = c(0.1, 0.3, 0.5),
                       fontsize = 1.0)
pdf(file.path(out_dir, "Cox_uni_forest.pdf"), width = 6, height = 3)
print(forest_uni)
dev.off()

# 3. 列线图 (Nomogram) - 预测1年、3年、5年生存概率
# 需要安装 rms 包
if (!require("rms")) install.packages("rms")
library(rms)
# 对多因素Cox模型进行转化，以便绘制列线图
dd <- datadist(merged)
options(datadist = "dd")

# 重新拟合Cox模型（使用rms的cph函数）
cox_nomo <- cph(Surv(OS.time, OS.status) ~ risk + Age + Grade + IDH_mutation_status,
                data = merged, surv = TRUE, x = TRUE, y = TRUE)

# 设定预测时间点（单位：天，假设OS.time是天）
# 您可以根据实际情况修改：1年=365，3年=1095，5年=1825
time_points <- c(365, 1095, 1825)

# 绘制列线图
pdf(file.path(out_dir, "Nomogram.pdf"), width = 10, height = 7)
plot(nomogram(cox_nomo, 
              fun = list(function(x) surv(365, x),
                         function(x) surv(1095, x),
                         function(x) surv(1825, x)),
              funlabel = c("1-year survival", "3-year survival", "5-year survival"),
              lp = FALSE,
              conf.int = FALSE))
dev.off()

cat("列线图已保存至:", file.path(out_dir, "Nomogram.pdf"), "\n")

# 4. 校准曲线 (Calibration) - 用于3年生存概率
# 需要加载 rms 包
# 计算3年生存概率的校准曲线
cal_3yr <- calibrate(cox_nomo, u = 1095, cmethod = "KM", m = 50, B = 200)
pdf(file.path(out_dir, "Calibration_3yr.pdf"), width = 6, height = 6)
plot(cal_3yr, xlab = "Predicted 3-year survival", ylab = "Observed 3-year survival",
     main = "Calibration curve for 3-year survival")
abline(0,1, col = "red", lty = 2)
dev.off()

cat("3年校准曲线已保存至:", file.path(out_dir, "Calibration_3yr.pdf"), "\n")

# ========================= 7. 合并风险评分与临床信息并保存 =========================
# 确保 merged 数据框中已有 risk (原始评分) 和 risk_norm (归一化评分)
if (!"risk_norm" %in% colnames(merged)) {
  risk_min <- min(merged$risk, na.rm = TRUE)
  risk_max <- max(merged$risk, na.rm = TRUE)
  merged$risk_norm <- (merged$risk - risk_min) / (risk_max - risk_min)
  cat("已计算 risk_norm (0-1 归一化)\n")
}

# 创建风险分组 (基于原始风险评分的中位数)
merged$risk_group <- ifelse(merged$risk > median(merged$risk, na.rm = TRUE), "High", "Low")

# 选择要保存的列：样本ID、风险原始值、风险归一化值、风险分组、以及重要的临床变量
# 临床变量根据您的实际列名调整，常见的有 OS.time, OS.status, Age, Grade, IDH_mutation_status
cols_to_save <- c("sample", "risk", "risk_norm", "risk_group", 
                  "OS.time", "OS.status", "Age", "Grade", "IDH_mutation_status")
# 只保留存在的列
cols_to_save <- intersect(cols_to_save, colnames(merged))

# 创建输出数据框
output_df <- merged[, cols_to_save, drop = FALSE]

# 保存为 CSV 文件
output_file <- file.path(out_dir, "NETS_signature_with_clinical.csv")
write.csv(output_df, output_file, row.names = FALSE)
cat("风险评分及临床信息已保存至:", output_file, "\n")

# ========================= 5. Kaplan-Meier 曲线 =============================
# 根据风险评分中位数分组
merged$risk_group <- ifelse(merged$risk > median(merged$risk, na.rm = TRUE), "High", "Low")
merged$risk_group <- factor(merged$risk_group, levels = c("Low", "High"))

# 拟合 KM 曲线
fit_km <- survfit(Surv(OS.time, OS.status) ~ risk_group, data = merged)

# 创建纯 KM 曲线图（无风险表、无 p 值、无置信区间）
p <- ggsurvplot(fit_km,
                data = merged,
                pval = TRUE,          # 不显示 p 值文本
                pval.method = TRUE,   # 不显示检验方法
                conf.int = FALSE,      # 不显示置信区间带
                risk.table = FALSE,    # 不显示风险表
                palette = c("orange", "blue"),
                xlab = "Time (days)",
                ylab = "Overall Survival Probability",
                title = "Kaplan-Meier Curve by Risk Group",
                legend.title = "Risk Group",
                legend.labs = c("Low Risk", "High Risk"))

# 只保存曲线图（不包含风险表）
ggsave("km_curve_clean.pdf", plot = p$plot, width = 6, height = 5)
ggsave("km_curve_clean.png", plot = p$plot, width = 6, height = 5, dpi = 300)
# ========================= 6. 保存结果 ======================================
save(cv_enet, pos_genes_sorted, stable_genes, final_genes, beta_final,
     risk_score, cox_uni, cox_multi, ph_test, file = "elasticnet_cox_results.RData")
cat("分析完成！最终签名基因已保存至 elasticnet_stable_signature.txt\n")

# ========================= timeROC (时间依赖ROC) =============================
# 安装并加载 timeROC
if (!require(timeROC)) install.packages("timeROC")
library(timeROC)

# 假设 merged_roc 已正确创建（包含 OS.time, OS.status, risk）
# 如果 merged_roc 未定义，先运行：
merged_roc <- na.omit(merged[, c("OS.time", "OS.status", "risk")])

# 定义时间点（1年、3年、5年，单位与 OS.time 一致，假设是天）
time_points <- c(365, 1095, 1825)

# 计算时间依赖ROC
roc_obj <- timeROC(T = merged_roc$OS.time,
                   delta = merged_roc$OS.status,
                   marker = merged_roc$risk,
                   cause = 1,
                   times = time_points,
                   iid = FALSE)   # iid = FALSE 避免小样本报错

# 输出AUC
print(round(roc_obj$AUC, 3))

# 绘制三条ROC曲线（叠加）
# 定义颜色
colors_roc <- c("#E41A1C", "#377EB8", "#4DAF4A")
line_width <- 2.5

# 建立空白画布 + 网格 + 对角线
plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
     xlab = "1 - Specificity", ylab = "Sensitivity",
     main = "", axes = FALSE)
box()
grid(col = "lightgray", lty = 1, lwd = 0.5)
axis(1, at = seq(0, 1, 0.2), cex.axis = 0.9)
axis(2, at = seq(0, 1, 0.2), las = 1, cex.axis = 0.9)
abline(a = 0, b = 1, lty = 2, col = "gray40", lwd = 1.2)

# 手动绘制三条 ROC 曲线（因为 TP/FP 是矩阵，用 [,i] 取列）
for (i in 1:3) {
  tpr <- roc_obj$TP[, i]   # 第 i 列，灵敏度
  fpr <- roc_obj$FP[, i]   # 第 i 列，1-特异性
  
  # 去除可能的 NA（边界情况）
  valid <- !is.na(tpr) & !is.na(fpr)
  tpr <- tpr[valid]
  fpr <- fpr[valid]
  
  # 按 fpr 排序后连线
  ord <- order(fpr)
  lines(fpr[ord], tpr[ord], col = colors_roc[i], lwd = line_width)
}

# 添加标题
title(main = "Time-dependent ROC Curves for Risk Score", line = 1.5)

# 添加图例（AUC 值从 roc_obj$AUC 获取）
auc_vals <- round(roc_obj$AUC, 3)
legend("bottomright", 
       legend = paste0(c("1-year", "3-year", "5-year"), " AUC = ", auc_vals),
       col = colors_roc, lty = 1, lwd = line_width,
       bty = "o", bg = "white", cex = 0.85, inset = c(0.02, 0.02))


# 保存图片
pdf("timeROC_optimized.pdf", width = 6.5, height = 6)

# 以下代码与您上面手动绘制的完全一致
colors_roc <- c("#E41A1C", "#377EB8", "#4DAF4A")
line_width <- 2.5

plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
     xlab = "1 - Specificity", ylab = "Sensitivity",
     main = "", axes = FALSE)
box()
grid(col = "lightgray", lty = 1, lwd = 0.5)
axis(1, at = seq(0, 1, 0.2), cex.axis = 0.9)
axis(2, at = seq(0, 1, 0.2), las = 1, cex.axis = 0.9)
abline(a = 0, b = 1, lty = 2, col = "gray40", lwd = 1.2)

for (i in 1:3) {
  tpr <- roc_obj$TP[, i]
  fpr <- roc_obj$FP[, i]
  valid <- !is.na(tpr) & !is.na(fpr)
  tpr <- tpr[valid]
  fpr <- fpr[valid]
  ord <- order(fpr)
  lines(fpr[ord], tpr[ord], col = colors_roc[i], lwd = line_width)
}

title(main = "Time-dependent ROC Curves for Risk Score", line = 1.5)

auc_vals <- round(roc_obj$AUC, 3)
legend("bottomright", 
       legend = paste0(c("1-year", "3-year", "5-year"), " AUC = ", auc_vals),
       col = colors_roc, lty = 1, lwd = line_width,
       bty = "o", bg = "white", cex = 0.85, inset = c(0.02, 0.02))

dev.off()
getwd()

# ========== 使用原始风险评分的稳定 DCA（1年/3年/5年） ==========
library(survival)
library(rmda)

# 清洗数据（省略，沿用 merged_clean）
merged_clean <- merged[!is.na(merged$OS.time) & !is.na(merged$OS.status), ]
merged_clean$OS.status <- as.numeric(merged_clean$OS.status)
merged_clean <- merged_clean[merged_clean$OS.status %in% c(0,1), ]

# Cox 模型（假设已拟合好 cox_multi）
if (!exists("cox_multi")) {
  cox_multi <- coxph(Surv(OS.time, OS.status) ~ risk + Age + Grade + IDH_mutation_status,
                     data = merged_clean)
}

# 选择预测时间点（1年=365，3年=1095，5年=1825）
t_pred <- 1095   # 3年，可根据需要修改

# ---------- 1. 预测每个样本的生存概率，并转为风险概率（可选，但这里直接用原始risk）----------
# 其实我们不需要预测概率，直接使用 merged_clean$risk 即可。
# 但为了与您的原始数据一致，假设 risk 列已经在 merged_clean 中。
df <- merged_clean[, c("risk", "OS.time", "OS.status")]
colnames(df)[1] <- "risk_marker"   # 避免与 risk 关键字冲突

# ---------- 2. 定义事件状态（不进行 keep 筛选，与参考代码一致）----------
df$event <- ifelse(df$OS.time <= t_pred & df$OS.status == 1, 1, 0)

cat("样本量:", nrow(df), " 事件数:", sum(df$event), "\n")

# ---------- 3. 运行决策曲线（使用原始 risk 评分作为标记）----------
# 关键：thresholds 范围不宜太大，根据您的 risk 实际取值范围设定
# 若 risk 评分大致在 0~2 之间，可尝试 thresholds = seq(0, 1, by = 0.05)
# 但为了与标准 DCA 图一致，通常阈值概率在 0~0.5 之间，这里我们用 0~0.5
# 但如果 risk 评分并非概率，这个阈值解释为“风险评分阈值”，也合理。
dca_model <- decision_curve(event ~ risk_marker, data = df,
                            family = binomial(link = "logit"),
                            thresholds = seq(0, 0.5, by = 0.02),   # 0.5 可能仍然过高，可改为 0.4
                            confidence.intervals = FALSE)

# ---------- 4. 绘图（自动包含 All 和 None 线）----------
pdf(paste0("DCA_risk_", t_pred/365, "year.pdf"), width = 6, height = 5)
plot_decision_curve(dca_model, 
                    curve.names = paste0("Risk score (", t_pred/365, "-year)"),
                    col = "#E41A1C", lwd = 2,
                    xlim = c(0, 0.5), ylim = c(-0.05, 0.5),
                    xlab = "Threshold (risk score)", 
                    ylab = "Net benefit",
                    confidence.intervals = FALSE)
dev.off()

