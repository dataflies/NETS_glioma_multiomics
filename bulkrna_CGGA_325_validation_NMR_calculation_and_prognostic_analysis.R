library(data.table)
library(survival)
library(survminer)
library(timeROC)
library(rmda)
library(ggplot2)
library(tidyverse)

# 输出目录
out_dir <- "E:/NETS/bulk rna数据/结果/54_validation"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 加载118样本训练得到的系数（请根据实际保存的文件名修改）
# 假设之前保存了 elasticnet_cox_results.RData 或类似文件
load("E:/NETS/elasticnet_cox_results.RData")  # 应包含 beta_final 或 coef_vec

# 定义6个目标基因
my_genes <- c("IFI30", "DCN", "DUSP1", "MS4A6A", "CRTAP", "TAGLN")

# 提取系数（如果 beta_final 存在，就用它；否则用 coef_vec）
if (exists("beta_final")) {
  coef_full <- beta_final
} else if (exists("coef_vec")) {
  coef_full <- coef_vec
} else {
  stop("未找到系数对象 beta_final 或 coef_vec")
}

# 提取8个基因的系数，缺失的补0
coef_use <- coef_full[my_genes]
coef_use[is.na(coef_use)] <- 0
cat("使用的基因及系数：\n")
print(coef_use)

# ==================== risk值 ====================

# 54样本表达矩阵文件（FPKM原始值，未log）
file_54_expr <- "E:/NETS/bulk rna数据/原数据/1_CGGA_54_Filtered/CGGA.mRNAseq_325_filtered_genes.txt"

expr_54_raw <- fread(file_54_expr, data.table = FALSE)
rownames(expr_54_raw) <- expr_54_raw[[1]]   # 第一列为基因名
expr_54_raw <- expr_54_raw[, -1, drop = FALSE]
expr_54 <- as.matrix(expr_54_raw)
mode(expr_54) <- "numeric"
expr_54[is.na(expr_54)] <- 0

# log2(FPKM+1) 转换，与训练集一致
expr_54_log <- log2(expr_54 + 1)

# 检查8个基因在表达矩阵中是否存在
use_genes <- intersect(my_genes, rownames(expr_54_log))
if (length(use_genes) == 0) stop("54样本中没有任何目标基因！")
cat("54样本中实际匹配到的基因:", use_genes, "\n")
cat("未匹配到的基因:", setdiff(my_genes, use_genes), "\n")

# 提取表达子矩阵（样本×基因）
expr_sig <- t(expr_54_log[use_genes, , drop = FALSE])   # 样本 × 基因
coef_use_sub <- coef_use[use_genes]

# 计算原始风险评分（线性加权和）
risk_raw <- as.matrix(expr_sig) %*% coef_use_sub
risk_raw <- as.numeric(risk_raw)
names(risk_raw) <- rownames(expr_sig)

# 是否标准化？为了与训练集尺度一致，通常需要对验证集使用相同的标准化参数（如果训练集做了scale，验证集也要用训练集的均值和标准差）
# 但风险评分本身是一个线性组合，如果不标准化，分组时用中位数即可。我们可以在验证集内部直接使用原始评分，后续分析使用中位数分组。
# 为了结果可读，可以选择标准化（但标准化参数应来自训练集，否则会偏移）。这里推荐直接使用原始评分进行后续生存分析，因为比例风险模型只关心相对顺序。
# 我们保留原始评分 risk_raw，用于COX回归。

# 也可以将评分缩放到0-1（可选）
risk_01 <- (risk_raw - min(risk_raw)) / (max(risk_raw) - min(risk_raw))

# 保存风险评分结果
result_54 <- data.frame(Sample = names(risk_raw), RiskScore = risk_raw, RiskScore_01 = risk_01)
fwrite(result_54, file.path(out_dir, "54_samples_risk_scores.csv"))

# ==================== Cox  ====================

# 读取临床数据
clinical_54 <- fread("E:/NETS/bulk rna数据/原数据/1_CGGA_54_Filtered/CGGA.mRNAseq_325_clinical_filtered.txt", data.table = FALSE)

# 查看列名，重命名必要的列
colnames(clinical_54)
# 假设有 Sample/CGGA_ID, OS, Censor (alive=0; dead=1) 等列
# 请根据实际列名修改以下代码
if ("CGGA_ID" %in% colnames(clinical_54)) {
  clinical_54$sample <- clinical_54$CGGA_ID
} else if ("Sample" %in% colnames(clinical_54)) {
  clinical_54$sample <- clinical_54$Sample
} else {
  clinical_54$sample <- clinical_54[,1]  # 第一列为样本ID
}

# 重命名生存时间和状态列
if ("OS" %in% colnames(clinical_54)) {
  colnames(clinical_54)[colnames(clinical_54) == "OS"] <- "OS.time"
}
if ("Censor (alive=0; dead=1)" %in% colnames(clinical_54)) {
  colnames(clinical_54)[colnames(clinical_54) == "Censor (alive=0; dead=1)"] <- "OS.status"
} else if ("Censor" %in% colnames(clinical_54)) {
  colnames(clinical_54)[colnames(clinical_54) == "Censor"] <- "OS.status"
}
# 确保OS.status是数值型
clinical_54$OS.status <- as.numeric(clinical_54$OS.status)

# 合并风险评分与临床数据
merged_54 <- merge(result_54, clinical_54, by.x = "Sample", by.y = "sample", all.x = TRUE)

# 移除生存数据缺失的样本
merged_54 <- merged_54[!is.na(merged_54$OS.time) & !is.na(merged_54$OS.status), ]
cat("有效样本数:", nrow(merged_54), "\n")

# 单因素Cox
cox_uni_54 <- coxph(Surv(OS.time, OS.status) ~ RiskScore, data = merged_54)
summary(cox_uni_54)

# 多因素Cox（根据可用临床变量调整，例如 Age, Grade, IDH_mutation_status）
clinical_vars <- c("Age", "Grade", "IDH_mutation_status")
available_vars <- intersect(clinical_vars, colnames(merged_54))
if (length(available_vars) > 0) {
  # 转换分类变量为因子
  if ("Grade" %in% available_vars) {
    merged_54$Grade <- factor(merged_54$Grade, levels = c("WHO II", "WHO III", "WHO IV"))
  }
  if ("IDH_mutation_status" %in% available_vars) {
    merged_54$IDH_mutation_status <- factor(merged_54$IDH_mutation_status, levels = c("Wildtype", "Mutant"))
  }
  formula_multi <- as.formula(paste("Surv(OS.time, OS.status) ~ RiskScore +", paste(available_vars, collapse = " + ")))
  cox_multi_54 <- coxph(formula_multi, data = merged_54)
  print(summary(cox_multi_54))
} else {
  cat("无可用临床变量，只做单因素\n")
  cox_multi_54 <- NULL
}

# ------------------- 验证集 KM 曲线（无风险表、无标注） -------------------

# 1. 分组（按 RiskScore 中位数）
merged_54$risk_group <- ifelse(merged_54$RiskScore > median(merged_54$RiskScore), "High", "Low")
# 明确因子顺序（低风险在前，高风险在后）
merged_54$risk_group <- factor(merged_54$risk_group, levels = c("Low", "High"))

# 2. 拟合 KM 曲线
fit_km <- survfit(Surv(OS.time, OS.status) ~ risk_group, data = merged_54)

# 3. 创建绘图对象（关闭所有辅助元素）
p <- ggsurvplot(fit_km,
                data = merged_54,
                pval = TRUE,          # 不显示 p 值
                pval.method = TRUE,   # 不显示检验方法
                conf.int = FALSE,      # 不显示置信区间阴影
                risk.table = FALSE,    # 不显示风险表
                palette = c("orange", "blue"),   # 低风险=橙色，高风险=蓝色
                xlab = "Time (days)",
                ylab = "Overall Survival Probability",
                title = "Kaplan-Meier Curve in Validation Set (54 samples)",
                legend.title = "Risk Group",
                legend.labs = c("Low Risk", "High Risk"))

# 4. 保存为 PDF（只输出曲线图，不含风险表）
pdf("KM_curve_validation_clean.pdf", width = 6, height = 5)
print(p$plot)   # p$plot 是纯曲线图，不含风险表
dev.off()

# 5. 可选：另存为 PNG 高清图
ggsave("KM_curve_validation_clean.png", plot = p$plot, width = 8, height = 6, dpi = 300)
# ==================== 森林图 + PH检验（54验证集） ====================
library(survminer)

# 设置输出目录（如果之前没有设置，请修改为实际路径）
if (!exists("out_dir")) {
  out_dir <- "E:/NETS/bulk rna数据/结果/54_validation"
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
}

# ---------- 1. 单因素Cox森林图 ----------
forest_uni <- ggforest(cox_uni_54, 
                       data = merged_54,
                       main = "Univariate Cox: Risk Score (Validation Set)",
                       cpositions = c(0.02, 0.22, 0.4),
                       fontsize = 1.0,
                       refLabel = "reference",
                       noDigits = 3)

# 保存单因素森林图
pdf(file.path(out_dir, "Univariate_Cox_forest_54.pdf"), width = 8, height = 4)
print(forest_uni)
dev.off()
png(file.path(out_dir, "Univariate_Cox_forest_54.png"), width = 2400, height = 1200, res = 300)
print(forest_uni)
dev.off()
cat("单因素森林图已保存至:", out_dir, "\n")

# ---------- 2. 多因素Cox森林图（如果存在） ----------
if (!is.null(cox_multi_54)) {
  forest_multi <- ggforest(cox_multi_54, 
                           data = merged_54,
                           main = "Multivariate Cox (Validation Set)",
                           cpositions = c(0.02, 0.22, 0.4),
                           fontsize = 1.0,
                           refLabel = "reference",
                           noDigits = 3)
  pdf(file.path(out_dir, "Multivariate_Cox_forest_54.pdf"), width = 8, height = 6)
  print(forest_multi)
  dev.off()
  png(file.path(out_dir, "Multivariate_Cox_forest_54.png"), width = 2400, height = 1800, res = 300)
  print(forest_multi)
  dev.off()
  cat("多因素森林图已保存至:", out_dir, "\n")
}

# ---------- 3. 比例风险假设检验 (PH检验) ----------
# 单因素模型PH检验
ph_uni <- cox.zph(cox_uni_54)
cat("\n========== 单因素Cox PH检验 ==========\n")
print(ph_uni)

# 绘制Schoenfeld残差图（单因素）
pdf(file.path(out_dir, "PH_test_uni_54.pdf"), width = 6, height = 5)
plot(ph_uni, main = "Schoenfeld Residuals - Univariate Cox")
abline(h = 0, col = "red", lty = 2)
dev.off()

# 多因素模型PH检验（如果存在）
if (!is.null(cox_multi_54)) {
  ph_multi <- cox.zph(cox_multi_54)
  cat("\n========== 多因素Cox PH检验 ==========\n")
  print(ph_multi)
  
  # 绘制Schoenfeld残差图（多因素，按变量分面）
  pdf(file.path(out_dir, "PH_test_multi_54.pdf"), width = 8, height = 6)
  plot(ph_multi, main = "Schoenfeld Residuals - Multivariate Cox")
  dev.off()
}

# 输出PH检验的p值汇总（便于快速查看）
if (!is.null(cox_multi_54)) {
  ph_p_uni <- round(ph_uni$table[1, "p"], 4)
  ph_p_multi <- round(ph_multi$table[, "p"], 4)
  cat("\nPH检验P值:\n")
  cat("  单因素 RiskScore: p =", ph_p_uni, ifelse(ph_p_uni < 0.05, " (违反PH假设)", " (满足PH假设)"), "\n")
  cat("  多因素各变量:\n")
  print(ph_p_multi)
} else {
  ph_p_uni <- round(ph_uni$table[1, "p"], 4)
  cat("\nPH检验P值 (单因素): RiskScore p =", ph_p_uni, ifelse(ph_p_uni < 0.05, " (违反PH假设)", " (满足PH假设)"), "\n")
}
set.seed(123)
ph_p_boot <- replicate(100, {
  idx <- sample(nrow(merged_54), replace = TRUE)
  cox_boot <- coxph(Surv(OS.time, OS.status) ~ RiskScore, data = merged_54[idx, ])
  cox.zph(cox_boot)$table[1, "p"]
})
quantile(ph_p_boot, probs = c(0.025, 0.5, 0.975))
# ---------- 多因素模型 PH 检验，输出与 118 批次完全一致的结果 ----------
if (!is.null(cox_multi_54)) {
  ph_multi <- cox.zph(cox_multi_54)
  
  # 1. 在控制台打印（与 118 相同的格式）
  cat("\n========== 多因素 Cox 比例风险假设检验 (PH检验) ==========\n")
  print(ph_multi)
  
  # 2. 同时将输出结果保存到文本文件（便于比较）
  sink(file.path(out_dir, "Multivariate_PH_test_54.txt"))
  cat("Multivariate Cox PH Assumption Test (Validation Set, n=54)\n")
  cat("----------------------------------------------------------\n")
  print(ph_multi)
  sink()  # 恢复输出
  
  # 3. 可选：提取为数据框，方便后续查看或导出为 CSV
  ph_table <- as.data.frame(ph_multi$table)
  ph_table$variable <- rownames(ph_table)
  rownames(ph_table) <- NULL
  ph_table <- ph_table[, c("variable", "chisq", "df", "p")]
  write.csv(ph_table, file.path(out_dir, "Multivariate_PH_test_54.csv"), row.names = FALSE)
  
  cat("\nPH检验结果已保存至:\n", 
      file.path(out_dir, "Multivariate_PH_test_54.txt"), "\n",
      file.path(out_dir, "Multivariate_PH_test_54.csv"), "\n")
  
  # 4. 绘制 Schoenfeld 残差图（与 118 类似，但 118 并未画图，可选）
  pdf(file.path(out_dir, "PH_test_multi_54.pdf"), width = 8, height = 6)
  plot(ph_multi, main = "Schoenfeld Residuals - Multivariate Cox (Validation Set)")
  dev.off()
}


# ==================== timeROC (时间依赖ROC) - 54验证集 ====================
library(timeROC)

# 准备数据（已存在）
merged_roc <- na.omit(merged_54[, c("OS.time", "OS.status", "RiskScore")])

# 定义时间点（1年、3年、5年，单位：天）
time_points <- c(365, 1095, 1825)

# 计算时间依赖ROC
roc_obj <- timeROC(T = merged_roc$OS.time,
                   delta = merged_roc$OS.status,
                   marker = merged_roc$RiskScore,
                   cause = 1,
                   times = time_points,
                   iid = FALSE)   # 避免小样本报错

# 输出AUC
cat("Time-dependent AUC (1/3/5 years):\n")
print(round(roc_obj$AUC, 3))

# ---------- 绘制与118批次相同风格的ROC曲线 ----------
# 颜色与线宽（匹配之前）
colors_roc <- c("#E41A1C", "#377EB8", "#4DAF4A")
line_width <- 2.5

# 输出PDF
pdf(file.path(out_dir, "timeROC_54.pdf"), width = 6.5, height = 6)

# 建立空白画布 + 网格 + 对角线
plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
     xlab = "1 - Specificity", ylab = "Sensitivity",
     main = "", axes = FALSE)
box()
grid(col = "lightgray", lty = 1, lwd = 0.5)
axis(1, at = seq(0, 1, 0.2), cex.axis = 0.9)
axis(2, at = seq(0, 1, 0.2), las = 1, cex.axis = 0.9)
abline(a = 0, b = 1, lty = 2, col = "gray40", lwd = 1.2)

# 手动绘制三条 ROC 曲线
for (i in 1:3) {
  tpr <- roc_obj$TP[, i]   # 灵敏度
  fpr <- roc_obj$FP[, i]   # 1-特异性
  valid <- !is.na(tpr) & !is.na(fpr)
  tpr <- tpr[valid]
  fpr <- fpr[valid]
  ord <- order(fpr)
  lines(fpr[ord], tpr[ord], col = colors_roc[i], lwd = line_width)
}

# 标题
title(main = "Time-dependent ROC Curves for Risk Score (Validation Set)", line = 1.5)

# 图例（AUC值取自roc_obj$AUC）
auc_vals <- round(roc_obj$AUC, 3)
legend("bottomright", 
       legend = paste0(c("1-year", "3-year", "5-year"), " AUC = ", auc_vals),
       col = colors_roc, lty = 1, lwd = line_width,
       bty = "o", bg = "white", cex = 0.85, inset = c(0.02, 0.02))

dev.off()

cat("时间依赖ROC曲线已保存至:", file.path(out_dir, "timeROC_54.pdf"), "\n")
# ==================== DCA ====================
# 用多因素Cox模型预测3年生存概率（若多因素模型不存在，则用单因素）
if (exists("cox_multi_54") && !is.null(cox_multi_54)) {
  model_dca <- cox_multi_54
} else {
  model_dca <- cox_uni_54
}

t_pred <- 1095
# 预测每个样本的3年生存概率
sf <- survfit(model_dca, newdata = merged_54, se.fit = FALSE)
time_pts <- sf$time
surv_prob_3yr <- sapply(1:nrow(merged_54), function(i) {
  idx <- which.min(abs(time_pts - t_pred))
  sf$surv[idx, i]
})
merged_54$prob_3yr <- surv_prob_3yr
merged_54$event_3yr <- ifelse(merged_54$OS.time <= t_pred & merged_54$OS.status == 1, 1, 0)

library(rmda)
dca_model <- decision_curve(event_3yr ~ prob_3yr, data = merged_54,
                            family = binomial(link = "logit"),
                            thresholds = seq(0.01, 0.5, by = 0.01),
                            confidence.intervals = FALSE)
pdf(file.path(out_dir, "DCA_54.pdf"), width = 6, height = 5)
plot_decision_curve(dca_model, curve.names = "NETs Risk Model (3-year)", 
                    xlab = "Threshold probability", ylab = "Net benefit", col = "blue", lwd = 2)
dev.off()

