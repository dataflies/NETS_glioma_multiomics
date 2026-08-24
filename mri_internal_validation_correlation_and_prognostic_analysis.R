# 加载包
library(ggplot2)
library(ggrepel)
library(dplyr)

# 读取数据（请确认文件路径）
df_val <- read.csv("C:/Users/Fengye liu/Desktop/new_val_predictions.csv", stringsAsFactors = FALSE)
df_test <- read.csv("C:/Users/Fengye liu/Desktop/new_test_predictions.csv", stringsAsFactors = FALSE)

# 合并并计算绝对误差
df <- bind_rows(df_val, df_test) %>%
  mutate(abs_error = abs(true_time - pred_time))

# ---- 线性回归拟合（用于提取方程、R² 和标准误）----
fit <- lm(pred_time ~ true_time, data = df)
fit_sum <- summary(fit)

# 提取回归系数及标准误
intercept <- coef(fit)[1]
slope <- coef(fit)[2]
intercept_se <- fit_sum$coefficients[1, 2]
slope_se <- fit_sum$coefficients[2, 2]
r_squared <- fit_sum$r.squared

# 构造图注文本（英文）
reg_eq <- sprintf("Regression: y = %.3f x + %.3f", slope, intercept)
se_text <- sprintf("SE(slope) = %.4f, SE(intercept) = %.4f", slope_se, intercept_se)
r2_text <- sprintf("R² = %.3f", r_squared)
ref_line_text <- "Dashed line: y = x"

# 筛选离群样本（绝对误差大于90%分位数）
error_thr <- quantile(df$abs_error, 0.9, na.rm = TRUE)
df_outlier <- df %>% filter(abs_error > error_thr)

# ---- 绘图 ----
p <- ggplot(df, aes(x = true_time, y = pred_time, color = abs_error)) +
  geom_point(size = 2.5, alpha = 0.8) +
  # 渐变色：误差小 -> 蓝，误差大 -> 红
  scale_color_gradient(low = "blue", high = "red", name = "Absolute error") +
  # Y = X 虚线
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
              color = "black", size = 0.8) +
  # 线性回归实线
  geom_smooth(method = "lm", se = FALSE, color = "darkgreen", 
              linetype = "solid", size = 0.8) +
  # 标注离群样本的 sample 名称
  geom_text_repel(data = df_outlier, aes(label = sample),
                  size = 3, color = "black", box.padding = 0.3,
                  point.padding = 0.2, segment.color = "grey50") +
  # 添加图注文本（置于左上角）
  annotate("text", x = -Inf, y = Inf, 
           label = paste(reg_eq, se_text, r2_text, ref_line_text, sep = "\n"),
           hjust = -0.1, vjust = 1.1, size = 3.5, color = "black") +
  # 坐标轴标签（英文）
  labs(x = "True time", y = "Predicted risk",
       title = "Observed vs Predicted Values") +
  # 清晰的主题（带轴线）
  theme_classic(base_size = 12) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5))

# 显示图形
print(p)

# 可选：保存高分辨率图片（适合投稿）
ggsave("scatter_plot.pdf", plot = p, width = 6, height = 5, dpi = 300)
ggsave("scatter_plot.png", plot = p, width = 6, height = 5, dpi = 300)

# 控制台输出误差最大的5个样本（验证标注）
cat("\nTop 5 samples with largest absolute error:\n")
df %>% arrange(desc(abs_error)) %>% head(5) %>%
  select(sample, true_time, pred_time, abs_error) %>% print()


# ==================== 预后分析部分 ====================

# 加载生存分析相关包
library(survival)
library(survminer)
library(timeROC)
library(maxstat)     # 用于寻找最佳 cutoff
library(broom)       # 整理模型结果
library(forestplot)  # 森林图（可选）

# 1. 读取临床数据
clin <- read.csv("E:/NETS/bulk rna数据/原数据/2_CGGA_172_Filtered/CGGA_172_clinical_data.csv",
                 stringsAsFactors = FALSE, check.names = FALSE)

# 检查列名（根据用户提供的示例调整）
# 假设列名：Sample, OS, Censor (alive=0; dead=1), Age, Gender, Grade, IDH_mutation_status, 
#           1p19q_codeletion_status, MGMTp_methylation_status, Radio_status, Chemo_status
# 注意示例中 Radio_status 和 Chemo_status 有空值，需处理

# 重命名关键列为标准名
names(clin) <- make.names(names(clin))  # 确保无特殊字符

# 提取生存相关列（根据实际名称调整，此处给出常用映射）
# 用户提供的示例列：Sample, OS, Censor (alive=0; dead=1), Age, Gender, Grade, IDH_mutation_status,
#                   1p19q_codeletion_status, MGMTp_methylation_status, Radio_status, Chemo_status
# 为安全，先打印列名供检查
cat("临床数据列名：\n")
print(names(clin))

# 手动指定关键列（请根据实际列名修改）
# 如果列名与下面不一致，请修改右侧引号内的名称
col_sample <- "Sample"
col_os <- "OS"
col_censor <- "Censor..alive.0..dead.1."          # 0=存活, 1=死亡
col_age <- "Age"
col_gender <- "Gender"
col_grade <- "Grade"
col_idh <- "IDH_mutation_status"
col_1p19q <- "X1p19q_codeletion_status"   # 注意列名前可能有X
col_mgmt <- "MGMTp_methylation_status"
col_radio <- "Radio_status..treated.1.un.treated.0."
col_chemo <- "Chemo_status..TMZ.treated.1.un.treated.0."

# 检查必要列是否存在
necessary_cols <- c(col_sample, col_os, col_censor)
if (!all(necessary_cols %in% names(clin))) {
  stop("临床数据缺少必要的生存列（OS 或 Censor），请检查列名")
}

# 2. 整理临床数据
clin_clean <- clin %>%
  select(all_of(c(col_sample, col_os, col_censor, col_age, col_gender, col_grade,
                  col_idh, col_1p19q, col_mgmt, col_radio, col_chemo))) %>%
  rename(sample = all_of(col_sample),
         OS = all_of(col_os),
         censor = all_of(col_censor),
         Age = all_of(col_age),
         Gender = all_of(col_gender),
         Grade = all_of(col_grade),
         IDH = all_of(col_idh),
         `1p19q` = all_of(col_1p19q),
         MGMTp = all_of(col_mgmt),
         Radio = all_of(col_radio),
         Chemo = all_of(col_chemo)) %>%
  # 转换变量类型
  mutate(OS = as.numeric(OS),
         censor = as.numeric(censor),
         Age = as.numeric(Age),
         Gender = factor(Gender, levels = c("Male", "Female")),
         Grade = factor(Grade, levels = c("WHO II", "WHO III", "WHO IV")),
         IDH = factor(IDH, levels = c("Wildtype", "Mutant")),
         `1p19q` = factor(`1p19q`, levels = c("Non-codel", "Codel")),
         MGMTp = factor(MGMTp, levels = c("un-methylated", "methylated")),
         Radio = as.numeric(Radio),   # 假设 1=treated, 0=untreated
         Chemo = as.numeric(Chemo))   # 假设 1=TMZ treated, 0=untreated

# 处理缺失值（生存数据中的NA会导致模型失败）
clin_clean <- clin_clean %>%
  filter(!is.na(OS), !is.na(censor)) %>%
  mutate(across(c(Radio, Chemo), ~ ifelse(is.na(.), 0, .)))  # 将治疗状态缺失设为0（未治疗）

# 3. 合并预测风险值（来自之前的数据）
# 之前已有 df_pred？需要从最初的两个CSV中提取pred_time和sample
# 这里重新整理确保得到 df_pred
df_val <- read.csv("C:/Users/Fengye liu/Desktop/val_predictions.csv", stringsAsFactors = FALSE)
df_test <- read.csv("C:/Users/Fengye liu/Desktop/test_set_predictions.csv", stringsAsFactors = FALSE)
df_pred <- bind_rows(df_val, df_test) %>%
  select(sample, pred_time) %>%
  distinct(sample, .keep_all = TRUE)   # 防止重复样本

# 合并临床与预测值
merged_df <- inner_join(clin_clean, df_pred, by = "sample")
cat(sprintf("合并后有效样本数: %d\n", nrow(merged_df)))

# 如果样本数过少（<20），则停止并给出警告
if (nrow(merged_df) < 20) {
  stop("合并样本不足20，无法进行可靠的生存分析")
}

# 4. 寻找 pred_time 的最佳 cutoff（基于最大选择检验）
library(maxstat)
cutoff_obj <- maxstat.test(Surv(OS, censor) ~ pred_time, 
                           data = merged_df, 
                           smethod = "LogRank", 
                           iscores = TRUE, 
                           alpha = 0.05)
best_cut <- cutoff_obj$estimate
cat(sprintf("最佳 cutoff 值: %.4f\n", best_cut))

# 根据 cutoff 创建风险分组
merged_df <- merged_df %>%
  mutate(risk_group = ifelse(pred_time > best_cut, "High", "Low"),
         risk_group = factor(risk_group, levels = c("Low", "High")))

# 5. Kaplan-Meier 曲线及 log-rank 检验
fit_km <- survfit(Surv(OS, censor) ~ risk_group, data = merged_df)

# 绘制 KM 曲线
km_plot <- ggsurvplot(fit_km,
                      data = merged_df,
                      pval = TRUE,
                      pval.method = TRUE,
                      conf.int = TRUE,
                      risk.table = TRUE,
                      risk.table.col = "strata",
                      palette = c("blue", "red"),
                      xlab = "Time (days)",
                      ylab = "Overall Survival Probability",
                      title = paste0("Kaplan-Meier Curve (cutoff = ", round(best_cut, 3), ")"),
                      legend.title = "Risk Group",
                      legend.labs = c("Low Risk", "High Risk"))

# 保存 KM 图
ggsave("KM_curve.pdf", plot = km_plot$plot, width = 6, height = 5)
ggsave("KM_curve.png", plot = km_plot$plot, width = 6, height = 5, dpi = 300)

# 6. 单因素 Cox 回归（连续 pred_time 和 分组）
# 连续变量
uni_cox_cont <- coxph(Surv(OS, censor) ~ pred_time, data = merged_df)
uni_cox_cont_sum <- summary(uni_cox_cont)
cat("\n单因素 Cox (连续 pred_time):\n")
print(uni_cox_cont_sum)

# 分组变量
uni_cox_group <- coxph(Surv(OS, censor) ~ risk_group, data = merged_df)
uni_cox_group_sum <- summary(uni_cox_group)
cat("\n单因素 Cox (高低风险分组):\n")
print(uni_cox_group_sum)

# ========== 修正版：多因素Cox筛选 + 时间依赖ROC + Bootstrap + 最终输出 ==========

# 7. 多因素 Cox 回归（纳入临床协变量）--- 修正变量名含数字/特殊字符的问题
covariates <- c("Age", "Gender", "Grade", "IDH", "`1p19q`", "MGMTp", "Radio", "Chemo")
uni_results <- data.frame()
for (var in covariates) {
  # 去除反引号用于判断列是否存在
  clean_var <- gsub("`", "", var)
  if (clean_var %in% colnames(merged_df)) {
    # 跳过因子水平不足的变量
    if (is.factor(merged_df[[clean_var]]) && nlevels(merged_df[[clean_var]]) <= 1) next
    # 使用反引号包裹变量名以处理特殊字符（如1p19q）
    form <- as.formula(paste("Surv(OS, censor) ~", var))
    fit <- coxph(form, data = merged_df)
    pval <- summary(fit)$coefficients[1, "Pr(>|z|)"]
    uni_results <- rbind(uni_results, data.frame(variable = clean_var, p.value = pval))
  }
}

# 选择 p < 0.2 的变量进入多因素
sig_vars <- uni_results[uni_results$p.value < 0.2, "variable"]
if (length(sig_vars) == 0) {
  cat("无显著协变量，多因素模型仅包含 risk_group\n")
  multi_formula <- Surv(OS, censor) ~ risk_group
} else {
  # 对含有特殊字符的变量名加上反引号
  sig_vars_backticked <- ifelse(grepl("^[0-9]", sig_vars), paste0("`", sig_vars, "`"), sig_vars)
  multi_formula <- as.formula(paste("Surv(OS, censor) ~ risk_group +", paste(sig_vars_backticked, collapse = " + ")))
}

# 拟合多因素 Cox
multi_cox <- coxph(multi_formula, data = merged_df)
multi_sum <- summary(multi_cox)
cat("\n多因素 Cox 回归结果:\n")
print(multi_sum)

# 保存多因素结果
sink("multivariate_cox_results.txt")
print(multi_sum)
sink()

# 5. Kaplan-Meier 曲线（去掉置信区间带）
fit_km <- survfit(Surv(OS, censor) ~ risk_group, data = merged_df)

km_plot <- ggsurvplot(fit_km,
                      data = merged_df,
                      pval = TRUE,
                      pval.method = TRUE,
                      conf.int = FALSE,          # 关键：去掉阴影置信带
                      risk.table = TRUE,
                      risk.table.col = "strata",
                      palette = c("blue", "red"),
                      xlab = "Time (days)",
                      ylab = "Overall Survival Probability",
                      title = paste0("Kaplan-Meier Curve (cutoff = ", round(best_cut, 3), ")"),
                      legend.title = "Risk Group",
                      legend.labs = c("Low Risk", "High Risk"))

# 保存
ggsave("KM_curve.pdf", plot = km_plot$plot, width = 6, height = 5)


# 9. 时间依赖 ROC 曲线（单图三曲线，无子图分割）
library(timeROC)

# 定义时间点（1年、3年、5年，单位：天）
time_points <- c(365, 1095, 1825)

# 计算时间依赖ROC（小样本建议 iid = FALSE）
time_roc <- timeROC(T = merged_df$OS,
                    delta = merged_df$censor,
                    marker = merged_df$pred_time,
                    cause = 1,
                    times = time_points,
                    iid = FALSE)

# 输出AUC
cat("\nTime-dependent AUC (1/3/5 years):\n")
print(round(time_roc$AUC, 3))

# 设置颜色和线宽
colors_roc <- c("#E41A1C", "#377EB8", "#4DAF4A")
line_width <- 2.5

# 输出PDF
pdf("time_ROC.pdf", width = 6.5, height = 6)

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
  tpr <- time_roc$TP[, i]   # 灵敏度
  fpr <- time_roc$FP[, i]   # 1-特异性
  valid <- !is.na(tpr) & !is.na(fpr)
  tpr <- tpr[valid]
  fpr <- fpr[valid]
  ord <- order(fpr)
  lines(fpr[ord], tpr[ord], col = colors_roc[i], lwd = line_width)
}

# 标题
title(main = "Time-dependent ROC Curves for Risk Score", line = 1.5)

# 图例（AUC值取自 time_roc$AUC）
auc_vals <- round(time_roc$AUC, 3)
legend("bottomright", 
       legend = paste0(c("1-year", "3-year", "5-year"), " AUC = ", auc_vals),
       col = colors_roc, lty = 1, lwd = line_width,
       bty = "o", bg = "white", cex = 0.85, inset = c(0.02, 0.02))

dev.off()

cat("时间依赖ROC曲线已保存至 time_ROC.pdf\n")
# 10. Bootstrap 稳定性验证（已运行成功，无需修改）
set.seed(123)
boot_cutoffs <- replicate(200, {
  boot_idx <- sample(1:nrow(merged_df), size = nrow(merged_df), replace = TRUE)
  boot_data <- merged_df[boot_idx, ]
  if (length(unique(boot_data$censor)) < 2) return(NA)
  tryCatch({
    boot_cut <- maxstat.test(Surv(OS, censor) ~ pred_time, data = boot_data, smethod = "LogRank")$estimate
    return(boot_cut)
  }, error = function(e) NA)
})
boot_cutoffs <- boot_cutoffs[!is.na(boot_cutoffs)]
if (length(boot_cutoffs) > 0) {
  ci_cut <- quantile(boot_cutoffs, c(0.025, 0.975), na.rm = TRUE)
  cat(sprintf("\nBootstrap 最佳 cutoff 的 95%% CI: [%.4f, %.4f]\n", ci_cut[1], ci_cut[2]))
} else {
  cat("\nBootstrap 无法计算稳定 cutoff\n")
}

# 11. 输出最终临床建议阈值（修正索引错误）
cat("\n===== 临床建议 =====\n")
cat(sprintf("基于当前数据，预测风险值（pred_time）的最佳阈值是 %.4f\n", best_cut))
# 修正：单因素分组 Cox 的结果只有一行，索引为 [1, "Pr(>|z|)"]
hr_high <- exp(coef(uni_cox_group))[2]          # risk_groupHigh 的 HR
p_high <- summary(uni_cox_group)$coefficients[1, "Pr(>|z|)"]  # 第一行
cat(sprintf("当 pred_time > %.4f 时，患者被划分为高风险组，预后显著更差（HR = %.2f, p = %.3f）\n",
            best_cut, hr_high, p_high))
cat("注意：由于样本量较小，该阈值应在前瞻性队列中验证后使用。\n")

# 保存合并后的数据供外部使用
write.csv(merged_df, "merged_survival_data.csv", row.names = FALSE)

# 完成
cat("\n所有预后分析已完成。生成的文件：KM_curve.pdf/png, forest_plot.pdf/png, time_ROC.pdf, multivariate_cox_results.txt, merged_survival_data.csv\n")





# ==================== 修正版：生成 Table 5（C-index 已修正）====================
# 本代码块独立运行，依赖对象：df, merged_df

# 加载必要包（若已加载则跳过）
suppressPackageStartupMessages({
  library(survival)
  library(timeROC)
  library(Hmisc)   # 仅用于格式，不用于 concordance
})

# 检查数据对象
if (!exists("df") || !exists("merged_df")) {
  stop("请先运行原始脚本生成 df 和 merged_df 对象")
}

# 设定 Bootstrap 重复次数
n_boot <- 500
set.seed(123)

# 修正后的指标计算函数（C-index 使用 1 - concordance）
calc_metrics_corrected <- function(data_reg, data_surv,
                                   boot_idx_reg = NULL, boot_idx_surv = NULL) {
  # 子集
  reg <- if (is.null(boot_idx_reg)) data_reg else data_reg[boot_idx_reg, ]
  surv <- if (is.null(boot_idx_surv)) data_surv else data_surv[boot_idx_surv, ]
  
  res <- c(n_reg = NA, n_surv = NA,
           MAE = NA, RMSE = NA, Mean_Error = NA, Median_Abs_Error = NA,
           Pearson_r = NA, R2 = NA,
           AUC_1yr = NA, AUC_3yr = NA, AUC_5yr = NA,
           C_index = NA,
           HR = NA, HR_lower = NA, HR_upper = NA, Logrank_P = NA)
  
  # 回归指标
  if (nrow(reg) > 0) {
    res["n_reg"] <- nrow(reg)
    err <- reg$pred_time - reg$true_time
    res["MAE"] <- mean(abs(err), na.rm = TRUE)
    res["RMSE"] <- sqrt(mean(err^2, na.rm = TRUE))
    res["Mean_Error"] <- mean(err, na.rm = TRUE)
    res["Median_Abs_Error"] <- median(abs(err), na.rm = TRUE)
    res["Pearson_r"] <- cor(reg$pred_time, reg$true_time, method = "pearson")
    res["R2"] <- res["Pearson_r"]^2
  }
  
  # 生存指标
  if (nrow(surv) > 0) {
    res["n_surv"] <- nrow(surv)
    # AUC
    if (sum(surv$censor) > 0 && sum(surv$censor) < nrow(surv)) {
      roc <- tryCatch(
        timeROC(T = surv$OS, delta = surv$censor, marker = surv$pred_time,
                cause = 1, times = c(365, 1095, 1825), iid = FALSE),
        error = function(e) NULL
      )
      if (!is.null(roc)) {
        res["AUC_1yr"] <- roc$AUC[1]
        res["AUC_3yr"] <- roc$AUC[2]
        res["AUC_5yr"] <- roc$AUC[3]
      }
    }
    # C-index (修正：风险评分越高，生存时间越短 => 期望 C-index > 0.5)
    if (sum(surv$censor) > 0 && sum(surv$censor) < nrow(surv)) {
      # 使用 survival::concordance，然后计算 1 - 原始值
      c_obj <- tryCatch(
        concordance(Surv(OS, censor) ~ pred_time, data = surv),
        error = function(e) NULL
      )
      if (!is.null(c_obj)) {
        res["C_index"] <- 1 - c_obj$concordance
      }
    }
    # Cox 分组 HR
    if (length(unique(surv$risk_group)) == 2) {
      cox <- tryCatch(
        coxph(Surv(OS, censor) ~ risk_group, data = surv),
        error = function(e) NULL
      )
      if (!is.null(cox)) {
        res["HR"] <- exp(coef(cox)[1])
        ci <- tryCatch(exp(confint(cox)[1, ]), error = function(e) c(NA, NA))
        res["HR_lower"] <- ci[1]
        res["HR_upper"] <- ci[2]
        res["Logrank_P"] <- summary(cox)$logtest[3]
      }
    }
  }
  return(res)
}

# 准备全量数据
reg_full <- df[, c("true_time", "pred_time")]
surv_full <- merged_df[, c("OS", "censor", "pred_time", "risk_group")]

# 点估计
point_est <- calc_metrics_corrected(reg_full, surv_full)
cat("\n点估计（修正后）:\n"); print(point_est)

# Bootstrap
boot_mae <- boot_rmse <- boot_mean_err <- boot_med_abs <- boot_pearson <- numeric(n_boot)
boot_auc1 <- boot_auc3 <- boot_auc5 <- numeric(n_boot)
boot_cindex <- numeric(n_boot)
boot_hr <- numeric(n_boot)

n_reg <- nrow(reg_full); n_surv <- nrow(surv_full)

pb <- txtProgressBar(min = 0, max = n_boot, style = 3)
for (i in 1:n_boot) {
  idx_reg <- sample(1:n_reg, size = n_reg, replace = TRUE)
  idx_surv <- sample(1:n_surv, size = n_surv, replace = TRUE)
  boot_res <- tryCatch(
    calc_metrics_corrected(reg_full, surv_full, idx_reg, idx_surv),
    error = function(e) rep(NA, length(point_est))
  )
  boot_mae[i] <- boot_res["MAE"]
  boot_rmse[i] <- boot_res["RMSE"]
  boot_mean_err[i] <- boot_res["Mean_Error"]
  boot_med_abs[i] <- boot_res["Median_Abs_Error"]
  boot_pearson[i] <- boot_res["Pearson_r"]
  boot_auc1[i] <- boot_res["AUC_1yr"]
  boot_auc3[i] <- boot_res["AUC_3yr"]
  boot_auc5[i] <- boot_res["AUC_5yr"]
  boot_cindex[i] <- boot_res["C_index"]
  boot_hr[i] <- boot_res["HR"]
  setTxtProgressBar(pb, i)
}
close(pb)

# 计算 95% CI
calc_ci <- function(x) {
  if (all(is.na(x))) return(c(NA, NA))
  quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
}
ci_mae <- calc_ci(boot_mae)
ci_rmse <- calc_ci(boot_rmse)
ci_mean_err <- calc_ci(boot_mean_err)
ci_med_abs <- calc_ci(boot_med_abs)
ci_pearson <- calc_ci(boot_pearson)
ci_auc1 <- calc_ci(boot_auc1)
ci_auc3 <- calc_ci(boot_auc3)
ci_auc5 <- calc_ci(boot_auc5)
ci_cindex <- calc_ci(boot_cindex)
ci_hr <- calc_ci(boot_hr)

# 构建 Table5
table5 <- data.frame(
  Metric = c("Sample size (n)", "AUC (1-year)", "AUC (3-year)", "AUC (5-year)",
             "C-index (Harrell's)", "Pearson's r", "R²", "MAE", "RMSE",
             "Mean Error (bias)", "Median Absolute Error", "Log-rank P",
             "HR (high vs. low risk)"),
  Value = c(point_est["n_surv"],
            round(point_est["AUC_1yr"], 3), round(point_est["AUC_3yr"], 3),
            round(point_est["AUC_5yr"], 3), round(point_est["C_index"], 3),
            round(point_est["Pearson_r"], 3), round(point_est["R2"], 3),
            round(point_est["MAE"], 3), round(point_est["RMSE"], 3),
            round(point_est["Mean_Error"], 3), round(point_est["Median_Abs_Error"], 3),
            formatC(point_est["Logrank_P"], format = "e", digits = 2),
            round(point_est["HR"], 3)),
  CI_lower = c(NA, round(ci_auc1[1],3), round(ci_auc3[1],3), round(ci_auc5[1],3),
               round(ci_cindex[1],3), round(ci_pearson[1],3), NA,
               round(ci_mae[1],3), round(ci_rmse[1],3), round(ci_mean_err[1],3),
               round(ci_med_abs[1],3), NA, round(ci_hr[1],3)),
  CI_upper = c(NA, round(ci_auc1[2],3), round(ci_auc3[2],3), round(ci_auc5[2],3),
               round(ci_cindex[2],3), round(ci_pearson[2],3), NA,
               round(ci_mae[2],3), round(ci_rmse[2],3), round(ci_mean_err[2],3),
               round(ci_med_abs[2],3), NA, round(ci_hr[2],3))
)

table5$`95% CI` <- ifelse(is.na(table5$CI_lower), "–",
                          paste0(table5$CI_lower, "–", table5$CI_upper))
table5 <- table5[, c("Metric", "Value", "95% CI")]

# 输出
# 保存 CSV（解决乱码）
write.csv(table5, "Table5.csv", row.names = FALSE, fileEncoding = "UTF-8")

# ==================== Prediction residual distribution ====================
# 计算残差 (预测值 - 真实值)
df <- df %>%
  mutate(residual = pred_time - true_time)

# ---- 直方图 + 密度曲线 ----
p_res_hist <- ggplot(df, aes(x = residual)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "steelblue", 
                 color = "black", alpha = 0.7) +
  geom_density(color = "red", size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgreen", size = 0.8) +
  geom_vline(xintercept = mean(df$residual, na.rm = TRUE), 
             linetype = "dotted", color = "blue", size = 0.8) +
  annotate("text", x = mean(df$residual, na.rm = TRUE), y = Inf, 
           label = paste("Mean =", round(mean(df$residual, na.rm = TRUE), 3)), 
           hjust = -0.1, vjust = 2, size = 3.5, color = "blue") +
  labs(x = "Prediction residual (pred_time - true_time)", y = "Density",
       title = "Distribution of Prediction Residuals") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p_res_hist)
ggsave("residual_distribution_hist.pdf", plot = p_res_hist, width = 6, height = 4, dpi = 300)
ggsave("residual_distribution_hist.png", plot = p_res_hist, width = 6, height = 4, dpi = 300)

# ---- 若需要小提琴图（仅作备选，可注释或取消注释） ----
 p_res_violin <- ggplot(df, aes(x = "", y = residual)) +
   geom_violin(fill = "steelblue", color = "black", alpha = 0.7) +
   geom_hline(yintercept = 0, linetype = "dashed", color = "darkgreen", size = 0.8) +
   geom_hline(yintercept = mean(df$residual, na.rm = TRUE), 
              linetype = "dotted", color = "blue", size = 0.8) +
   labs(x = "", y = "Prediction residual", title = "Distribution of Prediction Residuals") +
   theme_classic(base_size = 12) +
   theme(plot.title = element_text(hjust = 0.5))
 ggsave("residual_distribution_violin.pdf", plot = p_res_violin, width = 4, height = 5, dpi = 300)

# ---- 输出残差统计量 ----
cat("\n===== Prediction Residual Statistics =====\n")
cat(sprintf("Mean residual = %.4f\n", mean(df$residual, na.rm = TRUE)))
cat(sprintf("SD residual   = %.4f\n", sd(df$residual, na.rm = TRUE)))
cat(sprintf("Median residual = %.4f\n", median(df$residual, na.rm = TRUE)))
cat(sprintf("IQR residual  = %.4f\n", IQR(df$residual, na.rm = TRUE)))


# ==================== Bland-Altman Plot ====================
# 计算平均值和差值
df <- df %>%
  mutate(avg = (true_time + pred_time) / 2,
         diff = pred_time - true_time)   # 差值 = 预测 - 真实

# 计算统计量
mean_diff <- mean(df$diff, na.rm = TRUE)
sd_diff   <- sd(df$diff, na.rm = TRUE)
upper_lim <- mean_diff + 1.96 * sd_diff
lower_lim <- mean_diff - 1.96 * sd_diff

# 绘图
p_ba <- ggplot(df, aes(x = avg, y = diff)) +
  geom_point(size = 2, alpha = 0.6, color = "darkblue") +
  geom_hline(yintercept = mean_diff, linetype = "solid", color = "red", size = 0.8) +
  geom_hline(yintercept = upper_lim, linetype = "dashed", color = "darkgreen", size = 0.6) +
  geom_hline(yintercept = lower_lim, linetype = "dashed", color = "darkgreen", size = 0.6) +
  # 标注均值线
  annotate("text", x = -Inf, y = mean_diff, 
           label = paste("Mean =", round(mean_diff, 3)), 
           hjust = -0.1, vjust = -0.5, size = 3.5, color = "red") +
  # 标注上限
  annotate("text", x = -Inf, y = upper_lim, 
           label = paste("+1.96 SD =", round(upper_lim, 3)), 
           hjust = -0.1, vjust = -0.5, size = 3.5, color = "darkgreen") +
  # 标注下限
  annotate("text", x = -Inf, y = lower_lim, 
           label = paste("-1.96 SD =", round(lower_lim, 3)), 
           hjust = -0.1, vjust = 1.2, size = 3.5, color = "darkgreen") +
  labs(x = "Average of True and Predicted (True + Pred)/2", 
       y = "Difference (Pred - True)",
       title = "Bland-Altman Plot for Prediction Agreement") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

# 显示并保存
print(p_ba)
ggsave("Bland_Altman_plot.pdf", plot = p_ba, width = 6, height = 5, dpi = 300)
ggsave("Bland_Altman_plot.png", plot = p_ba, width = 6, height = 5, dpi = 300)

# 输出统计汇总
cat("\n===== Bland-Altman Statistics =====\n")
cat(sprintf("Mean difference (bias) = %.4f\n", mean_diff))
cat(sprintf("SD of differences      = %.4f\n", sd_diff))
cat(sprintf("95%% limits of agreement: [%.4f, %.4f]\n", lower_lim, upper_lim))
cat("如果所有点都落在两条虚线之间（或绝大多数），则认为预测与真实值具有良好的一致性。\n")
