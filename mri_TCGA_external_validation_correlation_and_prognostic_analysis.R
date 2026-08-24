# ============================================================================
# 修改：从两个文件分别读取预测值和真实值，通过样本ID合并
# ============================================================================
# 加载包
library(ggplot2)
library(ggrepel)
library(dplyr)
library(survival)
library(survminer)
library(timeROC)
library(maxstat)
library(broom)
library(forestplot)   # 可选
library(dplyr)

# ----------------------------------------------------------------------------
# 1. 读取并合并预测值和真实值（关闭列名自动转换避免编码报错）
# ----------------------------------------------------------------------------

# 预测结果文件
pred_df <- read.csv("C:/Users/Fengye liu/Desktop/predicted_risk_scores_v2_crop.csv",
                    stringsAsFactors = FALSE,
                    check.names = FALSE)

# 真实值文件
true_df <- read.csv("E:/NETS/bulk&影像 TCGA外部验证/原数据/1_TCGA_Filtered_Data/1_TCGA_TCGA_risk-图层.csv",
                    stringsAsFactors = FALSE,
                    check.names = FALSE)

# ----- 调试：查看列名（已确认，可保留或注释掉）-----
 cat("预测文件列名：\n"); print(names(pred_df))
 cat("真实值文件列名：\n"); print(names(true_df))

# ----- 根据实际列名重命名 -----
pred_df <- pred_df %>%
  rename(sample = pid, pred_time = predicted_risk) %>%
  select(sample, pred_time)

true_df <- true_df %>%
  rename(sample = sample_id, true_time = risk_score) %>%
  select(sample, true_time)

# 合并（内连接，只保留两个文件共有的样本）
df <- inner_join(true_df, pred_df, by = "sample") %>%
  distinct(sample, .keep_all = TRUE)   # 去重

cat("合并后样本数：", nrow(df), "\n")
head(df)

# 检查缺失值（若有则删除）
if (any(is.na(df$true_time)) || any(is.na(df$pred_time))) {
  warning("存在缺失值，将删除包含NA的行")
  df <- df %>% filter(!is.na(true_time), !is.na(pred_time))
}
# ----------------------------------------------------------------------------
# 2. 散点图（Observed vs Predicted）及回归诊断
# ----------------------------------------------------------------------------
df <- df %>%
  mutate(abs_error = abs(true_time - pred_time))

fit <- lm(pred_time ~ true_time, data = df)
fit_sum <- summary(fit)

intercept <- coef(fit)[1]
slope <- coef(fit)[2]
intercept_se <- fit_sum$coefficients[1, 2]
slope_se <- fit_sum$coefficients[2, 2]
r_squared <- fit_sum$r.squared

reg_eq <- sprintf("Regression: y = %.3f x + %.3f", slope, intercept)
se_text <- sprintf("SE(slope) = %.4f, SE(intercept) = %.4f", slope_se, intercept_se)
r2_text <- sprintf("R² = %.3f", r_squared)
ref_line_text <- "Dashed line: y = x"

error_thr <- quantile(df$abs_error, 0.9, na.rm = TRUE)
df_outlier <- df %>% filter(abs_error > error_thr)

p <- ggplot(df, aes(x = true_time, y = pred_time, color = abs_error)) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_color_gradient(low = "blue", high = "red", name = "Absolute error") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "darkgreen", linetype = "solid", size = 0.8) +
  geom_text_repel(data = df_outlier, aes(label = sample),
                  size = 3, color = "black", box.padding = 0.3,
                  point.padding = 0.2, segment.color = "grey50") +
  annotate("text", x = -Inf, y = Inf, 
           label = paste(reg_eq, se_text, r2_text, ref_line_text, sep = "\n"),
           hjust = -0.1, vjust = 1.1, size = 3.5, color = "black") +
  labs(x = "True time", y = "Predicted risk",
       title = "Observed vs Predicted Values") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5))

print(p)
ggsave("scatter_plot.pdf", plot = p, width = 6, height = 5, dpi = 300)
ggsave("scatter_plot.png", plot = p, width = 6, height = 5, dpi = 300)

cat("\nTop 5 samples with largest absolute error:\n")
df %>% arrange(desc(abs_error)) %>% head(5) %>%
  select(sample, true_time, pred_time, abs_error) %>% print()

# ============================================================================
# 3. 预后分析部分（生存分析）
# ============================================================================

# ============================================================================
# 3. 预后分析部分（生存分析）
# ============================================================================

# 3.1 读取临床数据（路径已更新，check.names = FALSE 防止列名编码问题）
clin <- read.csv("E:/NETS/bulk&影像 TCGA外部验证/原数据/1_TCGA_Filtered_Data/1_TCGA_TCGA_clinical_with_risk.csv",
                 stringsAsFactors = FALSE, check.names = FALSE)

cat("临床数据列名：\n")
print(names(clin))

# ----- 根据实际列名调整映射（关键修改）-----
col_sample <- "Patient_ID"                     # 样本ID列
col_os <- "OS"                                 # 生存时间
col_censor <- "Censor"                         # 生存状态（0=存活，1=死亡）
col_age <- "Age"
col_gender <- "Gender"
col_grade <- "Grade"
col_idh <- "IDH_mutation_status"
col_1p19q <- "1p19q_codeletion_status"         # 注意列名以数字开头，需用反引号包裹
# 以下列在数据中不存在，设为 NULL 或直接删除
col_mgmt <- NULL
col_radio <- NULL
col_chemo <- NULL

# 检查必要列是否存在
necessary_cols <- c(col_sample, col_os, col_censor)
if (!all(necessary_cols %in% names(clin))) {
  stop("临床数据缺少必要的生存列（OS 或 Censor），请检查列名")
}

# 3.2 整理临床数据（只选择存在的列，并处理缺失的协变量）
# 构建要选择的列名向量（只包括存在的列）
select_cols <- c(col_sample, col_os, col_censor, col_age, col_gender, col_grade,
                 col_idh, col_1p19q)
# 过滤掉不存在的列（如 col_mgmt 等为 NULL 则自动忽略）
select_cols <- select_cols[!sapply(select_cols, is.null)]
select_cols <- select_cols[select_cols %in% names(clin)]

clin_clean <- clin %>%
  select(all_of(select_cols)) %>%
  rename(sample = all_of(col_sample),
         OS = all_of(col_os),
         censor = all_of(col_censor),
         Age = all_of(col_age),
         Gender = all_of(col_gender),
         Grade = all_of(col_grade),
         IDH = all_of(col_idh),
         `1p19q` = all_of(col_1p19q)) %>%
  mutate(OS = as.numeric(OS),
         censor = as.numeric(censor),
         Age = as.numeric(Age),
         Gender = factor(Gender, levels = c("Male", "Female")),
         Grade = factor(Grade, levels = c("WHO II", "WHO III", "WHO IV")),
         IDH = factor(IDH, levels = c("Wildtype", "Mutant")),
         `1p19q` = factor(`1p19q`, levels = c("Non-codel", "Codel")))
# 注意：如果没有 MGMTp、Radio、Chemo，则不再处理它们

# 删除生存数据中的缺失值
clin_clean <- clin_clean %>%
  filter(!is.na(OS), !is.na(censor))

# 3.3 合并预测风险值（直接使用已合并的 df）
df_pred <- df %>% select(sample, pred_time) %>% distinct(sample, .keep_all = TRUE)

merged_df <- inner_join(clin_clean, df_pred, by = "sample")
cat(sprintf("合并后有效样本数: %d\n", nrow(merged_df)))

if (nrow(merged_df) < 20) {
  stop("合并样本不足20，无法进行可靠的生存分析")
}

# 3.3 寻找最佳 cutoff
cutoff_obj <- maxstat.test(Surv(OS, censor) ~ pred_time, 
                           data = merged_df, 
                           smethod = "LogRank", 
                           iscores = TRUE, 
                           alpha = 0.05)
best_cut <- cutoff_obj$estimate
cat(sprintf("最佳 cutoff 值: %.4f\n", best_cut))

merged_df <- merged_df %>%
  mutate(risk_group = ifelse(pred_time > best_cut, "High", "Low"),
         risk_group = factor(risk_group, levels = c("Low", "High")))

# 3.4 Kaplan-Meier 曲线
fit_km <- survfit(Surv(OS, censor) ~ risk_group, data = merged_df)

km_plot <- ggsurvplot(fit_km,
                      data = merged_df,
                      pval = TRUE,
                      pval.method = TRUE,
                      conf.int = FALSE,
                      risk.table = TRUE,
                      risk.table.col = "strata",
                      palette = c("blue", "red"),
                      xlab = "Time (days)",
                      ylab = "Overall Survival Probability",
                      title = paste0("Kaplan-Meier Curve (cutoff = ", round(best_cut, 3), ")"),
                      legend.title = "Risk Group",
                      legend.labs = c("Low Risk", "High Risk"))

ggsave("KM_curve.pdf", plot = km_plot$plot, width = 6, height = 5)
ggsave("KM_curve.png", plot = km_plot$plot, width = 6, height = 5, dpi = 300)

# 3.5 单因素 Cox
uni_cox_cont <- coxph(Surv(OS, censor) ~ pred_time, data = merged_df)
cat("\n单因素 Cox (连续 pred_time):\n")
print(summary(uni_cox_cont))

uni_cox_group <- coxph(Surv(OS, censor) ~ risk_group, data = merged_df)
cat("\n单因素 Cox (高低风险分组):\n")
print(summary(uni_cox_group))

# 3.6 多因素 Cox（筛选 p<0.2 的协变量，仅使用实际存在的列）
# 定义候选协变量（与 merged_df 中的列名一致）
candidate_vars <- c("Age", "Gender", "Grade", "IDH", "1p19q")

uni_results <- data.frame()

for (var in candidate_vars) {
  # 检查列是否存在
  if (!(var %in% colnames(merged_df))) next
  
  # 检查是否有足够的非缺失观测
  if (is.factor(merged_df[[var]])) {
    # 因子型：至少2个水平，且每个水平至少2个观测
    if (nlevels(merged_df[[var]]) <= 1) next
    if (any(table(merged_df[[var]]) < 2)) next
  } else {
    # 数值型：非缺失观测数 >= 5，且至少有2个不同值
    if (sum(!is.na(merged_df[[var]])) < 5) next
    if (length(unique(merged_df[[var]][!is.na(merged_df[[var]])])) < 2) next
  }
  
  # 构造公式：对以数字开头的列名（如1p19q）添加反引号
  var_formula <- if (grepl("^[0-9]", var)) paste0("`", var, "`") else var
  form <- as.formula(paste("Surv(OS, censor) ~", var_formula))
  
  # 尝试拟合，捕获错误
  fit <- tryCatch(
    coxph(form, data = merged_df),
    error = function(e) NULL
  )
  if (is.null(fit)) next
  
  pval <- summary(fit)$coefficients[1, "Pr(>|z|)"]
  uni_results <- rbind(uni_results, data.frame(variable = var, p.value = pval))
}

# 选择 p < 0.2 的变量进入多因素
sig_vars <- uni_results[uni_results$p.value < 0.2, "variable"]
if (length(sig_vars) == 0) {
  cat("无显著协变量，多因素模型仅包含 risk_group\n")
  multi_formula <- Surv(OS, censor) ~ risk_group
} else {
  # 对特殊列名添加反引号
  sig_vars_backticked <- ifelse(grepl("^[0-9]", sig_vars), paste0("`", sig_vars, "`"), sig_vars)
  multi_formula <- as.formula(paste("Surv(OS, censor) ~ risk_group +", 
                                    paste(sig_vars_backticked, collapse = " + ")))
}

# 拟合多因素 Cox
multi_cox <- coxph(multi_formula, data = merged_df)
cat("\n多因素 Cox 回归结果:\n")
print(summary(multi_cox))
sink("multivariate_cox_results.txt")
print(summary(multi_cox))
sink()

# 3.7 时间依赖 ROC
time_points <- c(365, 1095, 1825)
time_roc <- timeROC(T = merged_df$OS,
                    delta = merged_df$censor,
                    marker = merged_df$pred_time,
                    cause = 1,
                    times = time_points,
                    iid = FALSE)
cat("\nTime-dependent AUC (1/3/5 years):\n")
print(round(time_roc$AUC, 3))

colors_roc <- c("#E41A1C", "#377EB8", "#4DAF4A")
pdf("time_ROC.pdf", width = 6.5, height = 6)
plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
     xlab = "1 - Specificity", ylab = "Sensitivity",
     main = "", axes = FALSE)
box()
grid(col = "lightgray", lty = 1, lwd = 0.5)
axis(1, at = seq(0, 1, 0.2), cex.axis = 0.9)
axis(2, at = seq(0, 1, 0.2), las = 1, cex.axis = 0.9)
abline(a = 0, b = 1, lty = 2, col = "gray40", lwd = 1.2)
for (i in 1:3) {
  tpr <- time_roc$TP[, i]
  fpr <- time_roc$FP[, i]
  valid <- !is.na(tpr) & !is.na(fpr)
  tpr <- tpr[valid]
  fpr <- fpr[valid]
  ord <- order(fpr)
  lines(fpr[ord], tpr[ord], col = colors_roc[i], lwd = 2.5)
}
title(main = "Time-dependent ROC Curves for Risk Score", line = 1.5)
auc_vals <- round(time_roc$AUC, 3)
legend("bottomright", 
       legend = paste0(c("1-year", "3-year", "5-year"), " AUC = ", auc_vals),
       col = colors_roc, lty = 1, lwd = 2.5,
       bty = "o", bg = "white", cex = 0.85, inset = c(0.02, 0.02))
dev.off()

# 3.8 Bootstrap 验证 cutoff
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

# 3.9 输出临床建议
cat("\n===== 临床建议 =====\n")
cat(sprintf("基于当前数据，预测风险值（pred_time）的最佳阈值是 %.4f\n", best_cut))
hr_high <- exp(coef(uni_cox_group))[2]
p_high <- summary(uni_cox_group)$coefficients[1, "Pr(>|z|)"]
cat(sprintf("当 pred_time > %.4f 时，患者被划分为高风险组，预后显著更差（HR = %.2f, p = %.3f）\n",
            best_cut, hr_high, p_high))
cat("注意：由于样本量较小，该阈值应在前瞻性队列中验证后使用。\n")
write.csv(merged_df, "merged_survival_data.csv", row.names = FALSE)

# ============================================================================
# 4. 生成 Table 5（性能汇总表）
# ============================================================================
suppressPackageStartupMessages({
  library(survival)
  library(timeROC)
  library(Hmisc)
})

n_boot <- 500
set.seed(123)

calc_metrics_corrected <- function(data_reg, data_surv,
                                   boot_idx_reg = NULL, boot_idx_surv = NULL) {
  reg <- if (is.null(boot_idx_reg)) data_reg else data_reg[boot_idx_reg, ]
  surv <- if (is.null(boot_idx_surv)) data_surv else data_surv[boot_idx_surv, ]
  
  res <- c(n_reg = NA, n_surv = NA,
           MAE = NA, RMSE = NA, Mean_Error = NA, Median_Abs_Error = NA,
           Pearson_r = NA, R2 = NA,
           AUC_1yr = NA, AUC_3yr = NA, AUC_5yr = NA,
           C_index = NA,
           HR = NA, HR_lower = NA, HR_upper = NA, Logrank_P = NA)
  
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
  
  if (nrow(surv) > 0) {
    res["n_surv"] <- nrow(surv)
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
    if (sum(surv$censor) > 0 && sum(surv$censor) < nrow(surv)) {
      c_obj <- tryCatch(
        concordance(Surv(OS, censor) ~ pred_time, data = surv),
        error = function(e) NULL
      )
      if (!is.null(c_obj)) {
        res["C_index"] <- 1 - c_obj$concordance
      }
    }
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

reg_full <- df[, c("true_time", "pred_time")]
surv_full <- merged_df[, c("OS", "censor", "pred_time", "risk_group")]

point_est <- calc_metrics_corrected(reg_full, surv_full)
cat("\n点估计（修正后）:\n"); print(point_est)

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
write.csv(table5, "Table5.csv", row.names = FALSE, fileEncoding = "UTF-8")

# ============================================================================
# 5. 残差分布图
# ============================================================================
df <- df %>%
  mutate(residual = pred_time - true_time)

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

cat("\n===== Prediction Residual Statistics =====\n")
cat(sprintf("Mean residual = %.4f\n", mean(df$residual, na.rm = TRUE)))
cat(sprintf("SD residual   = %.4f\n", sd(df$residual, na.rm = TRUE)))
cat(sprintf("Median residual = %.4f\n", median(df$residual, na.rm = TRUE)))
cat(sprintf("IQR residual  = %.4f\n", IQR(df$residual, na.rm = TRUE)))

# ============================================================================
# 6. Bland-Altman 图
# ============================================================================
df <- df %>%
  mutate(avg = (true_time + pred_time) / 2,
         diff = pred_time - true_time)

mean_diff <- mean(df$diff, na.rm = TRUE)
sd_diff   <- sd(df$diff, na.rm = TRUE)
upper_lim <- mean_diff + 1.96 * sd_diff
lower_lim <- mean_diff - 1.96 * sd_diff

p_ba <- ggplot(df, aes(x = avg, y = diff)) +
  geom_point(size = 2, alpha = 0.6, color = "darkblue") +
  geom_hline(yintercept = mean_diff, linetype = "solid", color = "red", size = 0.8) +
  geom_hline(yintercept = upper_lim, linetype = "dashed", color = "darkgreen", size = 0.6) +
  geom_hline(yintercept = lower_lim, linetype = "dashed", color = "darkgreen", size = 0.6) +
  annotate("text", x = -Inf, y = mean_diff, 
           label = paste("Mean =", round(mean_diff, 3)), 
           hjust = -0.1, vjust = -0.5, size = 3.5, color = "red") +
  annotate("text", x = -Inf, y = upper_lim, 
           label = paste("+1.96 SD =", round(upper_lim, 3)), 
           hjust = -0.1, vjust = -0.5, size = 3.5, color = "darkgreen") +
  annotate("text", x = -Inf, y = lower_lim, 
           label = paste("-1.96 SD =", round(lower_lim, 3)), 
           hjust = -0.1, vjust = 1.2, size = 3.5, color = "darkgreen") +
  labs(x = "Average of True and Predicted (True + Pred)/2", 
       y = "Difference (Pred - True)",
       title = "Bland-Altman Plot for Prediction Agreement") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p_ba)
ggsave("Bland_Altman_plot.pdf", plot = p_ba, width = 6, height = 5, dpi = 300)
ggsave("Bland_Altman_plot.png", plot = p_ba, width = 6, height = 5, dpi = 300)

cat("\n===== Bland-Altman Statistics =====\n")
cat(sprintf("Mean difference (bias) = %.4f\n", mean_diff))
cat(sprintf("SD of differences      = %.4f\n", sd_diff))
cat(sprintf("95%% limits of agreement: [%.4f, %.4f]\n", lower_lim, upper_lim))
cat("如果所有点都落在两条虚线之间（或绝大多数），则认为预测与真实值具有良好的一致性。\n")

cat("\n全部分析完成！生成的文件包括：\n",
    "scatter_plot.pdf/png, KM_curve.pdf/png, time_ROC.pdf,\n",
    "multivariate_cox_results.txt, merged_survival_data.csv, Table5.csv,\n",
    "residual_distribution_hist.pdf/png, Bland_Altman_plot.pdf/png\n")

