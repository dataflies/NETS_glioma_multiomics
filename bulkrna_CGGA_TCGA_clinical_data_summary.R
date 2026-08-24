# ==========================================================
# Table 1: 临床基线特征表（三队列：mRNAseq_693、mRNAseq_325、TCGA）
# 数据: CGGA_172_clinical_data.csv + TCGA_clinical_with_risk.csv
# 输出: Excel (.xlsx) 和 CSV (.csv)
# ==========================================================

library(dplyr)
library(openxlsx)

# ========================= 1. 读取数据 =========================
# CGGA数据
cgga <- read.csv(
  "E:/NETS/bulk rna数据/原数据/2_CGGA_172_Filtered/CGGA_172_clinical_data.csv",
  header = TRUE, stringsAsFactors = FALSE, check.names = FALSE
)

# TCGA数据
tcga <- read.csv(
  "E:/NETS/bulk&影像 TCGA外部验证/原数据/1_TCGA_Filtered_Data/1_TCGA_TCGA_clinical_with_risk.csv",
  header = TRUE, stringsAsFactors = FALSE, check.names = FALSE
)

# ========================= 2. 数据预处理 =========================
# ---- CGGA预处理 ----
# 重命名列
colnames(cgga)[colnames(cgga) == "Censor (alive=0; dead=1)"] <- "Censor"
colnames(cgga)[colnames(cgga) == "Radio_status (treated=1;un-treated=0)"] <- "Radio_status"
colnames(cgga)[colnames(cgga) == "Chemo_status (TMZ treated=1;un-treated=0)"] <- "Chemo_status"

# 定义CGGA分组
cgga$Group <- ifelse(cgga$Batch == 1, "mRNAseq_693", 
                     ifelse(cgga$Batch == 2, "mRNAseq_325", NA))
cgga$Group <- factor(cgga$Group, levels = c("mRNAseq_693", "mRNAseq_325"))

# ---- TCGA预处理 ----
# 保留所需变量（与CGGA对齐）
tcga <- tcga[, c("Patient_ID", "Histology", "Grade", "Gender", "Age", 
                 "Censor", "IDH_mutation_status", "1p19q_codeletion_status")]
# 添加缺失的变量（TCGA无放疗/化疗/MGMTp信息）
tcga$Radio_status <- NA
tcga$Chemo_status <- NA
tcga$MGMTp_methylation_status <- NA

# 重命名列以统一
colnames(tcga)[colnames(tcga) == "Patient_ID"] <- "Patient_ID"  # 保留
# 注意：TCGA的Grade可能为 "Grade II", "Grade III", "Grade IV" 等，需映射为CGGA格式
tcga$Grade <- ifelse(grepl("II", tcga$Grade), "WHO II",
                     ifelse(grepl("III", tcga$Grade), "WHO III",
                            ifelse(grepl("IV", tcga$Grade), "WHO IV", NA)))
# 将TCGA的Histology转换为二元变量（GBM vs Non-GBM）
tcga$Histology_binary <- ifelse(tcga$Histology == "GBM", "GBM", "Non-GBM")

# 添加TCGA分组
tcga$Group <- "TCGA"

# ---- 合并数据集 ----
# 统一列顺序，并选取共同变量
common_vars <- c("Age", "Histology_binary", "Grade", "Gender", "Radio_status",
                 "Chemo_status", "IDH_mutation_status", "1p19q_codeletion_status",
                 "MGMTp_methylation_status", "Censor", "Group")
# CGGA需要添加Histology_binary（已有），并选取所需列
cgga$Histology_binary <- ifelse(cgga$Histology == "GBM", "GBM", "Non-GBM")
cgga_sub <- cgga[, common_vars]
tcga_sub <- tcga[, common_vars]

data <- rbind(cgga_sub, tcga_sub)

# 将字符型变量转换为因子（固定顺序）
data$Gender <- factor(data$Gender, levels = c("Female", "Male"))
data$Censor <- factor(data$Censor, levels = c(0, 1), labels = c("Alive", "Dead"))
data$Radio_status <- factor(data$Radio_status, levels = c(0, 1), labels = c("No", "Yes"))
data$Chemo_status <- factor(data$Chemo_status, levels = c(0, 1), labels = c("No", "Yes"))
data$IDH_mutation_status <- factor(data$IDH_mutation_status, levels = c("Wildtype", "Mutant"))
data$`1p19q_codeletion_status` <- factor(data$`1p19q_codeletion_status`, 
                                         levels = c("Non-codel", "Codel"))
data$MGMTp_methylation_status <- factor(data$MGMTp_methylation_status,
                                        levels = c("un-methylated", "methylated"))
data$Grade <- factor(data$Grade, levels = c("WHO II", "WHO III", "WHO IV"))
data$Histology_binary <- factor(data$Histology_binary, levels = c("Non-GBM", "GBM"))

# 将Group转换为因子，保持顺序
data$Group <- factor(data$Group, levels = c("mRNAseq_693", "mRNAseq_325", "TCGA"))

# ========================= 3. 统计函数（通用化） =========================
# 格式化连续变量 (均值±标准差)
fmt_cont <- function(x) {
  x <- x[!is.na(x)]
  if(length(x) == 0) return("NA")
  paste0(round(mean(x, na.rm = TRUE), 2), " ± ", round(sd(x, na.rm = TRUE), 2))
}

# 生成表格行（支持任意组数，P值计算仅针对前两组）
generate_rows <- function(var_name, data, group_var, compare_groups = c("mRNAseq_693", "mRNAseq_325")) {
  x <- data[[var_name]]
  if(is.null(x)) return(NULL)
  
  groups <- levels(group_var)
  # 列名：ALL + 各组
  col_names <- c("ALL", groups)
  
  # 如果是连续变量
  if(!is.factor(x)) {
    all_stat <- fmt_cont(x)
    stats <- sapply(groups, function(g) {
      idx <- group_var == g & !is.na(group_var)
      fmt_cont(x[idx])
    })
    # P值：比较前两组
    p <- NA
    if(length(compare_groups) == 2) {
      grp1 <- group_var == compare_groups[1] & !is.na(group_var)
      grp2 <- group_var == compare_groups[2] & !is.na(group_var)
      if(sum(grp1) > 1 & sum(grp2) > 1) {
        p <- tryCatch(t.test(x[grp1], x[grp2])$p.value, error = function(e) NA)
        if(!is.na(p)) p <- ifelse(p < 0.001, "<0.001", round(p, 3))
      }
    }
    df <- data.frame(
      Feature = var_name,
      SubFeature = "",
      ALL = all_stat,
      t(stats),
      P_value = p,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    # 重命名列
    colnames(df) <- c("Feature", "SubFeature", "ALL", groups, "P_value")
    return(df)
  }
  
  # ---- 分类变量 ----
  levels_original <- levels(x)
  
  # 初始化统计列表
  all_stats <- list()
  for (grp in c("ALL", groups)) {
    all_stats[[grp]] <- list()
    for (lev in c(levels_original, "Missing")) {
      all_stats[[grp]][[lev]] <- "0(0.0%)"
    }
  }
  
  # 统计每个分组
  for (grp in groups) {
    idx <- group_var == grp & !is.na(group_var)
    x_grp <- x[idx]
    total <- sum(idx)
    valid <- sum(!is.na(x_grp))
    tbl <- table(x_grp, useNA = "no")
    
    for (lev in levels_original) {
      count <- ifelse(lev %in% names(tbl), tbl[lev], 0)
      pct <- ifelse(valid > 0, round(count / valid * 100, 1), 0)
      all_stats[[grp]][[lev]] <- paste0(count, "(", pct, "%)")
    }
    missing_count <- sum(is.na(x_grp))
    missing_pct <- ifelse(total > 0, round(missing_count / total * 100, 1), 0)
    all_stats[[grp]][["Missing"]] <- paste0(missing_count, "(", missing_pct, "%)")
  }
  
  # ---- ALL列 ----
  x_all <- x
  total_all <- length(x_all)
  valid_all <- sum(!is.na(x_all))
  tbl_all <- table(x_all, useNA = "no")
  for (lev in levels_original) {
    count <- ifelse(lev %in% names(tbl_all), tbl_all[lev], 0)
    pct <- ifelse(valid_all > 0, round(count / valid_all * 100, 1), 0)
    all_stats[["ALL"]][[lev]] <- paste0(count, "(", pct, "%)")
  }
  missing_count_all <- sum(is.na(x_all))
  missing_pct_all <- ifelse(total_all > 0, round(missing_count_all / total_all * 100, 1), 0)
  all_stats[["ALL"]][["Missing"]] <- paste0(missing_count_all, "(", missing_pct_all, "%)")
  
  # ---- 组装数据框 ----
  df_rows <- list()
  # 第一行：特征名
  first_row <- data.frame(
    Feature = var_name,
    SubFeature = "",
    ALL = "",
    stringsAsFactors = FALSE
  )
  for (g in groups) {
    first_row[[g]] <- ""
  }
  first_row$P_value <- NA
  
  df_rows[[1]] <- first_row
  
  # 子行
  for (lev in c(levels_original, "Missing")) {
    row <- data.frame(
      Feature = "",
      SubFeature = lev,
      ALL = all_stats[["ALL"]][[lev]],
      stringsAsFactors = FALSE
    )
    for (g in groups) {
      row[[g]] <- all_stats[[g]][[lev]]
    }
    row$P_value <- NA
    df_rows[[length(df_rows) + 1]] <- row
  }
  
  df <- do.call(rbind, df_rows)
  
  # ---- 计算P值（仅前两组） ----
  if(length(compare_groups) == 2) {
    grp1 <- group_var == compare_groups[1] & !is.na(group_var)
    grp2 <- group_var == compare_groups[2] & !is.na(group_var)
    if(sum(grp1) > 1 & sum(grp2) > 1) {
      x1 <- x[grp1]
      x2 <- x[grp2]
      x1 <- x1[!is.na(x1)]
      x2 <- x2[!is.na(x2)]
      if(length(x1) > 0 & length(x2) > 0) {
        tbl <- table(c(x1, x2), rep(1:2, c(length(x1), length(x2))))
        if(any(chisq.test(tbl, simulate.p.value = TRUE)$expected < 5)) {
          p <- tryCatch(fisher.test(tbl, simulate.p.value = TRUE)$p.value, error = function(e) NA)
        } else {
          p <- tryCatch(chisq.test(tbl)$p.value, error = function(e) NA)
        }
        if(!is.na(p)) {
          df[1, "P_value"] <- ifelse(p < 0.001, "<0.001", round(p, 3))
        }
      }
    }
  }
  
  return(df)
}

# ========================= 4. 定义变量列表 =========================
var_list <- list(
  Age = "Age",
  Histology = "Histology_binary",
  Grade = "Grade",
  Gender = "Gender",
  Radiotherapy = "Radio_status",
  Chemotherapy = "Chemo_status",
  IDH = "IDH_mutation_status",
  `1p19q` = "1p19q_codeletion_status",
  MGMTp = "MGMTp_methylation_status"
)

var_labels <- c(
  "Age (years)",
  "Histology",
  "Grade",
  "Gender",
  "Radiotherapy",
  "Chemotherapy",
  "IDH",
  "1p19q codeletion",
  "MGMTp methylation"
)

# ========================= 5. 构建表格 =========================
# 计算各组样本量
n_all <- nrow(data)
n_693 <- sum(data$Group == "mRNAseq_693", na.rm = TRUE)
n_325 <- sum(data$Group == "mRNAseq_325", na.rm = TRUE)
n_tcga <- sum(data$Group == "TCGA", na.rm = TRUE)

# 初始化表格
final_table <- data.frame(
  Feature = character(), SubFeature = character(),
  ALL = character(), mRNAseq_693 = character(),
  mRNAseq_325 = character(), TCGA = character(),
  P_value = character(),
  stringsAsFactors = FALSE
)

# 添加样本量行
header_row <- data.frame(
  Feature = "Number of patients",
  SubFeature = "",
  ALL = as.character(n_all),
  mRNAseq_693 = as.character(n_693),
  mRNAseq_325 = as.character(n_325),
  TCGA = as.character(n_tcga),
  P_value = "",
  stringsAsFactors = FALSE
)
final_table <- rbind(final_table, header_row)

# 循环变量
for (i in seq_along(var_list)) {
  var_name <- var_list[[i]]
  label <- var_labels[i]
  df_var <- generate_rows(var_name, data, data$Group)
  if(!is.null(df_var) && nrow(df_var) > 0) {
    df_var[1, "Feature"] <- label
    final_table <- rbind(final_table, df_var)
  } else {
    cat("警告: 变量", label, "返回空数据框，已跳过\n")
  }
}

# ========================= 6. 导出 =========================
# 动态设置列名（含样本量）
colnames(final_table) <- c("Feature", "", "ALL", 
                           paste0("mRNAseq_693 (n=", n_693, ")"),
                           paste0("mRNAseq_325 (n=", n_325, ")"),
                           paste0("TCGA (n=", n_tcga, ")"),
                           "P-value")

write.xlsx(final_table, "Table1_Baseline_with_TCGA.xlsx", rowNames = FALSE)
write.csv(final_table, "Table1_Baseline_with_TCGA.csv", row.names = FALSE, na = "")

print(final_table)
cat("\nTable 1 已保存为 Table1_Baseline_with_TCGA.xlsx 和 .csv\n")
