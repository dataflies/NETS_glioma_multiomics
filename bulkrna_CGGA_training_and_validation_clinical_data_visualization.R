# ==========================================================
# 组合热图（含基因集评分行）
# ==========================================================

library(ComplexHeatmap)
library(circlize)
library(openxlsx)
library(dplyr)
library(grid)

# ========================= 1. 读取数据并转换空字符串为 NA =========================
cat("读取数据...\n")
clinical <- read.csv(
  "E:/NETS/bulk rna数据/原数据/2_CGGA_172_Filtered/CGGA_172_clinical_data.csv",
  header = TRUE, stringsAsFactors = FALSE, check.names = FALSE
)
clinical <- clinical %>% mutate(across(where(is.character), ~ ifelse(. == "", NA, .)))

colnames(clinical)[colnames(clinical) == "Censor (alive=0; dead=1)"] <- "Censor"
colnames(clinical)[colnames(clinical) == "Radio_status (treated=1;un-treated=0)"] <- "Radio_status"
colnames(clinical)[colnames(clinical) == "Chemo_status (TMZ treated=1;un-treated=0)"] <- "Chemo_status"

# 读取风险评分
risk_118 <- read.csv("E:/NETS/1_CGGA_118_risk_scores.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(risk_118)[1] <- "Sample"
risk_118 <- risk_118[, c("Sample", "risk_raw")]

risk_54 <- read.csv("E:/NETS/1_CGGA_54_risk_scores.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(risk_54)[1] <- "Sample"
risk_54 <- risk_54[, c("Sample", "risk_raw")]

# 读取 CE 增强数据
ce_data <- read.xlsx("E:/NETS/171.xlsx", sheet = 1)
colnames(ce_data)[1] <- "Sample"
ce_data <- ce_data[, c("Sample", "qh")]
ce_data$qh <- ifelse(ce_data$qh == "", NA, ce_data$qh)

# --------------------- 新增：读取基因集评分 ---------------------
# 假设您之前导出了两个 Excel 文件，列名为 Sample 和 UCell_score（或 gene_score）
score_118 <- read.xlsx("E:/NETS/GeneSignature_Score_118.xlsx")
score_54  <- read.xlsx("E:/NETS/GeneSignature_Score_54.xlsx")

# 统一列名，便于合并
colnames(score_118)[colnames(score_118) == "UCell_score"] <- "gene_score"
colnames(score_54)[colnames(score_54) == "UCell_score"] <- "gene_score"
# 如果您使用其他名称，请相应调整，确保列名为 "gene_score"

# ========================= 2. 整合数据 =========================
clinical$Batch <- ifelse(clinical$Batch == 1, "mRNAseq_693", 
                         ifelse(clinical$Batch == 2, "mRNAseq_325", NA))

# 依次合并风险评分、CE、基因集评分
clinical <- clinical %>%
  left_join(risk_118, by = "Sample") %>%
  left_join(risk_54, by = "Sample") %>%
  mutate(risk_raw = coalesce(risk_raw.x, risk_raw.y)) %>%
  select(-risk_raw.x, -risk_raw.y) %>%
  left_join(ce_data, by = "Sample") %>%
  left_join(score_118, by = "Sample") %>%
  left_join(score_54, by = "Sample") %>%
  mutate(gene_score = coalesce(gene_score.x, gene_score.y)) %>%
  select(-gene_score.x, -gene_score.y)

clinical <- clinical %>% mutate(across(where(is.character), ~ ifelse(. == "", NA, .)))

data_693 <- clinical %>% filter(Batch == "mRNAseq_693") %>% arrange(risk_raw)
data_325 <- clinical %>% filter(Batch == "mRNAseq_325") %>% arrange(risk_raw)

cat("mRNAseq_693 样本数:", nrow(data_693), "\n")
cat("mRNAseq_325 样本数:", nrow(data_325), "\n")

# ========================= 3. 构建单个热图的函数（已包含 gene_score 行） =========================
build_heatmap <- function(data, cohort_name) {
  # 新增 "gene_score" 到 features 和 feature_labels
  features <- c("risk_group", "risk_raw", "gene_score", "OS", "Age", "Gender", "Grade", 
                "1p19q_codeletion_status", "IDH_mutation_status", 
                "MGMTp_methylation_status", "qh")
  feature_labels <- c("Risk Stratification", "RiskScore", "Gene Signature Score", "OS", "Age", "Gender", 
                      "Grade", "1p/19q codeletion", "IDH", "MGMTp", "CE enhancement")
  
  med <- median(data$risk_raw, na.rm = TRUE)
  data$risk_group <- ifelse(data$risk_raw <= med, "Low", "High")
  
  mat <- matrix(NA, nrow = length(features), ncol = nrow(data))
  rownames(mat) <- features
  colnames(mat) <- data$Sample
  
  for (i in seq_along(features)) {
    f <- features[i]
    vals <- data[[f]]
    if (f == "risk_group") {
      mat[i, ] <- ifelse(vals == "Low", 0, 1)
    } else if (f == "Gender") {
      mat[i, ] <- ifelse(vals == "Female", 0, 1)
    } else if (f == "Grade") {
      mat[i, ] <- case_when(
        vals == "WHO II" ~ 0,
        vals == "WHO III" ~ 1,
        vals == "WHO IV" ~ 2,
        TRUE ~ NA_real_
      )
    } else if (f == "1p19q_codeletion_status") {
      mat[i, ] <- ifelse(is.na(vals), NA, ifelse(vals == "Codel", 1, 0))
    } else if (f == "IDH_mutation_status") {
      mat[i, ] <- ifelse(vals == "Wildtype", 0, 1)
    } else if (f == "MGMTp_methylation_status") {
      mat[i, ] <- ifelse(vals == "un-methylated", 0, 1)
    } else if (f == "qh") {
      mat[i, ] <- ifelse(is.na(vals), NA, ifelse(vals %in% c("黑", "极低含量"), 0, 1))
    } else if (f %in% c("OS", "risk_raw", "Age", "gene_score")) {
      mat[i, ] <- as.numeric(vals)
    }
  }
  
  # ---- 颜色映射 ----
  col_risk_cont <- colorRamp2(c(min(mat["risk_raw", ], na.rm=TRUE),
                                median(mat["risk_raw", ], na.rm=TRUE),
                                max(mat["risk_raw", ], na.rm=TRUE)),
                              c("#2C7BB6", "white", "#D7191C"))
  col_os <- colorRamp2(c(min(mat["OS", ], na.rm=TRUE),
                         median(mat["OS", ], na.rm=TRUE),
                         max(mat["OS", ], na.rm=TRUE)),
                       c("#7B3294", "white", "#FDB863"))
  col_age <- colorRamp2(c(min(mat["Age", ], na.rm=TRUE),
                          median(mat["Age", ], na.rm=TRUE),
                          max(mat["Age", ], na.rm=TRUE)),
                        c("#1A9B9B", "white", "#E78AC3"))
  # 新增基因集评分的颜色映射（可自定义颜色）
  col_gene_score <- colorRamp2(c(min(mat["gene_score", ], na.rm=TRUE),
                                 median(mat["gene_score", ], na.rm=TRUE),
                                 max(mat["gene_score", ], na.rm=TRUE)),
                               c("#2C7BB6", "white", "#D7191C"))
  
  col_risk_bin <- c("0" = "#2C7BB6", "1" = "#D7191C")
  col_gender <- c("0" = "#E78AC3", "1" = "#5D8DA8")
  col_grade <- c("0" = "#A6D96A", "1" = "#FDAE61", "2" = "#D73027")
  col_1p19q <- c("0" = "#A0C4E8", "1" = "#3F6A9B")
  col_idh <- c("0" = "#BDBDBD", "1" = "#AB82FF")
  col_mgmtp <- c("0" = "#FED976", "1" = "#4E9A8E")
  col_ce <- c("0" = "#E0E0E0", "1" = "#C51B8A")
  
  row_color_func <- list(
    risk_group = function(x) col_risk_bin[as.character(round(x))],
    risk_raw = function(x) col_risk_cont(x),
    gene_score = function(x) col_gene_score(x),
    OS = function(x) col_os(x),
    Age = function(x) col_age(x),
    Gender = function(x) col_gender[as.character(round(x))],
    Grade = function(x) col_grade[as.character(round(x))],
    `1p19q_codeletion_status` = function(x) col_1p19q[as.character(round(x))],
    IDH_mutation_status = function(x) col_idh[as.character(round(x))],
    MGMTp_methylation_status = function(x) col_mgmtp[as.character(round(x))],
    qh = function(x) col_ce[as.character(round(x))]
  )
  
  # 预计算颜色矩阵
  cell_colors <- matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
  for (i in 1:nrow(mat)) {
    row_name <- rownames(mat)[i]
    for (j in 1:ncol(mat)) {
      val <- mat[i, j]
      if (is.na(val)) {
        cell_colors[i, j] <- NA
      } else {
        cell_colors[i, j] <- row_color_func[[row_name]](val)
      }
    }
  }
  
  # 提取连续变量的范围（用于图例真实数值）
  range_risk <- range(mat["risk_raw", ], na.rm = TRUE)
  range_os   <- range(mat["OS", ], na.rm = TRUE)
  range_age  <- range(mat["Age", ], na.rm = TRUE)
  range_gene_score <- range(mat["gene_score", ], na.rm = TRUE)
  
  # ---- 构建热图（行间距加大） ----
  row_gap_val <- 0.3
  row_height <- 1.0
  n_row <- nrow(mat)
  ht_height <- unit(n_row * row_height + (n_row - 1) * row_gap_val, "cm")
  
  ht <- Heatmap(
    matrix = mat,
    name = "Features",
    col = colorRamp2(c(0, 1), c("white", "white")),
    rect_gp = gpar(type = "none"),
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_labels = feature_labels,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 10, fontface = "bold"),
    row_title = NULL,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_heatmap_legend = FALSE,
    width = unit(ncol(mat) * 0.25 + (ncol(mat)-1) * 0.1, "cm"),
    height = ht_height,
    border = FALSE,
    row_split = 1:n_row,
    row_gap = unit(row_gap_val, "cm"),
    column_gap = unit(0.1, "cm"),
    layer_fun = function(j, i, x, y, w, h, fill) {
      for (ri in seq_along(unique(i))) {
        row_idx <- unique(i)[ri]
        row_name <- rownames(mat)[row_idx]
        for (cj in seq_along(unique(j))) {
          col_idx <- unique(j)[cj]
          color <- cell_colors[row_idx, col_idx]
          if (is.na(color)) next
          if (row_name == "risk_group") {
            grid.rect(x = x[cj], y = y[ri], width = w[cj], height = h[ri],
                      gp = gpar(fill = color, col = NA))
          } else {
            grid.rect(x = x[cj], y = y[ri], width = w[cj], height = h[ri],
                      gp = gpar(fill = color, col = "black", lwd = 0.3))
          }
        }
      }
    }
  )
  
  list(ht = ht, 
       col_risk_cont = col_risk_cont,
       col_os = col_os,
       col_age = col_age,
       col_gene_score = col_gene_score,
       col_risk_bin = col_risk_bin,
       col_gender = col_gender,
       col_grade = col_grade,
       col_1p19q = col_1p19q,
       col_idh = col_idh,
       col_mgmtp = col_mgmtp,
       col_ce = col_ce,
       range_risk = range_risk,
       range_os = range_os,
       range_age = range_age,
       range_gene_score = range_gene_score)
}

# ========================= 4. 构建两个热图 =========================
res_693 <- build_heatmap(data_693, "mRNAseq_693")
res_325 <- build_heatmap(data_325, "mRNAseq_325")
ht_693 <- res_693$ht
ht_325 <- res_325$ht

# ========================= 5. 图例（连续变量显示真实数值且缩短轴） =========================
fmt <- function(x) round(x, 2)

legend_list <- list(
  Legend(title = "Risk Stratification", labels = c("Low", "High"), 
         legend_gp = gpar(fill = c("#2C7BB6", "#D7191C"))),
  Legend(title = "RiskScore", col_fun = res_693$col_risk_cont, direction = "vertical",
         at = res_693$range_risk, labels = fmt(res_693$range_risk),
         legend_height = unit(1, "cm")),
  # 新增基因集评分图例
  Legend(title = "Gene Signature Score", col_fun = res_693$col_gene_score, direction = "vertical",
         at = res_693$range_gene_score, labels = fmt(res_693$range_gene_score),
         legend_height = unit(1, "cm")),
  Legend(title = "OS", col_fun = res_693$col_os, direction = "vertical",
         at = res_693$range_os, labels = fmt(res_693$range_os),
         legend_height = unit(1, "cm")),
  Legend(title = "Age", col_fun = res_693$col_age, direction = "vertical",
         at = res_693$range_age, labels = fmt(res_693$range_age),
         legend_height = unit(1, "cm")),
  Legend(title = "Gender", labels = c("Female", "Male"),
         legend_gp = gpar(fill = c("#E78AC3", "#5D8DA8"))),
  Legend(title = "Grade", labels = c("WHO II", "WHO III", "WHO IV"),
         legend_gp = gpar(fill = c("#A6D96A", "#FDAE61", "#D73027"))),
  Legend(title = "1p/19q codeletion", labels = c("Non-codel", "Codel"),
         legend_gp = gpar(fill = c("#A0C4E8", "#3F6A9B"))),
  Legend(title = "IDH", labels = c("Wildtype", "Mutant"),
         legend_gp = gpar(fill = c("#BDBDBD", "#AB82FF"))),
  Legend(title = "MGMTp", labels = c("unmethylated", "methylated"),
         legend_gp = gpar(fill = c("#FED976", "#4E9A8E"))),
  Legend(title = "CE enhancement", labels = c("No", "Yes"),
         legend_gp = gpar(fill = c("#E0E0E0", "#C51B8A")))
)
legend_pack <- packLegend(list = legend_list, direction = "vertical", gap = unit(0.15, "cm"))

# ========================= 6. 绘制组合图 =========================
cat("生成组合热图（mRNAseq_693 + mRNAseq_325）...\n")
pdf("Heatmap_Combined_WithGeneScore.pdf", width = 40, height = 24)

draw(ht_693 + ht_325, 
     ht_gap = unit(2, "cm"),
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     legend_gap = unit(0.15, "cm"))

pushViewport(viewport(x = unit(0.95, "npc"), y = unit(0.5, "npc"),
                      width = unit(0.12, "npc"), height = unit(0.9, "npc")))
draw(legend_pack, x = unit(0, "npc"), y = unit(0.5, "npc"), just = "center")
upViewport()

grid.text("mRNAseq_693 cohort (n=118)", 
          x = unit(0.25, "npc"), y = unit(0.98, "npc"),
          gp = gpar(fontsize = 16, fontface = "bold"))
grid.text("mRNAseq_325 cohort (n=54)", 
          x = unit(0.72, "npc"), y = unit(0.98, "npc"),
          gp = gpar(fontsize = 16, fontface = "bold"))

dev.off()
cat("组合热图已保存为 Heatmap_Combined_WithGeneScore.pdf\n")
