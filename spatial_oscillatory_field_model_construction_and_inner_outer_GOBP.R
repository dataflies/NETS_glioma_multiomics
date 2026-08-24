# ============================================================
#  函数：识别NETs三重合病灶 + 红色六边形可视化
#  输入：Seurat对象（含三种NETs评分 + 空间坐标）
#  输出：更新后的Seurat对象 + 保存PDF图片
# ============================================================
# ============================================================
#  病灶检测最终版（三重合回归，min_lesion_size = 1）
#  支持单个spot作为独立病灶
# ============================================================

detect_lesions_final <- function(obj, 
                                 percentile_threshold = 0.85,
                                 sample_name = NULL,
                                 output_dir = NULL) {
  
  # ---- 1. 提取坐标 ----
  coords <- GetTissueCoordinates(obj, scale = "lowres")
  if (is.null(coords) || nrow(coords) == 0) {
    coords <- GetTissueCoordinates(obj, scale = "tissue")
  }
  if (is.null(coords)) stop("无法提取坐标")
  
  if (all(c("imagerow", "imagecol") %in% colnames(coords))) {
    colnames(coords)[1:2] <- c("x", "y")
  } else if (all(c("row", "col") %in% colnames(coords))) {
    colnames(coords)[1:2] <- c("x", "y")
  } else {
    coords <- coords[, 1:2]
    colnames(coords) <- c("x", "y")
  }
  
  # ---- 2. 自动查找NETs评分列 ----
  all_cols <- colnames(obj@meta.data)
  col_classic <- grep("Classic", all_cols, value = TRUE, ignore.case = TRUE)
  col_de <- grep("DE_204", all_cols, value = TRUE, ignore.case = TRUE)
  col_nmr <- grep("NMR", all_cols, value = TRUE, ignore.case = TRUE)
  
  score_cols <- c(
    if (length(col_classic) > 0) col_classic[1] else NULL,
    if (length(col_de) > 0) col_de[1] else NULL,
    if (length(col_nmr) > 0) col_nmr[1] else NULL
  )
  
  if (length(score_cols) < 3) {
    warning("找到的NETs评分列不足3个，跳过样本")
    obj$lesion_id <- NA
    obj$lesion_size <- NA
    return(obj)
  }
  
  cat("   阈值:", percentile_threshold, "\n")
  cat("   使用评分列:", paste(score_cols, collapse = ", "), "\n")
  
  # ---- 3. 三重合候选（3种评分都高，Top 15%） ----
  meta <- obj@meta.data
  is_high_list <- list()
  for (col in score_cols) {
    threshold <- quantile(meta[[col]], probs = percentile_threshold, na.rm = TRUE)
    is_high_list[[col]] <- meta[[col]] >= threshold
    cat("   ", col, "阈值:", round(threshold, 3), "\n")
  }
  is_triple_high <- Reduce(`&`, is_high_list)
  
  coords$is_candidate <- is_triple_high
  coords$spot_id <- rownames(coords)
  
  candidate_spots <- coords[coords$is_candidate, ]
  cat("   三重合候选spots:", nrow(candidate_spots), "\n")
  
  if (nrow(candidate_spots) == 0) {
    cat("   无候选spots\n")
    obj$lesion_id <- NA
    obj$lesion_size <- NA
    return(obj)
  }
  
  # ---- 4. 每个候选spot独立成为一个病灶（min_lesion_size = 1） ----
  # 直接给每个候选spot分配独立的Lesion ID
  coords$lesion_id <- NA
  coords$lesion_size <- NA
  
  for (i in 1:nrow(candidate_spots)) {
    spot <- candidate_spots$spot_id[i]
    coords[spot, "lesion_id"] <- paste0("Lesion_", i)
    coords[spot, "lesion_size"] <- 1
  }
  
  obj$lesion_id <- coords[colnames(obj), "lesion_id"]
  obj$lesion_size <- coords[colnames(obj), "lesion_size"]
  
  n_lesions <- nrow(candidate_spots)
  cat("   病灶数（含单spot）:", n_lesions, "\n")
  
  # ---- 5. 病灶统计 ----
  lesion_stats <- coords %>%
    filter(!is.na(lesion_id)) %>%
    group_by(lesion_id) %>%
    summarise(
      n_spots = n(),
      center_x = median(x),
      center_y = median(y),
      .groups = "drop"
    ) %>%
    arrange(desc(n_spots))
  
  # ---- 6. 绘制六边形图 ----
  # 为了可视化清晰，把病灶按大小分级：单spot vs 多spot
  df <- coords
  df$lesion_status <- ifelse(is.na(df$lesion_id), "Other", "NETs_Core")
  
  # 区分单spot病灶和成片病灶（用颜色深浅）
  df$lesion_size_cat <- ifelse(is.na(df$lesion_id), "Other",
                               ifelse(df$lesion_size == 1, "Single_Spot", "Multi_Spot"))
  
  # 计算六边形大小
  all_coords_mat <- as.matrix(coords[, c("x", "y")])
  nn_dist <- FNN::knn.dist(all_coords_mat, k = 2)[,2]
  side <- median(nn_dist, na.rm = TRUE) / sqrt(3) * 0.85
  if (is.na(side) || side < 0.1) side <- 10
  
  angles <- (60 * (0:5)) * pi / 180
  hex_list <- list()
  for (i in 1:nrow(df)) {
    xc <- df$x[i]; yc <- df$y[i]
    vertices <- data.frame(
      x = xc + side * cos(angles),
      y = yc + side * sin(angles)
    )
    vertices$group <- i
    vertices$status <- df$lesion_size_cat[i]
    hex_list[[i]] <- vertices
  }
  hex_df <- do.call(rbind, hex_list)
  
  p <- ggplot(hex_df, aes(x = x, y = y, fill = status, group = group)) +
    geom_polygon(color = NA, size = 0) +
    coord_fixed() +
    scale_fill_manual(
      values = c("Multi_Spot" = "#E41A1C", 
                 "Single_Spot" = "#E41A1C", 
                 "Other" = "#E0E0E0"),
      name = "Region",
      labels = c("Multi-spot lesion", "Single-spot lesion", "Other")
    ) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.key.width = unit(0.6, "cm"),
      legend.key.height = unit(0.5, "cm"),
      legend.text = element_text(size = 8),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    ) +
    labs(title = paste0(sample_name, " (", n_lesions, " lesions, threshold ", percentile_threshold, ")"))
  
  if (is.null(output_dir)) {
    output_dir <- getwd()
  }
  out_name <- paste0("Lesion_Map_Final_", sample_name)
  ggsave(file.path(output_dir, paste0(out_name, ".pdf")), 
         p, width = 8, height = 7, dpi = 300, device = "pdf")
  cat("   ✅ 已保存:", file.path(output_dir, paste0(out_name, ".pdf")), "\n")
  
  cat("   📊 病灶统计:\n")
  print(lesion_stats)
  
  # 保存单spot vs 多spot统计
  cat("   📊 单spot病灶:", sum(df$lesion_size == 1, na.rm = TRUE), "\n")
  cat("   📊 多spot病灶:", sum(df$lesion_size > 1, na.rm = TRUE), "\n")
  
  return(obj)
}


# ============================================================
#  批量加载空间数据 + 生成 seurat_list（包含NETs评分）
#  替代原来读取 seurat_list_with_NETs_scores.rds 的方式
# ============================================================
# 设置输出目录
detect_dir <- "E:/NETS/空转组数据/结果/Lesion_Analysis/Lesion_Detection"
if (!dir.exists(detect_dir)) dir.create(detect_dir, recursive = TRUE)

for (s in sample_names) {
  cat("\n========================================\n")
  cat("🔄 处理样本:", s, "\n")
  
  obj <- seurat_list[[s]]
  
  # 根据样本类型选择阈值
  if (grepl("T1", s)) {
    threshold <- 0.85
  } else {
    threshold <- 0.80  # inf和nec放宽0.05
  }
  
  seurat_list[[s]] <- detect_lesions_final(
    obj,
    percentile_threshold = threshold,
    sample_name = s,
    output_dir = detect_dir
  )
}

# 保存完整数据
saveRDS(seurat_list, file = file.path(detect_dir, "seurat_list_with_lesions_final.rds"))
cat("\n💾 已保存至:", file.path(detect_dir, "seurat_list_with_lesions_final.rds"), "\n")




# ============================================================
#  震荡场绘图函数（自检版）
#  统一断点，固定颜色，红色六边形震源，完整图例
# ============================================================

plot_oscillation_field <- function(coords, 
                                   lesion_id_vec,
                                   max_radius = 5,
                                   source_value = 5,
                                   decay_type = "linear",
                                   decay_rate = 1,
                                   sigma = 2,
                                   sample_name = NULL,
                                   output_dir = NULL,
                                   global_breaks = NULL,          # 必须指定，如 c(0, 2, 5, 9, 14, 20)
                                   grade_colors = NULL,           # 长度 = length(global_breaks)-1
                                   zero_color = "#D3D3D3",
                                   source_color = "#E41A1C") {    # 震源红色
  
  # ---- 1. 坐标检查 ----
  if (is.null(coords) || nrow(coords) == 0) stop("坐标为空")
  if (!all(c("x", "y") %in% colnames(coords))) {
    colnames(coords)[1:2] <- c("x", "y")
  }
  coords$spot_id <- rownames(coords)
  
  # ---- 2. 合并lesion_id ----
  coords$lesion_id <- lesion_id_vec[rownames(coords)]
  
  # ---- 3. 震源筛选 ----
  source_spots <- coords[!is.na(coords$lesion_id), ]
  if (nrow(source_spots) == 0) {
    warning("没有震源，跳过绘图")
    return(NULL)
  }
  cat("   震源数:", nrow(source_spots), "\n")
  
  # ---- 4. 计算震荡场 ----
  all_coords <- coords[, c("x", "y")]
  source_coords <- source_spots[, c("x", "y")]
  
  nn_dist <- FNN::knn.dist(as.matrix(all_coords), k = 2)[,2]
  avg_dist <- median(nn_dist, na.rm = TRUE)
  if (is.na(avg_dist) || avg_dist < 1) avg_dist <- 10
  
  osc_values <- rep(0, nrow(coords))
  
  for (i in 1:nrow(source_coords)) {
    sx <- source_coords$x[i]
    sy <- source_coords$y[i]
    
    dist_vec <- sqrt((all_coords$x - sx)^2 + (all_coords$y - sy)^2)
    dist_grid <- dist_vec / avg_dist
    
    within_range <- dist_grid <= max_radius
    if (sum(within_range) == 0) next
    
    if (decay_type == "linear") {
      contribution <- pmax(source_value - decay_rate * dist_grid, 0)
    } else if (decay_type == "gaussian") {
      contribution <- source_value * exp(-dist_grid^2 / (2 * sigma^2))
    } else {
      stop("decay_type 必须是 'linear' 或 'gaussian'")
    }
    
    contribution[!within_range] <- 0
    osc_values <- osc_values + contribution
  }
  
  coords$osc <- osc_values
  cat("   震荡值范围: [", round(min(osc_values), 2), ", ", round(max(osc_values), 2), "]\n")
  
  # ---- 5. 六边形基底（用于所有绘图） ----
  side <- median(nn_dist, na.rm = TRUE) / sqrt(3) * 0.85
  if (is.na(side) || side < 0.1) side <- 10
  
  angles <- (60 * (0:5)) * pi / 180
  hex_list <- list()
  for (i in 1:nrow(coords)) {
    xc <- coords$x[i]; yc <- coords$y[i]
    vertices <- data.frame(
      x = xc + side * cos(angles),
      y = yc + side * sin(angles)
    )
    vertices$group <- i
    vertices$osc <- coords$osc[i]
    hex_list[[i]] <- vertices
  }
  hex_df <- do.call(rbind, hex_list)
  
  # ---- 6. 渐变图（用于参考，但震源用红色六边形） ----
  # 为方便，渐变图也添加震源六边形（红色），但更清晰的方式是在分级图里展示。
  # 渐变图可保留，但用户主要关注分级图，所以这里只生成分级图即可。
  # 我们统一生成两种图：渐变（科学展示）和分级（离散展示），但分级图必须满足所有要求。
  
  # ---- 7. 分级图（统一断点，完整图例） ----
  # 7.1 检查断点
  if (is.null(global_breaks) || length(global_breaks) < 2) {
    stop("必须提供 global_breaks（至少两个断点）！")
  }
  # 确保断点覆盖所有非零值
  non_zero_vals <- osc_values[osc_values > 0]
  if (length(non_zero_vals) > 0) {
    min_val <- min(non_zero_vals)
    max_val <- max(non_zero_vals)
    # 如果断点第一个值 > min_val，则自动扩展
    if (global_breaks[1] > min_val) global_breaks[1] <- min_val - 1e-6
    if (global_breaks[length(global_breaks)] < max_val) global_breaks[length(global_breaks)] <- max_val + 1e-6
  }
  
  # 7.2 分配等级
  # 先默认所有点为 "Zero"
  level_labels <- c("Zero", paste0("L", 1:(length(global_breaks)-1)))
  coords$grade <- factor("Zero", levels = level_labels)
  
  # 对非零值进行 cut
  if (length(non_zero_vals) > 0) {
    cut_labels <- paste0("L", 1:(length(global_breaks)-1))
    cut_result <- cut(osc_values, breaks = global_breaks, labels = cut_labels, include.lowest = TRUE)
    coords$grade[osc_values > 0] <- as.character(cut_result[osc_values > 0])
  }
  
  # 确保 grade 是因子，包含所有可能级别（即使某级别缺失）
  coords$grade <- factor(coords$grade, levels = level_labels)
  
  # 7.3 颜色配置
  if (is.null(grade_colors)) {
    # 默认5级：橙黄绿青蓝
    default_colors <- c("#FF8C00", "#FFD700", "#32CD32", "#00CED1", "#4169E1")
    # 确保颜色数量与断点区间数匹配
    grade_colors <- default_colors[1:(length(global_breaks)-1)]
  }
  # 构建完整颜色向量（Zero + 各级）
  all_colors <- c(zero_color, grade_colors)
  names(all_colors) <- level_labels
  
  # 7.4 绘图数据
  hex_df_grade <- hex_df
  hex_df_grade$grade <- coords$grade[hex_df_grade$group]
  
  # ---- 8. 生成震源六边形多边形 ----
  # 对每个震源，生成一个红色六边形（覆盖在图上）
  source_hex_list <- list()
  if (nrow(source_spots) > 0) {
    source_side <- side * 0.9  # 略小一点避免覆盖全部
    for (i in 1:nrow(source_spots)) {
      xc <- source_spots$x[i]
      yc <- source_spots$y[i]
      vertices <- data.frame(
        x = xc + source_side * cos(angles),
        y = yc + source_side * sin(angles)
      )
      vertices$group <- paste0("source_", i)
      source_hex_list[[i]] <- vertices
    }
    source_hex_df <- do.call(rbind, source_hex_list)
  } else {
    source_hex_df <- NULL
  }
  
  # ---- 9. 绘制分级图 ----
  p_grade <- ggplot(hex_df_grade, aes(x = x, y = y, fill = grade, group = group)) +
    geom_polygon(color = NA, size = 0) +
    coord_fixed() +
    # 使用 scale_fill_manual 强制显示所有级别
    scale_fill_manual(
      values = all_colors,
      name = "Oscillation\nLevel",
      drop = FALSE,              # 保留未出现的级别
      labels = level_labels
    ) +
    # 叠加震源红色六边形
    geom_polygon(data = source_hex_df, 
                 aes(x = x, y = y, group = group),
                 fill = source_color, color = NA, size = 0,
                 inherit.aes = FALSE) +
    theme_void() +
    theme(legend.position = "right",
          legend.key.size = unit(0.6, "cm"),
          legend.text = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +
    labs(title = paste0(sample_name, " - Oscillation Grades (fixed breaks)"))
  
  # ---- 10. 渐变图（非必须，但保留） ----
  p_grad <- ggplot(hex_df, aes(x = x, y = y, fill = osc, group = group)) +
    geom_polygon(color = NA, size = 0) +
    coord_fixed() +
    scale_fill_gradient2(
      low = "#FFFFFF", 
      mid = "#66CC66", 
      high = "#004D00",
      midpoint = median(osc_values[osc_values > 0], na.rm = TRUE),
      name = "Oscillation\nValue"
    ) +
    geom_polygon(data = source_hex_df, 
                 aes(x = x, y = y, group = group),
                 fill = source_color, color = NA, size = 0,
                 inherit.aes = FALSE) +
    theme_void() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +
    labs(title = paste0(sample_name, " - Oscillation Field"))
  
  # ---- 11. 保存 ----
  if (is.null(output_dir)) output_dir <- getwd()
  ggsave(file.path(output_dir, paste0("Oscillation_Gradient_", sample_name, ".pdf")),
         p_grad, width = 8, height = 7, dpi = 300, device = "pdf")
  ggsave(file.path(output_dir, paste0("Oscillation_Grade_", sample_name, ".pdf")),
         p_grade, width = 8, height = 7, dpi = 300, device = "pdf")
  
  cat("   ✅ 已保存渐变图和分级图\n")
  
  return(osc_values)
}


# ============================================================
#  批量震荡可视化
#  保存分类
# ============================================================
# 输出目录
osc_dir <- "E:/NETS/空转组数据/结果/Lesion_Analysis/Oscillation_Field"
if (!dir.exists(osc_dir)) dir.create(osc_dir, recursive = TRUE)

# 全局断点（根据你的数据调整）
global_breaks <- c(0, 2, 5, 9, 14, 20, 30)
grade_colors <- c("#FF8C00", "#FFD700", "#32CD32", "#00CED1", "#4169E1", "#8A2BE2")

for (s in names(seurat_list)) {
  cat("\n========================================\n")
  cat("🔄 计算震荡场:", s, "\n")
  
  # ---- 容错：如果某个样本出错，跳过并继续 ----
  result <- tryCatch({
    obj <- seurat_list[[s]]
    
    # 检查是否有病灶
    if (!"lesion_id" %in% colnames(obj@meta.data) || all(is.na(obj$lesion_id))) {
      cat("   ⚠️ 无病灶，跳过\n")
      next
    }
    
    # 提取坐标
    coords <- GetTissueCoordinates(obj, scale = "lowres")
    if (is.null(coords) || nrow(coords) == 0) {
      coords <- GetTissueCoordinates(obj, scale = "tissue")
    }
    if (is.null(coords)) {
      cat("   ❌ 无法提取坐标，跳过\n")
      next
    }
    
    # 统一列名
    if (all(c("imagerow", "imagecol") %in% colnames(coords))) {
      colnames(coords)[1:2] <- c("x", "y")
    } else if (all(c("row", "col") %in% colnames(coords))) {
      colnames(coords)[1:2] <- c("x", "y")
    } else {
      coords <- coords[, 1:2]
      colnames(coords) <- c("x", "y")
    }
    
    # 提取 lesion_id
    lesion_vec <- obj$lesion_id[rownames(coords)]
    
    # 调用绘图函数
    osc_vals <- plot_oscillation_field(
      coords = coords,
      lesion_id_vec = lesion_vec,
      max_radius = 5,
      source_value = 5,
      decay_type = "linear",
      decay_rate = 1,
      sample_name = s,
      output_dir = osc_dir,
      global_breaks = global_breaks,
      grade_colors = grade_colors,
      zero_color = "#D3D3D3",
      source_color = "#E41A1C"
    )
    
    # 保存震荡值
    if (!is.null(osc_vals)) {
      obj$oscillation_value <- osc_vals
      names(obj$oscillation_value) <- colnames(obj)
      seurat_list[[s]] <- obj
    }
    
    cat("   ✅ 样本", s, "处理完成\n")
    
  }, error = function(e) {
    cat("   ❌ 样本", s, "出错，跳过:", e$message, "\n")
    return(NULL)
  })
}

# 保存更新后的数据
saveRDS(seurat_list, file = file.path(osc_dir, "seurat_list_with_oscillation_final.rds"))
cat("\n✅ 全部循环结束！\n")





# ============================================================
#  折线图 + 中心丰度 + 完整表格
#  横轴：L1~L6（存在） + Center（最高等级）
#  输出：PDF折线图、CSV表格（含所有级别和中心）
# ============================================================
rm(list = ls())
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

# ----- 1. 路径设置 -----
osc_dir <- "E:/NETS/空转组数据/结果/Lesion_Analysis/Oscillation_Field"
rctd_dir <- "E:/NETS/空转组数据/结果/RCTD_Analysis"
grad_dir <- "E:/NETS/空转组数据/结果/Lesion_Analysis/Gradient_Line"
if (!dir.exists(grad_dir)) dir.create(grad_dir, recursive = TRUE)

# ----- 2. 加载数据 -----
seurat_list <- readRDS(file.path(osc_dir, "seurat_list_with_oscillation_final.rds"))
sample_names <- names(seurat_list)

# ----- 3. 补充 RCTD 权重（如果缺失）-----
cat("\n===== 检查并补充 RCTD 权重 =====\n")
for (s in sample_names) {
  obj <- seurat_list[[s]]
  if (any(grepl("^proportion_", colnames(obj@meta.data)))) {
    cat(s, "已有权重列，跳过\n")
    next
  }
  cat("处理", s, "：加载 RCTD 结果...\n")
  rctd_file <- file.path(rctd_dir, paste0("RCTD_result_", s, ".rds"))
  if (!file.exists(rctd_file)) {
    cat("   ❌ 未找到文件，跳过\n")
    next
  }
  rctd_obj <- readRDS(rctd_file)
  weights <- rctd_obj@results$weights
  if (is.null(weights) || nrow(weights) == 0) {
    cat("   ⚠️ 权重矩阵为空，跳过\n")
    next
  }
  common <- intersect(rownames(weights), colnames(obj))
  if (length(common) == 0) {
    cat("   ❌ 无匹配 barcode，跳过\n")
    next
  }
  weights <- weights[common, , drop = FALSE]
  for (ct in colnames(weights)) {
    col_name <- paste0("proportion_", ct)
    vec <- rep(NA, ncol(obj))
    names(vec) <- colnames(obj)
    vec[common] <- weights[common, ct]
    obj[[col_name]] <- vec
  }
  seurat_list[[s]] <- obj
  cat("   ✅ 添加了", ncol(weights), "种细胞类型\n")
}
saveRDS(seurat_list, file.path(osc_dir, "seurat_list_with_weights.rds"))

# ----- 4. 定义颜色（严格匹配去前缀后的细胞名）-----
cell_colors <- c(
  "Neuron"          = "#6A3D9A", 
  "Vasc"            = "#1B9E77",
  "MES_Hyp"         = "#1F78B4", 
  "Mac"             = "#7570B3",
  "Oligo"           = "#FDBF6F", 
  "MES"             = "#E41A1C",
  "Prolif_Metab"    = "#D95F02", 
  "MES_Ast"         = "#B15928",
  "Reactive_Ast"    = "#8B4513", 
  "NPC"             = "#A6CEE3",
  "InflammatoryMac" = "#E7298A", 
  "Chromatin_reg"   = "#33A02C",
  "OPC"             = "#CAB2D6", 
  "AC"              = "#2E8B57",
  "General_Tumor"   = "#000000"
)

# ----- 5. 提取数据（Center = 震源）-----
grade_levels <- c("L1", "L2", "L3", "L4", "L5", "L6", "Center")
all_plot_data <- list()
all_table_data <- list()

for (s in sample_names) {
  cat("\n处理样本:", s, "\n")
  obj <- seurat_list[[s]]
  meta <- obj@meta.data
  
  if (!"oscillation_grade" %in% colnames(meta) || !"lesion_id" %in% colnames(meta)) {
    cat("   ⚠️ 缺少必要列，跳过\n")
    next
  }
  
  prop_cols <- grep("^proportion_", colnames(meta), value = TRUE)
  if (length(prop_cols) == 0) next
  
  # 提取细胞类型名称（去掉前缀）
  cell_types <- gsub("^proportion_", "", prop_cols)
  
  df_list <- list()
  
  for (g in grade_levels) {
    if (g == "Center") {
      spots <- !is.na(meta$lesion_id)
    } else {
      spots <- meta$oscillation_grade == g & is.na(meta$lesion_id)
    }
    
    if (sum(spots, na.rm = TRUE) > 0) {
      avg <- colMeans(meta[spots, prop_cols, drop = FALSE], na.rm = TRUE)
      # 命名并转为数据框
      names(avg) <- cell_types   # 关键：去掉前缀
      df <- as.data.frame(t(avg))
      df$Grade <- g
      df_list[[g]] <- df
      cat("   ", g, ":", sum(spots, na.rm = TRUE), "spots\n")
    } else {
      cat("   ", g, ": 无 spot，跳过\n")
    }
  }
  
  if (length(df_list) == 0) {
    cat("   ❌ 无有效数据，跳过\n")
    next
  }
  
  df_all <- do.call(rbind, df_list)
  # 确保 Grade 为因子，顺序固定
  df_all$Grade <- factor(df_all$Grade, levels = grade_levels)
  df_all <- df_all[order(df_all$Grade), ]
  
  # ---- 长格式（绘图用） ----
  # 直接使用 cell_types 作为列名，pivot_longer 会得到干净的 CellType
  df_long <- pivot_longer(df_all, 
                          cols = all_of(cell_types), 
                          names_to = "CellType", 
                          values_to = "Proportion")
  df_long <- df_long[!is.na(df_long$Proportion), ]
  df_long$Sample <- s
  all_plot_data[[s]] <- df_long
  
  # ---- 宽格式表格（导出用） ----
  table_wide <- df_all
  table_wide$Sample <- s
  table_wide$IsCenter <- table_wide$Grade == "Center"
  all_table_data[[s]] <- table_wide
}

# ----- 6. 绘制折线图（含 Center）-----
for (s in names(all_plot_data)) {
  df <- all_plot_data[[s]]
  if (nrow(df) == 0) next
  
  # 确保 Grade 顺序
  df$Grade <- factor(df$Grade, levels = grade_levels)
  
  # 检查是否有未定义颜色的细胞类型
  miss_ct <- setdiff(unique(df$CellType), names(cell_colors))
  if (length(miss_ct) > 0) {
    cat("警告：以下细胞类型未定义颜色，将用灰色：", paste(miss_ct, collapse=", "), "\n")
    for (ct in miss_ct) cell_colors[[ct]] <- "grey50"
  }
  
  df$Percent <- df$Proportion * 100
  
  p <- ggplot(df, aes(x = Grade, y = Percent, color = CellType, group = CellType)) +
    geom_line(size = 1.0) +
    geom_point(size = 2.5) +
    geom_text(aes(label = sprintf("%.1f", Percent)), 
              vjust = -0.6, size = 2.8, show.legend = FALSE) +
    scale_color_manual(values = cell_colors) +
    labs(title = paste0(s, " - Cell Composition (Center = Lesion Core)"),
         x = "Region / Oscillation Grade", y = "Average Proportion (%)") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "right",
          legend.key.size = unit(0.5, "cm"),
          axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold"),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  pdf_file <- file.path(grad_dir, paste0("Gradient_Line_Center_", s, ".pdf"))
  ggsave(pdf_file, p, width = 13, height = 8, dpi = 300, device = "pdf")
  cat("   ✅ 已保存:", pdf_file, "\n")
}

# ----- 7. 导出表格（含 Center）-----
combined_table <- do.call(rbind, all_table_data)
# 百分比转换（保留两位小数）
for (col in colnames(combined_table)) {
  if (!col %in% c("Grade", "Sample", "IsCenter")) {
    combined_table[[col]] <- round(combined_table[[col]] * 100, 2)
  }
}
# 调整列顺序
cell_col_names <- setdiff(colnames(combined_table), c("Grade", "Sample", "IsCenter"))
combined_table <- combined_table[, c("Grade", "Sample", "IsCenter", cell_col_names)]
combined_table <- combined_table[order(combined_table$Sample, combined_table$Grade), ]

write.csv(combined_table, file = file.path(grad_dir, "丰度表格_含Center.csv"), row.names = FALSE)

# 保存长格式
combined_long <- do.call(rbind, all_plot_data)
saveRDS(combined_long, file = file.path(grad_dir, "gradient_line_data_withCenter.rds"))

cat("\n✅ 所有折线图生成完成！\n")
cat("📁 输出目录:", grad_dir, "\n")
cat("   - Gradient_Line_Center_*.pdf  (折线图，Center = 震源)\n")
cat("   - 丰度表格_含Center.csv  (含 L1~L6 + Center 的丰度)\n")


rm(list = ls())
library(ggplot2)
library(dplyr)
library(tidyr)

# ----- 读取数据 -----
csv_path <- "E:/NETS/空转组数据/结果/Lesion_Analysis/Gradient_Line/丰度表格_含Center.csv"
df <- read.csv(csv_path)

# ----- 筛选 L1 和 Center -----
df_l1_center <- df %>%
  filter(Grade %in% c("L1", "Center")) %>%
  # 将 Grade 转为有序因子（确保 Center 在右侧）
  mutate(Grade = factor(Grade, levels = c("L1", "Center")))

# ----- 定义样本类型（手动标注，可根据样本名规则）-----
df_l1_center <- df_l1_center %>%
  mutate(SampleType = case_when(
    grepl("T1$", Sample) ~ "Tumor Core (T1)",
    grepl("inf$", Sample) ~ "Infiltrating Edge (inf)",
    grepl("nec$", Sample) ~ "Necrotic (nec)",
    TRUE ~ "Other"
  ))

# ----- 定义感兴趣的细胞类型（可按需调整顺序）-----
cell_types_of_interest <- c("MES_Ast", "MES", "MES_Hyp", "Mac", "Neutrophil", "Vasc", "AC", "OPC", "Oligo", "InflammatoryMac")

# ----- 转为长格式，筛选目标细胞类型 -----
df_long <- df_l1_center %>%
  select(Sample, SampleType, Grade, all_of(cell_types_of_interest)) %>%
  pivot_longer(cols = all_of(cell_types_of_interest), names_to = "CellType", values_to = "Proportion")

# ----- 绘制分面折线图（L1 vs Center）-----
p <- ggplot(df_long, aes(x = Grade, y = Proportion, group = Sample, color = Sample)) +
  geom_line(size = 1.2, alpha = 0.7) +
  geom_point(size = 2.5) +
  facet_grid(CellType ~ SampleType, scales = "free_y") +
  labs(title = "L1 vs Center: Cell Composition Change",
       x = "Region", y = "Average Proportion (%)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 10),
        axis.text.x = element_text(face = "bold"))

# 保存
ggsave(file.path("E:/NETS/空转组数据/结果/Lesion_Analysis/Gradient_Line", 
                 "L1_vs_Center_Comparison.pdf"), 
       p, width = 12, height = 14, dpi = 300)



# ============================================================
#  箱线图：基于震荡分级（内/外）的细胞丰度分布
#  每个样本单独输出 PDF
# ============================================================
# ============================================================
#  箱线图：内外比较（单图，横轴为细胞类型，内/外用颜色区分）
# ============================================================
rm(list = ls())
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(writexl)   # 用于导出xlsx

# ----- 路径 -----
osc_dir <- "E:/NETS/空转组数据/结果/Lesion_Analysis/Oscillation_Field"
grad_dir <- "E:/NETS/空转组数据/结果/Lesion_Analysis/Boxplots"
if (!dir.exists(grad_dir)) dir.create(grad_dir, recursive = TRUE)

seurat_list <- readRDS(file.path(osc_dir, "seurat_list_with_weights.rds"))
sample_names <- names(seurat_list)

# ----- 绘图函数（分组箱线图，横轴细胞类型，内外并排）-----
plot_inner_outer_boxplot <- function(obj, sample_name, output_dir) {
  meta <- obj@meta.data
  
  if (!"oscillation_grade" %in% colnames(meta)) {
    cat("   ⚠️ 缺少 oscillation_grade 列，跳过\n")
    return(NULL)
  }
  prop_cols <- grep("^proportion_", colnames(meta), value = TRUE)
  if (length(prop_cols) == 0) {
    cat("   ⚠️ 缺少 proportion_* 列，跳过\n")
    return(NULL)
  }
  
  valid <- !is.na(meta$oscillation_grade) & meta$oscillation_grade != "Zero"
  if (sum(valid) == 0) {
    cat("   ⚠️ 无有效分级，跳过\n")
    return(NULL)
  }
  
  grade_vec <- meta$oscillation_grade[valid]
  exist_levels <- sort(unique(grade_vec), decreasing = TRUE)
  n <- length(exist_levels)
  if (n < 2) {
    cat("   ⚠️ 级别数少于2，无法划分内外\n")
    return(NULL)
  }
  
  inner_levels <- exist_levels[1:floor(n/2)]
  outer_levels <- exist_levels[(floor(n/2)+1):n]
  group_vec <- ifelse(grade_vec %in% inner_levels, "Inner", "Outer")
  
  # 构建长格式数据
  df_list <- list()
  for (ct in prop_cols) {
    cell_name <- gsub("proportion_", "", ct)
    prop_vals <- meta[[ct]][valid]
    df_list[[cell_name]] <- data.frame(
      CellType = cell_name,
      Group = group_vec,
      Proportion = prop_vals,
      stringsAsFactors = FALSE
    )
  }
  df_all <- do.call(rbind, df_list)
  df_all$Group <- factor(df_all$Group, levels = c("Inner", "Outer"))
  
  # 计算统计量（用于标注和导出）
  stats_df <- df_all %>%
    group_by(CellType, Group) %>%
    summarise(
      mean_prop = mean(Proportion, na.rm = TRUE) * 100,
      sd_prop = sd(Proportion, na.rm = TRUE) * 100,
      n = n(),
      .groups = "drop"
    )
  
  # 计算p值（Wilcoxon检验）
  p_values <- df_all %>%
    group_by(CellType) %>%
    summarise(
      p = tryCatch(wilcox.test(Proportion ~ Group)$p.value, error = function(e) NA),
      .groups = "drop"
    ) %>%
    mutate(
      p_label = ifelse(is.na(p), "NA",
                       ifelse(p < 0.001, "***",
                              ifelse(p < 0.01, "**",
                                     ifelse(p < 0.05, "*", "n.s.")))),
      p_text = ifelse(is.na(p), "NA", paste0("p = ", round(p, 4)))
    )
  
  # 合并统计和p值
  stats_full <- left_join(stats_df, p_values, by = "CellType")
  
  # 为了绘图，需要计算每个细胞类型内外的最大y值用于放置标签
  max_y <- df_all %>%
    group_by(CellType) %>%
    summarise(max_y = max(Proportion * 100, na.rm = TRUE) * 1.1, .groups = "drop")
  stats_full <- left_join(stats_full, max_y, by = "CellType")
  
  # ---- 绘图（分组箱线图）----
  # 颜色：内层用深色，外层用浅色（两种颜色）
  inner_color <- "#2c3e50"   # 深灰蓝
  outer_color <- "#95a5a6"   # 浅灰
  # 或者使用蓝色系：inner = "#2b83ba", outer = "#abd9e9"
  
  p <- ggplot(df_all, aes(x = CellType, y = Proportion * 100, fill = Group)) +
    geom_boxplot(outlier.size = 0.5, width = 0.7, position = position_dodge(0.8)) +
    scale_fill_manual(values = c("Inner" = inner_color, "Outer" = outer_color),
                      name = "Group") +
    labs(
      title = paste0(sample_name, " (", n, " levels, inner: ", 
                     paste(inner_levels, collapse = ", "), ")"),
      y = "Proportion (%)",
      x = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.title.y = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "top"
    )
  
  # 添加均值标注（在箱线图上方）
  # 需要计算每个箱线图的位置（x坐标）来放置文本
  # 使用 position_dodge 后，x坐标偏移量为 ±0.225（0.8宽度的一半）
  # 我们可以直接在 geom_text 中通过 position = position_dodge(0.8) 自动定位
  p_labeled <- p +
    geom_text(
      data = stats_full,
      aes(x = CellType, y = max_y, 
          label = paste0(round(mean_prop, 1), "%"),
          group = Group),
      position = position_dodge(0.8),
      size = 2.5,
      vjust = -0.5
    ) +
    # 添加p值（在细胞类型名称上方，用星号或p文本）
    geom_text(
      data = p_values,
      aes(x = CellType, y = Inf, label = p_label),
      inherit.aes = FALSE,
      size = 4,
      vjust = 2,
      hjust = 0.5
    )
  
  # 保存PDF
  pdf_file <- file.path(output_dir, paste0("Boxplot_InnerOuter_", sample_name, ".pdf"))
  ggsave(pdf_file, p_labeled, width = 14, height = 8, dpi = 300, device = "pdf")
  cat("   ✅ 已保存:", pdf_file, "\n")
  
  # 导出统计数据为xlsx
  write_xlsx(stats_full, path = file.path(output_dir, paste0("Stats_", sample_name, ".xlsx")))
  cat("   ✅ 统计数据已导出:", file.path(output_dir, paste0("Stats_", sample_name, ".xlsx")), "\n")
  
  return(list(plot = p_labeled, stats = stats_full))
}

# ----- 批量执行 -----
for (s in sample_names) {
  cat("\n========================================\n")
  cat("🔄 处理样本:", s, "\n")
  obj <- seurat_list[[s]]
  result <- plot_inner_outer_boxplot(obj, sample_name = s, output_dir = grad_dir)
}
cat("\n✅ 所有样本箱线图生成完成！\n")







# ============================================================
#  双向 GOBP 分析：Inner vs Outer
#  取 Inner-up 和 Inner-down 各 Top5 词条，分别气泡图展示
# ============================================================

rm(list = ls())
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)

# ----- 1. 路径设置 -----
osc_dir <- "E:/NETS/空转组数据/结果/Lesion_Analysis/Oscillation_Field"
go_dir <- "E:/NETS/空转组数据/结果/Lesion_Analysis/GOBP_InnerOuter"
if (!dir.exists(go_dir)) dir.create(go_dir, recursive = TRUE)

seurat_list <- readRDS(file.path(osc_dir, "seurat_list_with_weights.rds"))
sample_names <- names(seurat_list)

# ----- 2. 修正分层函数（Center强制入Inner，Zero/NA剔除）-----
get_location_fixed <- function(obj) {
  grade_vec <- obj$oscillation_grade
  lesion_vec <- obj$lesion_id
  
  location <- rep(NA, length(grade_vec))
  all_levels <- sort(unique(grade_vec[!is.na(grade_vec) & grade_vec != "Zero"]))
  if (length(all_levels) == 0) return(location)
  
  sorted_levels <- sort(all_levels, decreasing = TRUE)
  n <- length(sorted_levels)
  inner_levels <- sorted_levels[1:floor(n/2)]
  outer_levels <- sorted_levels[(floor(n/2)+1):n]
  
  for (i in seq_along(grade_vec)) {
    g <- grade_vec[i]
    if (is.na(g) || g == "Zero") {
      location[i] <- NA
    } else if (g %in% inner_levels) {
      location[i] <- "Inner"
    } else if (g %in% outer_levels) {
      location[i] <- "Outer"
    }
  }
  # Center强制入Inner
  location[!is.na(lesion_vec)] <- "Inner"
  return(location)
}

# ----- 3. 批量执行（Up + Down 双向GOBP）-----
all_go_up <- list()
all_go_down <- list()

for (s in sample_names) {
  cat("\n========================================\n")
  cat("🔄 处理样本:", s, "\n")
  obj <- seurat_list[[s]]
  
  if (!"oscillation_grade" %in% colnames(obj@meta.data)) {
    cat("   ⚠️ 缺少 oscillation_grade 列，跳过\n")
    next
  }
  
  obj$location <- get_location_fixed(obj)
  cat("   Inner spots:", sum(obj$location == "Inner", na.rm = TRUE), "\n")
  cat("   Outer spots:", sum(obj$location == "Outer", na.rm = TRUE), "\n")
  cat("   NA (剔除):", sum(is.na(obj$location)), "\n")
  
  if (sum(obj$location %in% c("Inner", "Outer"), na.rm = TRUE) < 10) {
    cat("   ⚠️ 有效spots太少，跳过\n")
    next
  }
  
  obj_sub <- subset(obj, subset = !is.na(location))
  Idents(obj_sub) <- "location"
  
  cat("   ⏳ 运行 FindMarkers (Inner vs Outer)...\n")
  markers <- FindMarkers(
    obj_sub,
    ident.1 = "Inner",
    ident.2 = "Outer",
    assay = "Spatial",
    slot = "data",
    logfc.threshold = 0.25,
    test.use = "wilcox",
    min.pct = 0.1
  )
  
  # ---- Up（Inner上调） ----
  up_genes <- rownames(markers[markers$p_val_adj < 0.05 & markers$avg_log2FC > 0, ])
  cat("   ✅ Inner-up 基因:", length(up_genes), "\n")
  
  if (length(up_genes) >= 5) {
    ego_up <- enrichGO(up_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                       pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
    if (!is.null(ego_up) && nrow(ego_up@result) > 0) {
      go_up <- as.data.frame(ego_up@result)
      go_up$Sample <- s
      go_up$Direction <- "Inner_up"
      all_go_up[[s]] <- go_up
      cat("      ✅ GOBP 通路:", nrow(go_up), "\n")
    }
  }
  
  # ---- Down（Inner下调 = Outer上调） ----
  down_genes <- rownames(markers[markers$p_val_adj < 0.05 & markers$avg_log2FC < 0, ])
  cat("   ✅ Inner-down 基因:", length(down_genes), "\n")
  
  if (length(down_genes) >= 5) {
    ego_down <- enrichGO(down_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                         pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
    if (!is.null(ego_down) && nrow(ego_down@result) > 0) {
      go_down <- as.data.frame(ego_down@result)
      go_down$Sample <- s
      go_down$Direction <- "Inner_down"
      all_go_down[[s]] <- go_down
      cat("      ✅ GOBP 通路:", nrow(go_down), "\n")
    }
  }
}

# ----- 4. 合并 Up 和 Down 结果 -----
go_up_all <- do.call(rbind, all_go_up)
go_down_all <- do.call(rbind, all_go_down)

# ----- 5. 提取每个样本每个方向的 Top5 -----
extract_top5 <- function(df, direction_label) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df %>%
    group_by(Sample) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice_head(n = 5) %>%
    ungroup() %>%
    mutate(
      Direction = direction_label,
      # 截断过长的词条名
      Description_short = ifelse(nchar(Description) > 50,
                                 paste0(substr(Description, 1, 47), "..."),
                                 Description)
    )
}

top5_up <- extract_top5(go_up_all, "Inner_up")
top5_down <- extract_top5(go_down_all, "Inner_down")

# 合并，并确保每个样本都有方向标签
plot_data <- bind_rows(top5_up, top5_down)

if (is.null(plot_data) || nrow(plot_data) == 0) {
  cat("⚠️ 无足够数据绘图\n")
} else {
  # 确保样本顺序
  plot_data$Sample <- factor(plot_data$Sample, levels = sample_names)
  plot_data$Direction <- factor(plot_data$Direction, levels = c("Inner_up", "Inner_down"))
  
  # ---- 6. 气泡图（分面：Up / Down） ----
  p <- ggplot(plot_data, aes(x = Sample, y = Description_short)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "red", high = "blue", trans = "reverse") +
    scale_size_continuous(range = c(2, 8)) +
    facet_wrap(~ Direction, scales = "free_y", ncol = 2) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 8),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "right"
    ) +
    labs(
      title = "GOBP Top5: Inner vs Outer",
      x = NULL, y = NULL,
      size = "Gene Count", color = "p.adjust"
    )
  
  ggsave(file.path(go_dir, "GOBP_Top5_UpDown.pdf"), p, width = 14, height = 10, dpi = 300)
  cat("   ✅ 已保存:", file.path(go_dir, "GOBP_Top5_UpDown.pdf"), "\n")
}

# 保存完整结果
saveRDS(list(up = all_go_up, down = all_go_down), 
        file = file.path(go_dir, "GOBP_UpDown_Results.rds"))
cat("\n✅ GOBP 双向分析完成！\n")



# ============================================================
#  GOBP 可视化（对标单细胞 dotplot 风格）
#  输出：每个样本的 dotplot + 汇总气泡图
# ============================================================
# ============================================================
#  修正版 GOBP 分析：logFC > 1，up/down 分别处理
#  输出：每个样本 4 个 Excel 文件 + 汇总可视化
# ============================================================

rm(list = ls())
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(writexl)  # 用于导出 Excel

# ----- 1. 路径设置 -----
osc_dir <- "E:/NETS/空转组数据/结果/Lesion_Analysis/Oscillation_Field"
go_dir <- "E:/NETS/空转组数据/结果/Lesion_Analysis/GOBP_InnerOuter"
if (!dir.exists(go_dir)) dir.create(go_dir, recursive = TRUE)

seurat_list <- readRDS(file.path(osc_dir, "seurat_list_with_weights.rds"))
sample_names <- names(seurat_list)

# ----- 2. 分层函数（保持原样）-----
get_location_fixed <- function(obj) {
  grade_vec <- obj$oscillation_grade
  lesion_vec <- obj$lesion_id
  
  location <- rep(NA, length(grade_vec))
  
  all_levels <- sort(unique(grade_vec[!is.na(grade_vec) & grade_vec != "Zero"]))
  if (length(all_levels) == 0) return(location)
  
  sorted_levels <- sort(all_levels, decreasing = TRUE)
  n <- length(sorted_levels)
  inner_levels <- sorted_levels[1:floor(n/2)]
  outer_levels <- sorted_levels[(floor(n/2)+1):n]
  
  for (i in seq_along(grade_vec)) {
    g <- grade_vec[i]
    if (is.na(g) || g == "Zero") {
      location[i] <- NA
    } else if (g %in% inner_levels) {
      location[i] <- "Inner"
    } else if (g %in% outer_levels) {
      location[i] <- "Outer"
    }
  }
  
  is_center <- !is.na(lesion_vec)
  location[is_center] <- "Inner"
  
  return(location)
}

# ----- 3. 批量分析（up + down 分别处理）-----
MIN_GENES_FOR_GO <- 10  # 最少基因数阈值
all_results <- list()

for (s in sample_names) {
  cat("\n========================================\n")
  cat("🔄 处理样本:", s, "\n")
  obj <- seurat_list[[s]]
  
  if (!"oscillation_grade" %in% colnames(obj@meta.data)) {
    cat("   ⚠️ 缺少 oscillation_grade 列，跳过\n")
    next
  }
  
  obj$location <- get_location_fixed(obj)
  
  cat("   Inner spots:", sum(obj$location == "Inner", na.rm = TRUE), "\n")
  cat("   Outer spots:", sum(obj$location == "Outer", na.rm = TRUE), "\n")
  cat("   NA (剔除):", sum(is.na(obj$location)), "\n")
  
  if (sum(obj$location %in% c("Inner", "Outer"), na.rm = TRUE) < 10) {
    cat("   ⚠️ 有效spots太少，跳过\n")
    next
  }
  
  obj_sub <- subset(obj, subset = !is.na(location))
  Idents(obj_sub) <- "location"
  
  # ---- 差异表达分析（logfc.threshold = 1） ----
  cat("   ⏳ 运行 FindMarkers (logFC > 1)...\n")
  markers <- FindMarkers(
    obj_sub,
    ident.1 = "Inner",
    ident.2 = "Outer",
    assay = "Spatial",
    slot = "data",
    logfc.threshold = 1,        # 关键：改为 1
    test.use = "wilcox",
    min.pct = 0.1
  )
  
  # ---- 提取 up 和 down 基因 ----
  sig_markers <- markers[markers$p_val_adj < 0.05, ]
  genes_up <- rownames(sig_markers[sig_markers$avg_log2FC > 0, ])
  genes_down <- rownames(sig_markers[sig_markers$avg_log2FC < 0, ])
  
  cat("   ✅ 上调基因（Inner高表达）:", length(genes_up), "\n")
  cat("   ✅ 下调基因（Outer高表达）:", length(genes_down), "\n")
  
  # ---- 保存差异基因 Excel ----
  deg_list <- list(
    Up_Genes = data.frame(Gene = genes_up, sig_markers[sig_markers$avg_log2FC > 0, ]),
    Down_Genes = data.frame(Gene = genes_down, sig_markers[sig_markers$avg_log2FC < 0, ])
  )
  write_xlsx(deg_list, path = file.path(go_dir, paste0("DEGs_", s, ".xlsx")))
  cat("   💾 已保存差异基因表格\n")
  
  # ---- 处理 up 基因的 GOBP ----
  go_up <- NULL
  if (length(genes_up) >= MIN_GENES_FOR_GO) {
    cat("   ⏳ 运行 GOBP (Up genes)...\n")
    go_up <- tryCatch({
      enrichGO(
        gene = genes_up,
        OrgDb = org.Hs.eg.db,
        keyType = "SYMBOL",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = TRUE
      )
    }, error = function(e) NULL)
  }
  
  # ---- 处理 down 基因的 GOBP ----
  go_down <- NULL
  if (length(genes_down) >= MIN_GENES_FOR_GO) {
    cat("   ⏳ 运行 GOBP (Down genes)...\n")
    go_down <- tryCatch({
      enrichGO(
        gene = genes_down,
        OrgDb = org.Hs.eg.db,
        keyType = "SYMBOL",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = TRUE
      )
    }, error = function(e) NULL)
  }
  
  # ---- 保存 GOBP 结果（Excel） ----
  go_list <- list()
  if (!is.null(go_up) && nrow(go_up@result) > 0) {
    go_list$GOBP_Up <- as.data.frame(go_up@result)
    cat("   ✅ Up GOBP:", nrow(go_list$GOBP_Up), "个通路\n")
  }
  if (!is.null(go_down) && nrow(go_down@result) > 0) {
    go_list$GOBP_Down <- as.data.frame(go_down@result)
    cat("   ✅ Down GOBP:", nrow(go_list$GOBP_Down), "个通路\n")
  }
  
  if (length(go_list) > 0) {
    write_xlsx(go_list, path = file.path(go_dir, paste0("GOBP_", s, ".xlsx")))
    cat("   💾 已保存 GOBP 结果\n")
  } else {
    cat("   ⚠️ 无显著 GOBP 通路\n")
  }
  
  # 存储结果用于后续可视化
  all_results[[s]] <- list(
    genes_up = genes_up,
    genes_down = genes_down,
    go_up = go_up,
    go_down = go_down,
    n_up = length(genes_up),
    n_down = length(genes_down)
  )
}

# ----- 4. 汇总统计（展示哪些样本通过阈值）-----
cat("\n========================================\n")
cat("📊 差异基因汇总:\n")
summary_df <- data.frame(
  Sample = names(all_results),
  Up_DEGs = sapply(all_results, function(x) x$n_up),
  Down_DEGs = sapply(all_results, function(x) x$n_down),
  GO_Up = sapply(all_results, function(x) ifelse(!is.null(x$go_up), nrow(x$go_up@result), 0)),
  GO_Down = sapply(all_results, function(x) ifelse(!is.null(x$go_down), nrow(x$go_down@result), 0))
)
print(summary_df)
write_xlsx(summary_df, path = file.path(go_dir, "GOBP_Summary.xlsx"))

# ----- 5. 可视化：每个样本的 dotplot（up 和 down 分别）-----
vis_dir <- file.path(go_dir, "Visualization")
if (!dir.exists(vis_dir)) dir.create(vis_dir, recursive = TRUE)

for (s in names(all_results)) {
  res <- all_results[[s]]
  
  # ---- Up genes GOBP dotplot ----
  if (!is.null(res$go_up) && nrow(res$go_up@result) > 0) {
    df <- as.data.frame(res$go_up@result)
    df <- df[order(df$p.adjust), ]
    n_show <- min(15, nrow(df))
    plot_df <- df[1:n_show, ]
    plot_df$Description <- factor(plot_df$Description, levels = rev(plot_df$Description))
    
    p_up <- ggplot(plot_df, aes(x = Count, y = Description, size = Count, color = p.adjust)) +
      geom_point(alpha = 0.9) +
      scale_color_gradient(low = "#E41A1C", high = "#377EB8", trans = "reverse", name = "p.adjust") +
      scale_size_continuous(range = c(3, 10), name = "Gene Count") +
      theme_bw(base_size = 11) +
      theme(
        axis.text.y = element_text(size = 9, face = "bold"),
        axis.text.x = element_text(size = 9),
        axis.title = element_text(size = 10),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
      ) +
      labs(title = paste0(s, " - GO BP (Up in Inner)"), x = "Gene Count", y = NULL)
    
    ggsave(file.path(vis_dir, paste0("GOBP_dotplot_Up_", s, ".pdf")), 
           p_up, width = 10, height = 5 + n_show * 0.2, dpi = 300, device = "pdf")
  }
  
  # ---- Down genes GOBP dotplot ----
  if (!is.null(res$go_down) && nrow(res$go_down@result) > 0) {
    df <- as.data.frame(res$go_down@result)
    df <- df[order(df$p.adjust), ]
    n_show <- min(15, nrow(df))
    plot_df <- df[1:n_show, ]
    plot_df$Description <- factor(plot_df$Description, levels = rev(plot_df$Description))
    
    p_down <- ggplot(plot_df, aes(x = Count, y = Description, size = Count, color = p.adjust)) +
      geom_point(alpha = 0.9) +
      scale_color_gradient(low = "#2C3E50", high = "#95A5A6", trans = "reverse", name = "p.adjust") +
      scale_size_continuous(range = c(3, 10), name = "Gene Count") +
      theme_bw(base_size = 11) +
      theme(
        axis.text.y = element_text(size = 9, face = "bold"),
        axis.text.x = element_text(size = 9),
        axis.title = element_text(size = 10),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
      ) +
      labs(title = paste0(s, " - GO BP (Up in Outer)"), x = "Gene Count", y = NULL)
    
    ggsave(file.path(vis_dir, paste0("GOBP_dotplot_Down_", s, ".pdf")), 
           p_down, width = 10, height = 5 + n_show * 0.2, dpi = 300, device = "pdf")
  }
}

cat("\n========================================\n")
cat("✅ 所有 GOBP 分析完成！\n")
cat("📁 输出目录:", go_dir, "\n")
cat("   - DEGs_*.xlsx (up/down 基因表格)\n")
cat("   - GOBP_*.xlsx (up/down 通路表格)\n")
cat("   - GOBP_Summary.xlsx (汇总统计)\n")
cat("   - Visualization/GOBP_dotplot_*.pdf (dotplot)\n")
cat("========================================\n")


