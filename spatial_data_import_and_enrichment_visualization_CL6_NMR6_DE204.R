#=======================
# 解包 + 整理（完整版）
#=======================

rm(list = ls())

base_dir <- "E:/NETS/空转组数据/原数据"

# 注意：这里的名称必须和你的压缩文件名完全一致！
# 如果你的压缩包叫 GBM_ZH1007_inf.tar.gz，就用 "GBM_ZH1007_inf"
sample_names <- c(
  "GBM_ZH1007inf",
  "GBM_ZH1007nec",
  "GBM_ZH916inf",
  "GBM_ZH881T1",
  "GBM_ZH1019T1",
  "GBM_ZH916T1",
  "GBM_ZH881inf"
)

setwd(base_dir)

# 创建报告
report <- data.frame(sample = character(), status = character(), stringsAsFactors = FALSE)

for (s in sample_names) {
  
  cat("\n========================================\n")
  cat("正在处理样本:", s, "\n")
  
  tar_file <- file.path(base_dir, paste0(s, ".tar.gz"))  # 压缩包路径
  sample_dir <- file.path(base_dir, s)                   # 解压目标文件夹
  
  # ============================================================
  # 第一步：解压 .tar.gz 文件
  # ============================================================
  if (file.exists(tar_file)) {
    cat("   📦 找到压缩包:", basename(tar_file), "\n")
    
    # 如果目标文件夹已存在，先删除（避免旧文件干扰）
    if (dir.exists(sample_dir)) {
      cat("   ⚠️ 目标文件夹已存在，正在删除旧版本...\n")
      unlink(sample_dir, recursive = TRUE)
    }
    
    # 执行解压
    cat("   ⏳ 正在解压，请稍候...\n")
    untar_result <- untar(tar_file, exdir = base_dir, tar = "internal")
    
    if (untar_result == 0) {
      cat("   ✅ 解压成功！\n")
    } else {
      cat("   ❌ 解压失败，错误代码:", untar_result, "\n")
      report <- rbind(report, data.frame(sample = s, status = "解压失败"))
      next
    }
    
  } else {
    cat("   ❌ 未找到压缩包:", tar_file, "\n")
    cat("   💡 请检查文件名是否与压缩包完全一致（包括大小写和下划线）\n")
    report <- rbind(report, data.frame(sample = s, status = "压缩包缺失"))
    next
  }
  
  # ============================================================
  # 第二步：检查解压后的文件夹
  # ============================================================
  if (!dir.exists(sample_dir)) {
    cat("   ❌ 解压后未找到文件夹:", sample_dir, "\n")
    report <- rbind(report, data.frame(sample = s, status = "解压后文件夹缺失"))
    next
  }
  
  # 列出解压后的内容
  file_list <- list.files(sample_dir)
  cat("   📂 解压后包含", length(file_list), "个文件/文件夹\n")
  if (length(file_list) > 0) {
    cat("      前5个:", paste(head(file_list, 5), collapse = ", "), "\n")
  }
  
  # ============================================================
  # 第三步：处理可能的双层嵌套（常见于Zenodo数据）
  # ============================================================
  inner_dir <- file.path(sample_dir, s)
  if (dir.exists(inner_dir)) {
    cat("   📂 发现双层嵌套，正在提取内层文件...\n")
    inner_files <- list.files(inner_dir, full.names = TRUE)
    for (f in inner_files) {
      file.rename(f, file.path(sample_dir, basename(f)))
    }
    unlink(inner_dir, recursive = TRUE)
    cat("   ✅ 已提取", length(inner_files), "个文件到外层\n")
  } else {
    cat("   ✅ 无双层嵌套\n")
  }
  
  # ============================================================
  # 第四步：删除 __MACOSX
  # ============================================================
  mac_dir <- file.path(sample_dir, "__MACOSX")
  if (dir.exists(mac_dir)) {
    unlink(mac_dir, recursive = TRUE)
    cat("   🗑️ 已删除 __MACOSX\n")
  }
  
  # ============================================================
  # 第五步：提取前缀
  # ============================================================
  prefix <- tolower(gsub("GBM_", "", s))
  
  # ============================================================
  # 第六步：创建 spatial 文件夹并移动文件
  # ============================================================
  spatial_dir <- file.path(sample_dir, "spatial")
  dir.create(spatial_dir, recursive = TRUE, showWarnings = FALSE)
  
  src_names <- c(
    paste0(prefix, "_tissue_positions_list.csv"),
    paste0(prefix, "_scalefactors_json.json"),
    paste0(prefix, "_tissue_lowres_image.png"),
    paste0(prefix, "_aligned_fiducials.jpg"),
    paste0(prefix, "_detected_tissue_image.jpg")
  )
  dst_names <- c(
    "tissue_positions_list.csv",
    "scalefactors_json.json",
    "tissue_lowres_image.png",
    "aligned_fiducials.jpg",
    "detected_tissue_image.jpg"
  )
  
  moved_count <- 0
  for (i in seq_along(src_names)) {
    src_path <- file.path(sample_dir, src_names[i])
    dst_path <- file.path(spatial_dir, dst_names[i])
    if (file.exists(src_path)) {
      if (file.exists(dst_path)) file.remove(dst_path)
      file.rename(src_path, dst_path)
      moved_count <- moved_count + 1
    }
  }
  cat("   ✅ 已移动", moved_count, "个文件到 spatial/ 文件夹\n")
  
  # ============================================================
  # 第七步：重命名主 H5 文件
  # ============================================================
  h5_src <- file.path(sample_dir, paste0(prefix, "_filtered_feature_bc_matrix.h5"))
  h5_dst <- file.path(sample_dir, "filtered_feature_bc_matrix.h5")
  if (file.exists(h5_src)) {
    if (file.exists(h5_dst)) file.remove(h5_dst)
    file.rename(h5_src, h5_dst)
    cat("   ✅ 已重命名主H5文件\n")
  } else {
    # 有些样本可能是老版本，用文件夹而非H5
    h5_folder <- file.path(sample_dir, paste0(prefix, "_filtered_feature_bc_matrix"))
    if (dir.exists(h5_folder)) {
      # 如果是文件夹格式，重命名文件夹
      h5_dst_folder <- file.path(sample_dir, "filtered_feature_bc_matrix")
      if (dir.exists(h5_dst_folder)) unlink(h5_dst_folder, recursive = TRUE)
      file.rename(h5_folder, h5_dst_folder)
      cat("   ✅ 已重命名主H5文件夹（非H5格式）\n")
    } else {
      cat("   ⚠️ 未找到H5文件或H5文件夹\n")
    }
  }
  
  # ============================================================
  # 第八步：最终校验
  # ============================================================
  h5_check <- file.exists(h5_dst) || dir.exists(file.path(sample_dir, "filtered_feature_bc_matrix"))
  spatial_check <- dir.exists(spatial_dir) && file.exists(file.path(spatial_dir, "tissue_positions_list.csv"))
  
  if (h5_check && spatial_check) {
    cat("   ✅ 样本", s, "整理完成！\n")
    report <- rbind(report, data.frame(sample = s, status = "成功"))
  } else {
    cat("   ❌ 样本", s, "整理失败，缺少核心文件\n")
    report <- rbind(report, data.frame(sample = s, status = "失败"))
  }
}

cat("\n========================================\n")
cat("📊 处理报告：\n")
print(report)
cat("========================================\n")



# ============================================================
#  GSE237183 空间转录组：完整预处理流程（标准化+降维+聚类）
#  优化要点：
#    1. 自动匹配所有已解压的样本文件夹（不需要手动写sample_names）
#    2. 添加内存清理（gc()），防止7个样本累积占满内存
#    3. 保存图片时统一使用 result_dir
#    4. 增加低质量spoT的详细过滤报告
# ============================================================

rm(list = ls())
library(Seurat)
library(ggplot2)
library(patchwork)

# ----- 1. 设置路径 -----
base_dir <- "E:/NETS/空转组数据/原数据"
result_dir <- "E:/NETS/空转组数据/结果"
processed_dir <- file.path(base_dir, "processed_rds")

# 创建结果文件夹（如果不存在）
if (!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)
if (!dir.exists(processed_dir)) dir.create(processed_dir, recursive = TRUE)

# ----- 2. 自动识别所有已解压的样本文件夹（避免手动遗漏） -----
all_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)
# 只保留符合 GBM_* 格式的文件夹（排除 processed_rds 等）
sample_names <- all_dirs[grepl("^GBM_", all_dirs)]

cat("🔍 自动发现以下", length(sample_names), "个样本文件夹：\n")
print(sample_names)
cat("\n")

# 记录每个样本的处理状态
report <- data.frame(sample = character(), 
                     spots_raw = integer(), 
                     spots_kept = integer(),
                     clusters = integer(),
                     status = character(),
                     stringsAsFactors = FALSE)

# ----- 3. 循环处理每个样本 -----
for (s in sample_names) {
  
  cat("\n========================================\n")
  cat("🔄 正在处理样本:", s, "\n")
  cat("========================================\n")
  
  sample_dir <- file.path(base_dir, s)
  
  # ---------- 3.1 读取数据（增加容错） ----------
  obj <- tryCatch({
    Load10X_Spatial(
      data.dir = sample_dir,
      filename = "filtered_feature_bc_matrix.h5"
    )
  }, error = function(e) {
    cat("   ❌ 读取失败:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(obj)) {
    report <- rbind(report, data.frame(
      sample = s, spots_raw = NA, spots_kept = NA, 
      clusters = NA, status = "读取失败"
    ))
    next
  }
  
  # 添加样本标签和线粒体百分比
  obj$sample_origin <- s
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  raw_spots <- ncol(obj)
  cat("📊 原始数据：", raw_spots, "个spots,", nrow(obj), "个基因\n")
  
  # ---------- 3.2 质控图（保存到结果文件夹） ----------
  cat("📊 生成质控图...\n")
  p_qc <- VlnPlot(obj, 
                  features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"),
                  ncol = 3, pt.size = 0.1) +
    patchwork::plot_annotation(title = paste0("QC - ", s))
  
  ggsave(
    filename = file.path(result_dir, paste0("QC_", s, ".png")),
    plot = p_qc,
    width = 12, height = 6, dpi = 300
  )
  
  # ---------- 3.3 过滤低质量spots（针对坏死区适当放宽） ----------
  # 对于nec样本，线粒体比例可能天然偏高，阈值适当放宽到25%
  mito_threshold <- ifelse(grepl("nec", s), 25, 20)
  
  obj <- subset(obj, 
                subset = nFeature_Spatial > 200 & 
                  nCount_Spatial > 500 & 
                  percent.mt < mito_threshold)
  
  kept_spots <- ncol(obj)
  cat("🧹 过滤后剩余：", kept_spots, "个spots (线粒体阈值:", mito_threshold, "%)\n")
  
  # 如果过滤后spots太少（<100），跳过后续分析
  if (kept_spots < 100) {
    cat("   ⚠️ 剩余spots过少，跳过该样本的聚类分析\n")
    report <- rbind(report, data.frame(
      sample = s, spots_raw = raw_spots, spots_kept = kept_spots,
      clusters = NA, status = "spots过少"
    ))
    next
  }
  
  # ---------- 3.4 标准化、降维、聚类 ----------
  cat("⚙️ 正在进行标准化和降维...\n")
  
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, nfeatures = 2000)
  obj <- ScaleData(obj, vars.to.regress = "percent.mt")
  obj <- RunPCA(obj, npcs = 30)
  
  # 判断PCA维度（如果spots很少，用较少维度）
  n_dims <- min(30, floor(kept_spots / 5))
  cat("   使用", n_dims, "个PC维度\n")
  
  obj <- FindNeighbors(obj, dims = 1:n_dims)
  obj <- FindClusters(obj, resolution = 0.5)
  
  # UMAP降维（用于可视化参考）
  obj <- RunUMAP(obj, dims = 1:n_dims)
  
  n_clusters <- length(unique(obj$seurat_clusters))
  cat("✅ 聚类完成，共", n_clusters, "个群\n")
  
  # ---------- 3.5 保存空间聚类图 ----------
  cat("🖼️ 绘制空间聚类图（HE透明度0.4，保存为PDF）...\n")
  
  p_cluster <- SpatialDimPlot(
    obj,
    label = TRUE,
    label.size = 3,
    repel = TRUE,
    image.alpha = 0.4          # HE图透明度，0~1，越小越淡
  ) +
    patchwork::plot_annotation(
      title = paste0(s, " - Spatial Clusters (", n_clusters, " clusters)")
    )
  
  ggsave(
    filename = file.path(result_dir, paste0("SpatialClusters_", s, ".pdf")),
    plot = p_cluster,
    width = 10,                # 宽度10英寸
    height = 8,                # 高度8英寸
    device = "pdf"             # 输出PDF格式
  )
  
  cat("✅ 已保存PDF:", file.path(result_dir, paste0("SpatialClusters_", s, ".pdf")), "\n")
  
  # ---------- 3.6 保存处理后的R对象 ----------
  saveRDS(obj, file = file.path(processed_dir, paste0("processed_", s, ".rds")))
  cat("💾 处理后的数据已保存至:", file.path(processed_dir, paste0("processed_", s, ".rds")), "\n")
  
  # 记录报告
  report <- rbind(report, data.frame(
    sample = s, spots_raw = raw_spots, spots_kept = kept_spots,
    clusters = n_clusters, status = "成功"
  ))
  
  # ---------- 3.7 清理内存 ----------
  rm(obj)
  gc()
  
  cat("✅ 样本", s, "全部处理完成！\n")
}

# ----- 4. 输出总结报告 -----
cat("\n========================================\n")
cat("📊 处理报告汇总\n")
cat("========================================\n")
print(report)

# 统计成功/失败
success_count <- sum(report$status == "成功")
total_count <- nrow(report)
cat("\n✅ 成功:", success_count, "/", total_count, "个样本\n")

if (success_count == total_count) {
  cat("🎉 所有样本预处理完毕！\n")
} else {
  cat("⚠️ 有", total_count - success_count, "个样本处理失败，请检查日志\n")
}

cat("\n📁 结果文件位置：\n")
cat("  - 质控图: ", result_dir, "\n")
cat("  - 处理后的RDS: ", processed_dir, "\n")
cat("========================================\n")


# ============================================================
#  空间NETs评分可视化（两种模式）
#  1. 圆点 + 半透明HE背景
#  2. 纯六边形拼接图（无背景）
#  输出：每个样本两个PDF文件
# ============================================================

rm(list = ls())
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)      # 用于数据操作
library(tidyr)      # 用于长格式转换
library(FNN)  

# ----- 1. 路径设置 -----
base_dir <- "E:/NETS/空转组数据/原数据"
result_dir <- "E:/NETS/空转组数据/结果"
sc_results_dir <- "E:/NETS/单细胞数据/结果"
processed_rds_dir <- file.path(base_dir, "processed_rds")   # 存放预处理后的rds

# 确保结果目录存在
if (!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)

# ----- 2. 加载已处理的空间对象（或从原始文件夹重新读取，这里直接加载rds）-----
# 假设你已经用预处理脚本生成了 processed_*.rds 文件并放在 processed_rds_dir 中
sample_names <- c(
  "GBM_ZH1007inf",
  "GBM_ZH1007nec",
  "GBM_ZH916inf",
  "GBM_ZH881T1",
  "GBM_ZH1019T1",
  "GBM_ZH916T1",
  "GBM_ZH881inf"
)

cat("📂 加载预处理空间数据...\n")
seurat_list <- list()
for (s in sample_names) {
  rds_file <- file.path(processed_rds_dir, paste0("processed_", s, ".rds"))
  if (file.exists(rds_file)) {
    seurat_list[[s]] <- readRDS(rds_file)
    cat("   ✅", s, "加载成功 (", ncol(seurat_list[[s]]), "spots )\n")
  } else {
    # 如果rds不存在，尝试从原始文件夹直接读取
    raw_dir <- file.path(base_dir, s)
    if (dir.exists(raw_dir)) {
      cat("   ⚠️ 未找到rds文件，直接从原始文件夹读取:", s, "\n")
      obj <- Load10X_Spatial(data.dir = raw_dir, filename = "filtered_feature_bc_matrix.h5")
      obj$sample_origin <- s
      obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
      # 简单过滤（与原预处理保持一致）
      mito_threshold <- ifelse(grepl("nec", s), 25, 20)
      obj <- subset(obj, subset = nFeature_Spatial > 200 & nCount_Spatial > 500 & percent.mt < mito_threshold)
      seurat_list[[s]] <- obj
      cat("   ✅", s, "读取并过滤完成 (", ncol(obj), "spots )\n")
    } else {
      cat("   ❌ 样本", s, "未找到，跳过\n")
    }
  }
}

# 移除未加载的样本
seurat_list <- seurat_list[!sapply(seurat_list, is.null)]
if (length(seurat_list) == 0) stop("没有可用的空间样本，请检查路径！")

# ----- 3. 定义基因集（经典6、204、NMR）-----
cat("\n🧬 加载基因集...\n")

# 3.1 经典6基因
genes_classic <- c("PTGS2", "MME", "SLC2A3", "PADI4", "ELANE", "MPO")
genes_classic <- genes_classic[genes_classic %in% rownames(seurat_list[[1]])]
cat("   ✅ 经典6基因集:", length(genes_classic), "个\n")

# 3.2 204基因集（从单细胞DE结果读取）
de_file <- file.path(sc_results_dir, "DE_up_genes_filtered.csv")
genes_204 <- NULL
if (file.exists(de_file)) {
  df <- read.csv(de_file)
  if ("gene" %in% colnames(df)) genes_204 <- as.character(df$gene) else genes_204 <- as.character(df[,1])
  genes_204 <- genes_204[genes_204 %in% rownames(seurat_list[[1]])]
  cat("   ✅ 204基因集:", length(genes_204), "个\n")
} else {
  cat("   ⚠️ 未找到204基因集文件，跳过\n")
}

# 3.3 NMR 6基因集（弹性网络结果）
enet_file <- "E:/NETS/elasticnet_cox_results.RData"
genes_nmr <- NULL
if (file.exists(enet_file)) {
  load(enet_file)  # 假设有 final_genes
  if (exists("final_genes")) {
    genes_nmr <- final_genes
    genes_nmr <- genes_nmr[genes_nmr %in% rownames(seurat_list[[1]])]
    cat("   ✅ NMR基因集:", length(genes_nmr), "个\n")
  }
} else {
  cat("   ⚠️ 未找到弹性网络结果文件，跳过\n")
}

# 组合所有基因集（只保留非空）
gene_sets <- list(
  "Classic_6" = genes_classic,
  "DE_204"    = genes_204,
  "NMR_6"     = genes_nmr
)
gene_sets <- gene_sets[!sapply(gene_sets, is.null) & sapply(gene_sets, length) > 0]
cat("\n📊 共使用", length(gene_sets), "个基因集\n")

# ----- 4. 计算或加载NETs评分 -----
# 检查是否已有保存的包含评分的列表
score_file <- file.path(base_dir, "seurat_list_with_NETs_scores.rds")
if (file.exists(score_file)) {
  cat("\n📂 发现已保存的评分结果，直接加载...\n")
  seurat_list <- readRDS(score_file)
} else {
  cat("\n⚙️ 正在计算NETs评分（使用UCell）...\n")
  if (!require("UCell", quietly = TRUE)) install.packages("UCell")
  library(UCell)
  
  for (s in names(seurat_list)) {
    obj <- seurat_list[[s]]
    for (set_name in names(gene_sets)) {
      genes <- gene_sets[[set_name]]
      if (length(genes) > 1) {
        obj <- AddModuleScore_UCell(obj, features = list(genes), 
                                    name = paste0("NETs_", set_name, "_"))
      }
    }
    seurat_list[[s]] <- obj
    cat("   ✅", s, "打分完成\n")
  }
  # 保存评分结果，供后续使用
  saveRDS(seurat_list, file = score_file)
  cat("💾 评分结果已保存至:", score_file, "\n")
}


# ----- 5. 绘图函数 -----

# 5.1 圆点叠加 HE（保持不变）
plot_overlay <- function(obj, features, image.alpha = 0.15, title_suffix = "") {
  p <- SpatialFeaturePlot(obj, 
                          features = features, 
                          ncol = length(features),
                          pt.size.factor = 1.6,
                          alpha = c(0.1, 1),
                          image.alpha = image.alpha) +
    patchwork::plot_annotation(title = paste0(unique(obj$sample_origin), title_suffix))
  return(p)
}

# 5.2 平顶六边形拼接（离散三色）
plot_hex <- function(obj, features, title_suffix = "") {
  # ---- 坐标 ----
  coords <- GetTissueCoordinates(obj, scale = "lowres")
  if (all(c("imagerow", "imagecol") %in% colnames(coords))) {
    xcol <- "imagerow"; ycol <- "imagecol"
  } else if (all(c("row", "col") %in% colnames(coords))) {
    xcol <- "row"; ycol <- "col"
  } else {
    xcol <- colnames(coords)[1]; ycol <- colnames(coords)[2]
  }
  
  # ---- 评分 ----
  meta <- obj@meta.data
  score_cols <- intersect(features, colnames(meta))
  if (length(score_cols) == 0) stop("无匹配评分列！")
  
  df <- coords
  for (col in score_cols) df[[col]] <- meta[[col]]
  df_long <- df %>% pivot_longer(cols = all_of(score_cols), names_to = "Score", values_to = "Value")
  
  # ---- 六边形边长（间隙系数 0.80） ----
  if (!require("FNN", quietly = TRUE)) install.packages("FNN")
  library(FNN)
  coords_mat <- as.matrix(coords[, c(xcol, ycol)])
  nn_dist <- FNN::knn.dist(coords_mat, k = 2)[,2]
  side <- median(nn_dist, na.rm = TRUE) / sqrt(3) * 0.80   # 明显空隙
  if (is.na(side) || side < 0.1) side <- 5
  
  # ---- 平顶六边形顶点 ----
  angles <- ( 60 * (0:5)) * pi / 180
  
  hex_list <- list()
  for (i in 1:nrow(df_long)) {
    xc <- df_long[[xcol]][i]
    yc <- df_long[[ycol]][i]
    vertices <- data.frame(
      x = xc + side * cos(angles),
      y = yc + side * sin(angles)
    )
    vertices$group <- i
    vertices$Score <- df_long$Score[i]
    vertices$Value <- df_long$Value[i]
    hex_list[[i]] <- vertices
  }
  hex_df <- do.call(rbind, hex_list)
  
  # ---- 每个评分独立连续渐变（白-黄-橙），中点 = (min+max)/2 ----
  unique_scores <- unique(hex_df$Score)
  plot_list <- list()
  for (sc in unique_scores) {
    sub_df <- hex_df[hex_df$Score == sc, ]
    v_min <- min(sub_df$Value, na.rm = TRUE)
    v_max <- max(sub_df$Value, na.rm = TRUE)
    v_mid <- (v_min + v_max) / 2   # 核心改动：值域中点
    
    p <- ggplot(sub_df, aes(x = x, y = y, fill = Value)) +
      geom_polygon(aes(group = group), color = NA, size = 0) +
      coord_fixed() +
      scale_fill_gradient2(
        low = "#FFFFFF",      # 白色
        mid = "#FDBF2F",      # 金黄色
        high = "#E31A1C",     # 橙红色
        midpoint = v_mid,
        limits = c(v_min, v_max)
      ) +
      theme_void() +
      theme(legend.position = "bottom",
            legend.key.width = unit(1.5, "cm"),
            plot.title = element_text(hjust = 0.5)) +
      labs(title = sc)
    plot_list[[sc]] <- p
  }
  
  combined <- wrap_plots(plot_list, ncol = length(plot_list))
  combined <- combined + plot_annotation(title = paste0(unique(obj$sample_origin), title_suffix))
  return(combined)
}

# ----- 6. 批量生成 -----
cat("\n📄 生成离散三色 PDF...\n")
for (s in names(seurat_list)) {
  obj <- seurat_list[[s]]
  score_cols <- grep("NETs", colnames(obj@meta.data), value = TRUE, ignore.case = TRUE)
  if (length(score_cols) == 0) next
  
  cat("\n🔄 处理样本:", s, "\n")
  cat("   评分列:", paste(score_cols, collapse = ", "), "\n")
  
  p1 <- plot_overlay(obj, score_cols, image.alpha = 0.15, title_suffix = " (Overlay)")
  ggsave(file.path(result_dir, paste0("NETs_Overlay_", s, ".pdf")),
         p1, width = 5 * length(score_cols), height = 6, dpi = 300, device = "pdf")
  
  p2 <- plot_hex(obj, score_cols, title_suffix = " (Hex)")
  ggsave(file.path(result_dir, paste0("NETs_Hex_", s, ".pdf")),
         p2, width = 5 * length(score_cols), height = 6, dpi = 300, device = "pdf")
  
  cat("   ✅ 已保存 Overlay 和 Hex PDF\n")
}
cat("\n🎉 全部完成！\n")

