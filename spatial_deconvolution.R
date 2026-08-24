rm(list = ls())
library(Seurat)
library(readxl)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(pheatmap)

# ----- 1. 路径设置 -----
sc_data_dir <- "E:/NETS/单细胞数据/结果"
sc_obj_path <- file.path(sc_data_dir, "cgga_scRNA_processed_final.rds")
mp_file <- "E:/NETS/空转组数据/原数据/mmc2.xlsx"
out_dir <- "E:/NETS/空转组数据/结果/Reference_Annotation"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ----- 2. 加载数据 -----
cat("📂 加载原始单细胞Seurat对象...\n")
sc_obj <- readRDS(sc_obj_path)
cat("   ✅ 细胞数:", ncol(sc_obj), "\n")
cat("   ✅ 当前聚类数:", length(unique(sc_obj$seurat_clusters)), "\n")

# ----- 3. 提高分辨率重新聚类（resolution = 0.9）-----
cat("⚙️ 提高分辨率重新聚类（resolution = 0.9）...\n")
sc_obj <- FindClusters(sc_obj, resolution = 0.9)
new_clusters <- length(unique(sc_obj$seurat_clusters))
cat("   ✅ 新聚类数:", new_clusters, "\n")

# ----- 4. 读取GBM MPs（14个）-----
cat("\n📂 读取GBM MPs基因列表...\n")
mp_df <- read_excel(mp_file, sheet = "GBM_MPs")
mp_list <- list()
for (col_name in colnames(mp_df)) {
  genes <- mp_df[[col_name]]
  genes <- genes[!is.na(genes) & genes != ""]
  clean_name <- gsub(" \\(malig\\.\\)", "", col_name)
  clean_name <- gsub(" ", "_", clean_name)
  clean_name <- gsub("\\.", "_", clean_name)
  mp_list[[clean_name]] <- genes
}
cat("   ✅ 加载了", length(mp_list), "个GBM MPs\n")

# ----- 5. 构建自定义中性粒细胞标志集-----
cat("\n🧬 构建自定义中性粒细胞标志集...\n")
neutrophil_genes <- c(
  "FCGR3B", "CSF3R", "CXCR2", "PROK2",
  "ELANE", "MPO", "CTSG", "AZU1", "PRTN3", "MMP9",
  "S100A8", "S100A9", "S100A12",
  "CXCL8", "IL1R2", "CXCR1"
)
neutrophil_genes <- neutrophil_genes[neutrophil_genes %in% rownames(sc_obj)]
cat("   ✅ 中性粒细胞标志集:", length(neutrophil_genes), "个基因\n")
ref_list <- mp_list
ref_list[["Neutrophil_Custom"]] <- neutrophil_genes

# ----- 6. 计算每个新cluster的Top50差异基因-----
cat("\n⚙️ 计算每个cluster的Top50差异基因...\n")
Idents(sc_obj) <- "seurat_clusters"
all_markers <- FindAllMarkers(
  sc_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)
top50_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 50) %>%
  ungroup()
cluster_genes_list <- split(top50_markers$gene, top50_markers$cluster)
cat("   ✅ 共计算", length(cluster_genes_list), "个cluster的Top50\n")

# ----- 7. 富集分析：超几何检验 -----
cat("\n🔬 执行富集分析...\n")
background_genes <- rownames(sc_obj)
N <- length(background_genes)
enrichment_results <- list()

for (cl in names(cluster_genes_list)) {
  cluster_genes <- cluster_genes_list[[cl]]
  k <- length(cluster_genes)
  cl_res <- data.frame(
    Cluster = cl,
    Reference = names(ref_list),
    Overlap = NA,
    P_value = NA,
    P_adjust = NA
  )
  for (ref_name in names(ref_list)) {
    ref_genes <- ref_list[[ref_name]]
    ref_genes <- ref_genes[ref_genes %in% background_genes]
    M <- length(ref_genes)
    if (M == 0) { cl_res[cl_res$Reference == ref_name, "P_value"] <- 1; next }
    overlap_genes <- intersect(cluster_genes, ref_genes)
    x <- length(overlap_genes)
    cl_res[cl_res$Reference == ref_name, "Overlap"] <- x
    p_val <- phyper(x - 1, M, N - M, k, lower.tail = FALSE)
    cl_res[cl_res$Reference == ref_name, "P_value"] <- p_val
  }
  cl_res$P_adjust <- p.adjust(cl_res$P_value, method = "fdr")
  enrichment_results[[cl]] <- cl_res
}
enrichment_df <- do.call(rbind, enrichment_results)
write.csv(enrichment_df, file.path(out_dir, "cluster_enrichment_scores.csv"), row.names = FALSE)

# ----- 8. 分配最终标签 -----
cat("\n🏷️ 分配最终标签...\n")
final_labels <- setNames(rep("Unknown", new_clusters), as.character(0:(new_clusters-1)))

for (cl in names(cluster_genes_list)) {
  cl_data <- enrichment_results[[cl]]
  sig <- cl_data[cl_data$P_adjust < 0.05 & cl_data$Overlap >= 2, ]
  if (nrow(sig) == 0) {
    final_labels[cl] <- "Unknown"
    next
  }
  # 中性粒细胞优先
  if ("Neutrophil_Custom" %in% sig$Reference) {
    neutro_sig <- sig[sig$Reference == "Neutrophil_Custom", ]
    if (neutro_sig$P_adjust < 0.01) {
      final_labels[cl] <- "Neutrophil"
      next
    }
  }
  best_hit <- sig[which.min(sig$P_adjust), ]
  label <- best_hit$Reference
  if (label == "Neutrophil_Custom") label <- "Neutrophil"
  final_labels[cl] <- label
}
cat("   ✅ 标签分配完成\n")
print(table(final_labels))

# ----- 9. 映射到Seurat对象（先映射，再补充标注）-----
cat("\n🔄 映射到Seurat对象...\n")
cluster_ids <- as.character(sc_obj$seurat_clusters)
annot_vec <- final_labels[cluster_ids]
names(annot_vec) <- colnames(sc_obj)
sc_obj <- AddMetaData(sc_obj, metadata = annot_vec, col.name = "reference_annotation")
cat("   ✅ 映射完成\n")

# ----- 9.1 补充标注：处理Unknown中的肿瘤细胞（此时reference_annotation已存在）-----
cat("\n🔧 补充标注：处理Unknown中的肿瘤细胞...\n")

# 定义通用肿瘤标志物
general_tumor_markers <- c(
  # 核心干性
  "SOX2", "OLIG2", "NES", "PROM1", "STAT3",
  # 增殖
  "MKI67", "TOP2A",
  # 亚型相关 (经典型、间充质型)
  "EGFR", "CD44", "CHI3L1", "VIM", "MET", "PDGFRA",
  # 微环境与侵袭
  "FN1", "MMP9", "TNC",
  # 其他
  "HIF1A", "TP53")
general_tumor_markers <- general_tumor_markers[general_tumor_markers %in% rownames(sc_obj)]

# 计算每个细胞的肿瘤标志物平均表达
tumor_marker_score <- colMeans(GetAssayData(sc_obj, assay = "RNA", layer = "data")[general_tumor_markers, ])
threshold <- median(tumor_marker_score, na.rm = TRUE)

# 找出当前为Unknown的细胞
unknown_cells <- colnames(sc_obj)[sc_obj$reference_annotation == "Unknown"]
if (length(unknown_cells) > 0) {
  # 检查这些细胞是否高表达肿瘤标志物
  for (cell in unknown_cells) {
    if (tumor_marker_score[cell] > threshold) {
      sc_obj$reference_annotation[cell] <- "General_Tumor"
    }
  }
  cat("   ✅ 已将", sum(tumor_marker_score[unknown_cells] > threshold), "个Unknown细胞标记为 General_Tumor\n")
} else {
  cat("   ✅ 没有Unknown细胞需要处理\n")
}

# 查看更新后的分布
cat("\n📊 更新后的参考注释分布:\n")
print(table(sc_obj$reference_annotation, useNA = "ifany"))

# ----- 10. 诊断可视化 -----
# 10.1 UMAP
pdf(file.path(out_dir, "UMAP_Reference_Annotation.pdf"), width = 12, height = 8)
p <- DimPlot(sc_obj, group.by = "reference_annotation", label = TRUE, repel = TRUE) +
  ggtitle("Reference Annotation (based on Top50 marker enrichment)")
print(p)
dev.off()

# 10.2 重叠热图
overlap_matrix <- matrix(0, nrow = length(cluster_genes_list), ncol = length(ref_list),
                         dimnames = list(names(cluster_genes_list), names(ref_list)))
for (cl in names(cluster_genes_list)) {
  cl_data <- enrichment_results[[cl]]
  for (ref in names(ref_list)) {
    overlap_matrix[cl, ref] <- cl_data[cl_data$Reference == ref, "Overlap"]
  }
}
pdf(file.path(out_dir, "Cluster_Reference_Overlap_Heatmap.pdf"), width = 14, height = 10)
pheatmap(overlap_matrix, 
         main = "Overlap count between cluster Top50 and reference sets",
         color = colorRampPalette(c("white", "steelblue", "darkred"))(100),
         display_numbers = TRUE,
         fontsize_number = 6,
         cluster_rows = TRUE,
         cluster_cols = TRUE)
dev.off()

# ----- 11. 导出用于RCTD的Reference（修复GetAssayData语法）-----
# Seurat v5 使用 layer 参数
ref_counts <- GetAssayData(sc_obj, assay = "RNA", layer = "data")
ref_meta <- sc_obj@meta.data
keep_cells <- rownames(ref_meta)[!ref_meta$reference_annotation %in% c("Unknown", NA)]
ref_counts_filtered <- ref_counts[, keep_cells]
ref_meta_filtered <- ref_meta[keep_cells, ]

saveRDS(list(counts = ref_counts_filtered, metadata = ref_meta_filtered), 
        file = file.path(out_dir, "RCTD_Reference_GBM_MPs.rds"))

# 保存完整对象
saveRDS(sc_obj, file = file.path(out_dir, "sc_obj_annotated_with_reference.rds"))

cat("\n========================================\n")
cat("🎉 参考注释构建完成！\n")
cat("📁 输出目录:", out_dir, "\n")
cat("   - UMAP_Reference_Annotation.pdf\n")
cat("   - Cluster_Reference_Overlap_Heatmap.pdf\n")
cat("   - cluster_enrichment_scores.csv\n")
cat("   - RCTD_Reference_GBM_MPs.rds\n")
cat("   - sc_obj_annotated_with_reference.rds\n")
cat("========================================\n")
# ============================================================
#  完整 RCTD 解卷积（使用修正合并逻辑）并绘图
#  适用于所有 5 个样本
# ============================================================

cat("\n🔄 开始对全部样本进行 RCTD 解卷积（修正版）...\n")

# 循环处理每个样本
for (s in sample_names) {
  cat("\n========================================\n")
  cat("处理样本 (RCTD 修正):", s, "\n")
  
  spatial_obj <- seurat_list[[s]]
  
  # 提取数据和坐标
  spatial_counts <- GetAssayData(spatial_obj, assay = "Spatial", layer = "counts")
  coords <- GetTissueCoordinates(spatial_obj, cols = c("imagerow", "imagecol"))
  if (ncol(coords) > 2) coords <- coords[, 1:2, drop = FALSE]
  colnames(coords)[1:2] <- c("x", "y")
  if (is.null(rownames(coords))) rownames(coords) <- colnames(spatial_obj)
  
  spatial_rna <- SpatialRNA(coords, spatial_counts)
  cat("   空间数据包含", ncol(spatial_rna@counts), "个spots\n")
  
  # 运行 RCTD
  myRCTD <- create.RCTD(spatial_rna, reference, max_cores = 1)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")
  
  weights <- myRCTD@results$weights
  cat("   ✅ RCTD 完成，丰度矩阵维度:", dim(weights), "\n")
  
  # ----- 修正合并逻辑（已验证有效） -----
  weights_df <- as.data.frame(as.matrix(weights))
  common_barcodes <- intersect(rownames(weights_df), colnames(spatial_obj))
  weights_df <- weights_df[common_barcodes, , drop = FALSE]
  
  for (ct in colnames(weights_df)) {
    col_name <- paste0("proportion_", ct)
    vec <- rep(NA, ncol(spatial_obj))
    names(vec) <- colnames(spatial_obj)
    vec[common_barcodes] <- weights_df[common_barcodes, ct]
    spatial_obj[[col_name]] <- vec
  }
  
  # 验证（打印非 NA 数量）
  test_col <- paste0("proportion_", colnames(weights_df)[1])
  if (test_col %in% colnames(spatial_obj@meta.data)) {
    vals <- spatial_obj[[test_col]][, 1]
    cat("   ✅ 写入验证: 非NA值数 =", sum(!is.na(vals)), "/", length(vals), "\n")
  }
  
  seurat_list[[s]] <- spatial_obj
  cat("   ✅ 样本", s, "解卷积完成\n")
}

# 保存完整结果（覆盖之前的文件）
saveRDS(seurat_list, file = file.path(base_dir, "seurat_list_with_deconvolution_rctd_corrected.rds"))
cat("\n💾 修正后的结果已保存至:", file.path(base_dir, "seurat_list_with_deconvolution_rctd_corrected.rds"), "\n")

# ----- 重新生成空间丰度图（PDF） -----
cat("\n📄 生成细胞类型空间丰度图（PDF）...\n")

for (s in sample_names) {
  obj <- seurat_list[[s]]
  prop_cols <- grep("^proportion_", colnames(obj@meta.data), value = TRUE)
  
  if (length(prop_cols) == 0) {
    cat("   ⚠️", s, "无比例列，跳过绘图\n")
    next
  }
  
  # 按平均丰度排序，取前 6 种细胞类型
  avg_prop <- sapply(prop_cols, function(x) mean(obj[[x]][, 1], na.rm = TRUE))
  top_cols <- names(sort(avg_prop, decreasing = TRUE))[1:min(6, length(prop_cols))]
  cat("   ", s, "绘制:", paste(gsub("proportion_", "", top_cols), collapse = ", "), "\n")
  
  # 绘制空间特征图（淡背景）
  p <- SpatialFeaturePlot(
    obj,
    features = top_cols,
    ncol = 3,
    pt.size.factor = 1.5,
    alpha = c(0.1, 1),
    image.alpha = 0.15
  ) +
    patchwork::plot_annotation(title = paste0(s, " - RCTD Cell Type Proportions"))
  
  pdf_file <- file.path(result_dir, paste0("Deconv_RCTD_corrected_", s, ".pdf"))
  ggsave(
    filename = pdf_file,
    plot = p,
    width = 15,
    height = 6 * ceiling(length(top_cols) / 3),
    dpi = 300,
    device = "pdf"
  )
  cat("   ✅ 已保存:", pdf_file, "\n")
}

cat("\n========================================\n")
cat("🎉 所有样本 RCTD 解卷积（修正版）完成！\n")
cat("请查看结果文件夹中的 Deconv_RCTD_corrected_*.pdf 文件。\n")
cat("========================================\n")

# ============================================================
#  检查并修正Reference（确保使用原始counts）
#  输出：修正后的 reference_list.rds
# ============================================================

rm(list = ls())
library(Seurat)

ref_dir <- "E:/NETS/空转组数据/结果/Reference_Annotation"
out_dir <- "E:/NETS/空转组数据/结果/RCTD_Analysis"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ----- 1. 加载完整的注释对象 -----
cat("📂 加载注释对象...\n")
sc_obj <- readRDS(file.path(ref_dir, "sc_obj_annotated_with_reference.rds"))

# ----- 2. 提取原始counts（整数）-----
cat("📊 提取原始counts（layer = counts）...\n")
ref_counts <- GetAssayData(sc_obj, assay = "RNA", layer = "counts")

# 过滤掉Unknown和NA
keep_cells <- colnames(ref_counts)[!sc_obj$reference_annotation %in% c("Unknown", NA)]
ref_counts <- ref_counts[, keep_cells]
ref_labels <- sc_obj$reference_annotation[keep_cells]
ref_labels <- factor(ref_labels, levels = unique(ref_labels))

cat("   ✅ Reference cells:", length(ref_labels), "\n")
cat("   ✅ Cell types:\n")
print(table(ref_labels))

# ----- 3. 验证是否为整数-----
if (all(ref_counts == round(ref_counts))) {
  cat("   ✅ 确认：counts为整数矩阵\n")
} else {
  cat("   ❌ 警告：counts不是整数！请检查数据来源\n")
}

# ----- 4. 保存修正后的Reference-----
reference <- list(
  counts = ref_counts,
  cell_types = ref_labels,
  n_cells = ncol(ref_counts),
  n_types = length(unique(ref_labels)),
  type_counts = table(ref_labels)
)

saveRDS(reference, file = file.path(out_dir, "reference_for_RCTD.rds"))
cat("\n💾 已保存:", file.path(out_dir, "reference_for_RCTD.rds"), "\n")

# ============================================================
#  RCTD 空间解卷积（使用修正后的Reference）
#  整合你的原始代码风格 + 新的Reference
# ============================================================

rm(list = ls())
library(Seurat)
library(spacexr)          # RCTD包
library(ggplot2)
library(patchwork)

# ----- 1. 路径设置 -----
base_dir <- "E:/NETS/空转组数据/原数据"
result_dir <- "E:/NETS/空转组数据/结果"
rctd_dir <- file.path(result_dir, "RCTD_Analysis")
if (!dir.exists(rctd_dir)) dir.create(rctd_dir, recursive = TRUE)

# ----- 2. 加载空间数据（含NETs评分）-----
cat("📂 加载空间数据...\n")
seurat_list <- readRDS(file.path(base_dir, "seurat_list_with_NETs_scores.rds"))
sample_names <- names(seurat_list)
cat("   ✅ 加载了", length(sample_names), "个样本\n")

# ----- 3. 加载修正后的Reference（原始counts）-----
cat("\n📂 加载Reference...\n")
reference <- readRDS(file.path(rctd_dir, "reference_for_RCTD.rds"))
ref_counts <- reference$counts
ref_labels <- reference$cell_types
cat("   ✅ Reference:", ncol(ref_counts), "个细胞,", length(unique(ref_labels)), "种类型\n")
print(table(ref_labels))

# 创建RCTD Reference对象
reference_obj <- Reference(
  counts = ref_counts,
  cell_types = ref_labels
)

# ----- 4. RCTD解卷积（循环处理每个样本）-----
cat("\n🔄 开始对全部样本进行 RCTD 解卷积...\n")

for (s in sample_names) {
  cat("\n========================================\n")
  cat("处理样本:", s, "\n")
  
  spatial_obj <- seurat_list[[s]]
  
  # ----- 4.1 提取原始counts（整数）-----
  spatial_counts <- GetAssayData(spatial_obj, assay = "Spatial", layer = "counts")
  
  # ----- 4.2 提取坐标-----
  coords <- GetTissueCoordinates(spatial_obj, cols = c("imagerow", "imagecol"))
  if (is.null(coords) || nrow(coords) == 0) {
    # 尝试其他列名
    coords <- GetTissueCoordinates(spatial_obj)
  }
  if (is.null(coords)) {
    cat("   ❌ 无法提取坐标，跳过此样本\n")
    next
  }
  
  # 确保只有两列x,y
  if (ncol(coords) >= 2) {
    coords <- coords[, 1:2, drop = FALSE]
  }
  colnames(coords)[1:2] <- c("x", "y")
  if (is.null(rownames(coords))) rownames(coords) <- colnames(spatial_obj)
  
  cat("   空间数据:", ncol(spatial_counts), "个spots,", nrow(spatial_counts), "个基因\n")
  
  # ----- 4.3 创建SpatialRNA对象-----
  spatial_rna <- SpatialRNA(
    coords = coords,
    counts = spatial_counts
  )
  
  # ----- 4.4 运行RCTD-----
  cat("   ⏳ 运行RCTD（可能需要几分钟）...\n")
  myRCTD <- create.RCTD(
    spatialRNA = spatial_rna,
    reference = reference_obj,
    max_cores = 1,
    test_mode = FALSE
  )
  myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")
  
  # 提取权重矩阵（每个spot的细胞类型丰度）
  weights <- myRCTD@results$weights
  cat("   ✅ RCTD完成，丰度矩阵维度:", dim(weights), "\n")
  
  # ----- 4.5 写入Seurat对象（你的原始方法，已验证有效）-----
  weights_df <- as.data.frame(as.matrix(weights))
  common_barcodes <- intersect(rownames(weights_df), colnames(spatial_obj))
  weights_df <- weights_df[common_barcodes, , drop = FALSE]
  
  for (ct in colnames(weights_df)) {
    col_name <- paste0("proportion_", ct)
    vec <- rep(NA, ncol(spatial_obj))
    names(vec) <- colnames(spatial_obj)
    vec[common_barcodes] <- weights_df[common_barcodes, ct]
    spatial_obj[[col_name]] <- vec
  }
  
  # 验证
  test_col <- paste0("proportion_", colnames(weights_df)[1])
  if (test_col %in% colnames(spatial_obj@meta.data)) {
    vals <- spatial_obj[[test_col]][, 1]
    cat("   ✅ 写入验证: 非NA值数 =", sum(!is.na(vals)), "/", length(vals), "\n")
  }
  
  # 更新列表
  seurat_list[[s]] <- spatial_obj
  
  # ----- 4.6 保存单个RCTD结果-----
  saveRDS(myRCTD, file = file.path(rctd_dir, paste0("RCTD_result_", s, ".rds")))
  cat("   💾 已保存:", file.path(rctd_dir, paste0("RCTD_result_", s, ".rds")), "\n")
}

# ----- 5. 保存完整结果-----
saveRDS(seurat_list, file = file.path(rctd_dir, "seurat_list_with_deconvolution_rctd.rds"))
cat("\n💾 完整解卷积结果已保存至:", file.path(rctd_dir, "seurat_list_with_deconvolution_rctd.rds"), "\n")

# ----- 6. 生成空间丰度图-----
cat("\n📄 生成细胞类型空间丰度图（PDF）...\n")

for (s in sample_names) {
  obj <- seurat_list[[s]]
  prop_cols <- grep("^proportion_", colnames(obj@meta.data), value = TRUE)
  
  if (length(prop_cols) == 0) {
    cat("   ⚠️", s, "无比例列，跳过绘图\n")
    next
  }
  
  # 按平均丰度排序，取前6种细胞类型
  avg_prop <- sapply(prop_cols, function(x) mean(obj[[x]][, 1], na.rm = TRUE))
  top_cols <- names(sort(avg_prop, decreasing = TRUE))[1:min(6, length(prop_cols))]
  
  # 提取短名称（去掉"proportion_"前缀）
  short_names <- gsub("proportion_", "", top_cols)
  cat("   ", s, "绘制:", paste(short_names, collapse = ", "), "\n")
  
  # 绘制空间特征图（淡背景）
  p <- SpatialFeaturePlot(
    obj,
    features = top_cols,
    ncol = 3,
    pt.size.factor = 1.5,
    alpha = c(0.1, 1),
    image.alpha = 0.15
  ) +
    patchwork::plot_annotation(title = paste0(s, " - RCTD Cell Type Proportions"))
  
  pdf_file <- file.path(rctd_dir, paste0("Deconv_RCTD_", s, ".pdf"))
  ggsave(
    filename = pdf_file,
    plot = p,
    width = 15,
    height = 6 * ceiling(length(top_cols) / 3),
    dpi = 300,
    device = "pdf"
  )
  cat("   ✅ 已保存:", pdf_file, "\n")
}

cat("\n========================================\n")
cat("🎉 所有样本 RCTD 解卷积完成！\n")
cat("📁 输出目录:", rctd_dir, "\n")
cat("   - Deconv_RCTD_*.pdf (空间丰度图)\n")
cat("   - RCTD_result_*.rds (原始RCTD输出)\n")
cat("   - seurat_list_with_deconvolution_rctd.rds (完整数据)\n")
cat("========================================\n")

# ============================================================
#  独立可视化脚本：MP硬分类六边形图 + RCTD全丰度单图 + 汇总表格
#  运行前提：已完成RCTD，存在 seurat_list_with_deconvolution_rctd.rds
#  自动检测并计算MP打分（如果未计算）
# ============================================================

rm(list = ls())
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(readxl)
library(FNN)

# ----- 1. 路径设置 -----
base_dir <- "E:/NETS/空转组数据/原数据"
result_dir <- "E:/NETS/空转组数据/结果"
rctd_dir <- file.path(result_dir, "RCTD_Analysis")
if (!dir.exists(rctd_dir)) dir.create(rctd_dir, recursive = TRUE)

# MP基因集文件
mp_file <- file.path(base_dir, "mmc2.xlsx")

# ----- 2. 加载RCTD解卷积结果 -----
cat("📂 加载RCTD解卷积结果...\n")
seurat_list <- readRDS(file.path(rctd_dir, "seurat_list_with_deconvolution_rctd.rds"))
sample_names <- names(seurat_list)
cat("   ✅ 加载了", length(sample_names), "个样本\n")

# ----- 3. 检测并计算MP打分（如果尚未计算）-----
cat("\n🔍 检查MP打分是否存在...\n")
first_obj <- seurat_list[[sample_names[1]]]

if (!"dominant_MP" %in% colnames(first_obj@meta.data)) {
  cat("   ⚠️ 未检测到MP打分，开始计算...\n")
  
  # 读取MP基因集
  if (!file.exists(mp_file)) {
    stop("MP基因集文件不存在: ", mp_file)
  }
  mp_df <- read_excel(mp_file, sheet = "GBM_MPs")
  mp_list <- list()
  for (col_name in colnames(mp_df)) {
    genes <- mp_df[[col_name]]
    genes <- genes[!is.na(genes) & genes != ""]
    clean_name <- gsub(" \\(malig\\.\\)", "", col_name)
    clean_name <- gsub(" ", "_", clean_name)
    clean_name <- gsub("\\.", "_", clean_name)
    mp_list[[clean_name]] <- genes
  }
  cat("   ✅ 加载了", length(mp_list), "个GBM MPs\n")
  
  # 对每个样本计算MP打分
  for (s in names(seurat_list)) {
    obj <- seurat_list[[s]]
    cat("   🔄 计算MP打分:", s, "\n")
    
    # 清除旧MP列（避免重复）
    old_mp <- grep("^MP_|^MPscore_|^dominant_MP", colnames(obj@meta.data), value = TRUE)
    if (length(old_mp) > 0) {
      obj@meta.data <- obj@meta.data[, !colnames(obj@meta.data) %in% old_mp]
    }
    
    for (mp_name in names(mp_list)) {
      genes <- mp_list[[mp_name]]
      genes <- genes[genes %in% rownames(obj)]
      if (length(genes) > 2) {
        obj <- AddModuleScore(
          obj,
          features = list(genes),
          name = paste0("MP_", mp_name, "_"),
          assay = "Spatial",
          search = TRUE
        )
      }
    }
    
    # 提取MP打分列，计算dominant
    mp_score_cols <- grep("^MP_", colnames(obj@meta.data), value = TRUE)
    if (length(mp_score_cols) > 0) {
      score_mat <- as.matrix(obj@meta.data[, mp_score_cols])
      clean_names <- gsub("^MP_", "", mp_score_cols)
      clean_names <- gsub("_1$", "", clean_names)
      colnames(score_mat) <- clean_names
      
      max_idx <- apply(score_mat, 1, which.max)
      obj$dominant_MP <- clean_names[max_idx]
      obj$MP_max_score <- apply(score_mat, 1, max, na.rm = TRUE)
      
      for (col in clean_names) {
        obj[[paste0("MPscore_", col)]] <- score_mat[, col]
      }
    }
    seurat_list[[s]] <- obj
  }
  
  # 保存带MP打分的版本
  saveRDS(seurat_list, file = file.path(rctd_dir, "seurat_list_with_RCTD_and_MPscores.rds"))
  cat("   💾 已保存带MP打分的完整数据\n")
} else {
  cat("   ✅ 已存在MP打分，直接使用\n")
}
# 保存MP基因列表（紧接在读取mp_df之后，计算打分之前）
saveRDS(mp_list, file = file.path(rctd_dir, "MP_gene_list.rds"))

# ----- 4. 定义颜色（14种MP + 12种RCTD，实色）-----
mp_colors <- c(
  "Neuron"          = "#6A3D9A", "Vasc"            = "#1B9E77",
  "MES_Hyp"         = "#1F78B4", "Mac"             = "#7570B3",
  "Oligo"           = "#FDBF6F", "MES"             = "#E41A1C",
  "Prolif_Metab"    = "#D95F02", "MES_Ast"         = "#B15928",
  "Reactive_Ast"    = "#8B4513", "NPC"             = "#A6CEE3",
  "InflammatoryMac" = "#E7298A", "Chromatin_reg"   = "#33A02C",
  "OPC"             = "#CAB2D6", "AC"              = "#2E8B57"
)

celltype_colors <- c(
  "InflammatoryMac" = "#D95F02", "Vasc"  = "#1B9E77", "Mac"   = "#7570B3",
  "General_Tumor"   = "#E7298A", "MES"   = "#E41A1C", "Neutrophil" = "#FF7F00",
  "Neuron"          = "#6A3D9A", "OPC"   = "#A6CEE3", "Oligo" = "#FDBF6F",
  "MES_Ast"         = "#B15928", "AC"    = "#33A02C", "MES_Hyp" = "#1F78B4"
)

# ----- 5. MP六边形硬分类图（全图例）-----
cat("\n🟢 生成MP六边形硬分类图（全图例）...\n")

for (s in names(seurat_list)) {
  obj <- seurat_list[[s]]
  if (!"dominant_MP" %in% colnames(obj@meta.data)) next
  
  coords <- GetTissueCoordinates(obj, scale = "lowres")
  if (is.null(coords) || nrow(coords) == 0) {
    coords <- GetTissueCoordinates(obj, scale = "tissue")
  }
  if (is.null(coords)) next
  
  # 坐标列名
  if (all(c("imagerow", "imagecol") %in% colnames(coords))) {
    xcol <- "imagerow"; ycol <- "imagecol"
  } else if (all(c("row", "col") %in% colnames(coords))) {
    xcol <- "row"; ycol <- "col"
  } else {
    xcol <- colnames(coords)[1]; ycol <- colnames(coords)[2]
  }
  
  df <- coords
  df$dominant <- obj$dominant_MP
  df <- df[!is.na(df$dominant), ]
  if (nrow(df) == 0) next
  
  # 六边形参数
  coords_mat <- as.matrix(df[, c(xcol, ycol)])
  if (nrow(coords_mat) < 3) next
  nn_dist <- FNN::knn.dist(coords_mat, k = 2)[,2]
  side <- median(nn_dist, na.rm = TRUE) / sqrt(3) * 0.85
  if (is.na(side) || side < 0.1) side <- 5
  
  angles <- (60 * (0:5)) * pi / 180
  hex_list <- list()
  for (i in 1:nrow(df)) {
    xc <- df[[xcol]][i]; yc <- df[[ycol]][i]
    vertices <- data.frame(
      x = xc + side * cos(angles),
      y = yc + side * sin(angles)
    )
    vertices$group <- i
    vertices$dominant <- df$dominant[i]
    hex_list[[i]] <- vertices
  }
  hex_df <- do.call(rbind, hex_list)
  
  p <- ggplot(hex_df, aes(x = x, y = y, fill = dominant, group = group)) +
    geom_polygon(color = NA, size = 0) +
    coord_fixed() +
    scale_fill_manual(values = mp_colors, name = "Dominant MP", drop = FALSE) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.key.width = unit(0.6, "cm"),
      legend.key.height = unit(0.5, "cm"),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 8),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    labs(title = paste0(s, " - MP Dominant (Hex)"))
  
  pdf_file <- file.path(rctd_dir, paste0("Hex_MP_Dominant_", s, ".pdf"))
  ggsave(pdf_file, p, width = 10, height = 8, dpi = 300, device = "pdf")
  cat("   ✅ 已保存:", basename(pdf_file), "\n")
}

# ----- 6. RCTD全12种细胞丰度单图（一页全部展示）-----
cat("\n🟢 生成RCTD全12种细胞丰度单图...\n")

for (s in names(seurat_list)) {
  obj <- seurat_list[[s]]
  prop_cols <- grep("^proportion_", colnames(obj@meta.data), value = TRUE)
  if (length(prop_cols) == 0) next
  
  # 按颜色顺序排序
  ct_order <- intersect(names(celltype_colors), gsub("proportion_", "", prop_cols))
  prop_cols <- paste0("proportion_", ct_order)
  
  plot_list <- list()
  for (col_name in prop_cols) {
    ct_name <- gsub("proportion_", "", col_name)
    max_val <- max(obj[[col_name]][, 1], na.rm = TRUE)
    if (is.na(max_val) || max_val == 0) max_val <- 1
    
    p <- SpatialFeaturePlot(
      obj,
      features = col_name,
      pt.size.factor = 3.2,          # 【改】点的大小，原1.2 → 1.8
      alpha = c(0.1, 1),
      image.alpha = 0.15,
      stroke = 0
    ) +
      scale_fill_gradient(
        low = "white",
        high = celltype_colors[ct_name],
        limits = c(0, max_val),
        name = "Proportion",
        guide = guide_colorbar(
          frame.colour = "black",
          ticks.colour = "black",
          barwidth = 10,            # 【改】标尺宽度，原0.8 → 1.2
          barheight = 1              # 【改】标尺高度，原4 → 6
        )
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        legend.position = "bottom"
      ) +
      ggtitle(ct_name)
    
    plot_list[[length(plot_list) + 1]] <- p
  }
  
  n_plots <- length(plot_list)
  if (n_plots == 0) next
  n_col <- 4
  n_row <- ceiling(n_plots / n_col)
  
  combined <- wrap_plots(plot_list, ncol = n_col) +
    plot_annotation(
      title = paste0(s, " - RCTD Full Proportions"),
      theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    )
  
  pdf_file <- file.path(rctd_dir, paste0("RCTD_Full12_", s, ".pdf"))
  # 【改】适当增大画布，避免点挤在一起
  ggsave(pdf_file, combined, width = 6 * n_col, height = 6 * n_row, dpi = 300, device = "pdf")
  cat("   ✅ 已保存:", basename(pdf_file), "\n")
}

# ----- 7. 汇总表格（MP平均打分 + RCTD平均丰度）-----
cat("\n🟢 生成汇总表格...\n")

summary_list <- list()
for (s in names(seurat_list)) {
  obj <- seurat_list[[s]]
  row_data <- data.frame(Sample = s)
  
  mp_cols <- grep("^MPscore_", colnames(obj@meta.data), value = TRUE)
  if (length(mp_cols) > 0) {
    mp_avg <- colMeans(obj@meta.data[, mp_cols], na.rm = TRUE)
    names(mp_avg) <- gsub("MPscore_", "MP_", names(mp_avg))
    row_data <- cbind(row_data, t(mp_avg))
  }
  
  prop_cols <- grep("^proportion_", colnames(obj@meta.data), value = TRUE)
  if (length(prop_cols) > 0) {
    prop_avg <- colMeans(obj@meta.data[, prop_cols], na.rm = TRUE)
    names(prop_avg) <- gsub("proportion_", "RCTD_", names(prop_avg))
    row_data <- cbind(row_data, t(prop_avg))
  }
  
  summary_list[[s]] <- row_data
}

summary_df <- do.call(rbind, summary_list)
if ("MP_MES_Hyp" %in% colnames(summary_df)) {
  summary_df <- summary_df[order(summary_df$MP_MES_Hyp, decreasing = TRUE), ]
}
rownames(summary_df) <- NULL

csv_file <- file.path(rctd_dir, "Summary_MP_RCTD_Averages.csv")
write.csv(summary_df, csv_file, row.names = FALSE)
cat("   ✅ 汇总表格已保存:", csv_file, "\n")

cat("\n📊 汇总表格预览（前6行）:\n")
print(head(summary_df))

cat("\n========================================\n")
cat("🎉 所有输出生成完毕！\n")
cat("📁 输出目录:", rctd_dir, "\n")
cat("   ├── Hex_MP_Dominant_*.pdf   (MP六边形硬分类，全图例)\n")
cat("   ├── RCTD_Full12_*.pdf       (RCTD全12种细胞丰度)\n")
cat("   └── Summary_MP_RCTD_Averages.csv (汇总表格)\n")
cat("========================================\n")


