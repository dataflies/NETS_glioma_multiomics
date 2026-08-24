# 加载必要的包
library(Seurat)
library(tidyverse)
library(patchwork)

# 设置路径
data_dir <- "E:/NETS/单细胞数据/原数据/原1-CGGA.scRNAseq_6148.count.matrix.tsv/CGGA.scRNAseq_6148.count.matrix.tsv"
clinical_file <- "E:/NETS/单细胞数据/原数据/原1-CGGA.scRNAseq_6148_clinical.20220104.txt/CGGA.scRNA_6148_clinical.20220104.txt"
out_dir <- "E:/NETS/单细胞数据/结果"

# 创建输出目录
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ------------------------------
# 1. 读取数据
# ------------------------------
cat("正在读取表达矩阵...\n")
# 读取表达矩阵（基因名作为行名，细胞名作为列名，制表符分隔）
counts <- read.table(data_dir, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
cat("表达矩阵维度：", dim(counts), "\n")

cat("正在读取临床信息...\n")
clinical <- read.table(clinical_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cat("临床信息维度：", dim(clinical), "\n")
print(head(clinical))

# ------------------------------
# 2. 创建Seurat对象
# ------------------------------
cat("创建Seurat对象...\n")
# 注意：Seurat要求表达矩阵行为基因，列为细胞。counts已经满足。
seurat_obj <- CreateSeuratObject(counts = counts, project = "CGGA_scRNA", min.cells = 3, min.features = 200)

# 查看初始过滤后细胞数
cat("初始过滤后细胞数：", ncol(seurat_obj), "\n")

# ------------------------------
# 3. 添加临床样本信息（完整版）
# ------------------------------
cat("添加临床信息（基于 Sample.point）...\n")

# 提取细胞名中的前缀（例如 "S2P1"）
cell_names <- colnames(seurat_obj)
cell_prefix <- str_extract(cell_names, "^S\\d+P\\d+")
# 加上 "G" 得到临床 Sample.point 格式
sample_point <- paste0("G", cell_prefix)

cat("示例匹配：", head(sample_point), "\n")

# 从临床文件中提取需要的列
# 首先确认 sample.type 列的确切名称（从您的输出看是 sample.type）
# 如果列名包含点，可能是 "sample.type"，直接使用
clinical_sub <- clinical %>% 
  select(Sample.point, Gender, Age, WHO.Grade, Histology, IDH1.R132.mutation, MGMT.promoter.methylation) %>%
  distinct(Sample.point, .keep_all = TRUE)

# 为每个细胞添加基础临床信息
meta <- data.frame(row.names = cell_names, sample_point = sample_point, stringsAsFactors = FALSE)
meta <- meta %>% 
  left_join(clinical_sub, by = c("sample_point" = "Sample.point"))

# 添加到Seurat对象
seurat_obj <- AddMetaData(seurat_obj, metadata = meta)

# 单独添加 sample.type 列（使用原始临床数据，不进行去重，因为每个 sample.point 对应唯一 sample.type）
# 建立 sample.point -> sample.type 的映射（取第一个匹配值，确保唯一）
sample_type_map <- clinical %>% 
  select(Sample.point, Sample.type) %>%
  distinct(Sample.point, .keep_all = TRUE)

# 根据已有的 sample_point 列匹配 sample.type
sample_type_for_cells <- sample_type_map$Sample.type[match(seurat_obj$sample_point, sample_type_map$Sample.point)]

# 添加到 Seurat 对象
seurat_obj$Sample.type <- sample_type_for_cells

# 检查结果
cat("sample.type 匹配结果（前10个）：\n")
print(head(seurat_obj$Sample.type, 10))
cat("sample.type 分布：\n")
print(table(seurat_obj$Sample.type, useNA = "ifany"))
# ------------------------------
# 4. 质控：线粒体基因比例
# ------------------------------
cat("计算线粒体基因比例...\n")
# 假设人类基因，线粒体基因以"MT-"开头（注意大小写）
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# 可视化质控指标（可选，保存为PDF）
pdf(file.path(out_dir, "qc_violin.pdf"), width = 12, height = 6)
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()

# 过滤细胞：保留基因数200-6000，线粒体比例<20%
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
cat("质控后细胞数：", ncol(seurat_obj), "\n")
# 查看基因名中可能包含 "MT" 的（不区分大小写）
mt_genes <- grep("MT", rownames(seurat_obj), value = TRUE, ignore.case = TRUE)
head(mt_genes, 20)
# ------------------------------
# 5. 标准化与高变基因
# ------------------------------
cat("标准化...\n")
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 1e4)

cat("寻找高变基因...\n")
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# 展示高变基因（可选）
top10 <- head(VariableFeatures(seurat_obj), 10)
pdf(file.path(out_dir, "variable_features.pdf"), width = 8, height = 6)
plot1 <- VariableFeaturePlot(seurat_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()

# ------------------------------
# 6. 数据缩放与PCA
# ------------------------------
cat("数据缩放...\n")
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)

cat("PCA降维...\n")
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), npcs = 50)

# 可视化PCA
pdf(file.path(out_dir, "pca_elbow.pdf"), width = 6, height = 4)
ElbowPlot(seurat_obj, ndims = 50)
dev.off()

# ------------------------------
# 7. 聚类与UMAP
# ------------------------------
cat("计算邻域图...\n")
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.6)

cat("UMAP降维...\n")
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, min.dist = 0.5, spread = 1)

# 绘制UMAP
pdf(file.path(out_dir, "umap_clusters.pdf"), width = 8, height = 6)
DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE) + ggtitle("UMAP by clusters")
dev.off()

# ------------------------------
# 8. 细胞类型注释（基于差异基因热图 + 手动判断）
# ------------------------------

# 确保聚类身份为 seurat_clusters
Idents(seurat_obj) <- "seurat_clusters"

# 1. 计算每个簇的差异表达基因（如果尚未计算，重新运行）
#   如果已经存在 markers 对象，可跳过此步；这里为保险重新计算
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, file.path(out_dir, "cluster_markers.csv"), row.names = FALSE)

# 2. 提取每个簇的 top 10 差异基因（按 avg_log2FC 排序）
library(dplyr)
top10_per_cluster <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  ungroup()

# 查看 top10 基因列表
print(top10_per_cluster %>% select(cluster, gene, avg_log2FC, pct.1), n = 200)

# 3. 获取这些 top 基因在簇中的平均表达（用于热图）
top_genes <- unique(top10_per_cluster$gene)
avg_exp_top <- AverageExpression(seurat_obj, features = top_genes, assay = "RNA", slot = "data")$RNA

# 4. 绘制热图（按行缩放，展示相对表达）
library(pheatmap)
# 第八步：绘制 top10 标记基因热图（优化版）
library(pheatmap)
library(RColorBrewer)

# 准备表达矩阵（已经在代码中计算好 avg_exp_top）
# 确保 avg_exp_top 是一个矩阵，行为基因，列为细胞群
avg_exp_top <- AverageExpression(seurat_obj, features = top_genes, assay = "RNA", slot = "data")$RNA

# 设置美观的颜色方案（蓝白红渐变，适合展示相对表达）
my_colors <- colorRampPalette(c("#2166ac", "white", "#b2182b"))(100)

# 绘制并保存热图
pdf(file.path(out_dir, "top10_marker_heatmap_optimized.pdf"), width = 14, height = 12)

pheatmap(avg_exp_top,
         scale = "row",                     # 按行标准化（Z-score）
         cluster_rows = TRUE,               # 对基因聚类
         cluster_cols = TRUE,               # 对细胞群聚类
         clustering_method = "ward.D2",     # 聚类方法（可选 "complete", "ward.D2" 等）
         treeheight_row = 40,               # 基因聚类树的高度（像素）
         treeheight_col = 30,               # 细胞群聚类树的高度
         color = my_colors,                 # 自定义颜色
         main = "Top 10 marker genes per cluster",   # 标题
         fontsize_row = 7,                  # 行标签字体大小
         fontsize_col = 10,                 # 列标签字体大小
         angle_col = 45,                    # 列名旋转角度（避免重叠）
         display_numbers = FALSE,           # 是否显示数值（数据量大时不推荐）
         number_color = "black",            # 数值颜色
         fontsize_number = 2,               # 数值字体大小
         border_color = "grey60",           # 单元格边框颜色
         legend = TRUE,                     # 显示图例
         legend_breaks = c(-2, 0, 2),       # 图例刻度（根据实际数据范围调整）
         legend_labels = c("Low", "Mean", "High"),
         filename = NA                      # 已在 pdf() 中指定，此处不重复
)

dev.off()
# 查看 Other 簇中肿瘤标志基因的平均表达
tumor_markers <- c("EGFR", "SOX2", "OLIG2", "PTPRZ1", "VIM", "GFAP")
present_markers <- intersect(tumor_markers, rownames(seurat_obj))
if (length(present_markers) > 0) {
  avg_tumor <- AverageExpression(seurat_obj, features = present_markers, assays = "RNA", slot = "data")$RNA
  print(avg_tumor)
}

# ------------------------------
# 导出每个细胞群的 top10 差异基因（带完整统计量）
# ------------------------------

# 确保 markers 已存在（如果之前已计算，可直接使用）
# markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# 按群组提取 top10（按 avg_log2FC 降序）
top10_full <- markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%   # 每个群内按 log2FC 排序
  slice_head(n = 10) %>%                            # 取前10
  ungroup() %>%
  select(cluster, gene, avg_log2FC, p_val, p_val_adj, pct.1, pct.2)  # 选择需要的列

# 显示前几行（预览）
print(head(top10_full, 20))

# 保存为 CSV 文件（兼容 Excel）
write.csv(top10_full, file.path(out_dir, "top10_markers_per_cluster_full.csv"), row.names = FALSE)

# 可选：同时保存为 Excel 格式（需要 openxlsx 包）
if (requireNamespace("openxlsx", quietly = TRUE)) {
  openxlsx::write.xlsx(top10_full, file = file.path(out_dir, "top10_markers_per_cluster_full.xlsx"))
} else {
  message("未安装 openxlsx 包，仅保存 CSV 文件。")
}

# 打印每个群的 top 基因列表（便于快速查阅）
cat("\n===== 各细胞群 Top10 差异基因列表 =====\n")
for (cl in unique(top10_full$cluster)) {
  cat("\n集群 ", cl, ":\n")
  sub <- top10_full[top10_full$cluster == cl, ]
  print(sub[, c("gene", "avg_log2FC", "p_val_adj", "pct.1")])
}

# 5. 基于热图和 top10 基因，手动判断每个簇的细胞类型
#   参考胶质瘤常见细胞类型（简化命名）：
#   - Tumor        : EGFR, SOX2, OLIG2, PTPRZ1, VIM, GFAP, TNC 等
#   - Macrophage   : CD68, C1QA, C1QB, CSF1R, CD163, MSR1
#   - T/NK         : CD3D, CD3E, CD8A, CD4, NKG7, KLRB1
#   - B cell       : MS4A1, CD79A
#   - Neutrophil   : FCGR3B, S100A8, S100A9, CXCR2
#   - Endothelial  : PECAM1, VWF, CLDN5
#   - Fibroblast   : COL1A1, DCN, ACTA2
#   - Oligodendrocytes : OLIG1, MOG, MBP, OPALIN
#   - Other        : 无法明确归类的簇

# 根据您之前的热图输出和 top10 基因，我们提出以下映射建议（可自行调整）：
# 注意：簇号是字符型，请根据实际簇数修改（当前为0-17）
# 手动映射（基于 top10 基因分析）
manual_mapping <- c(
  "0" = "Macrophage",
  "1" = "Oligodendrocytes",
  "2" = "Tumor",          # 高表达 PTPRZ1, VIM, SOX2
  "3" = "Macrophage",
  "4" = "Tumor",          # 极高 PTPRZ1, VIM
  "5" = "Fibroblast",
  "6" = "Tumor",          # 高 PTPRZ1, VIM
  "7" = "Tumor",          # 增殖性肿瘤
  "8" = "Neutrophil",
  "9" = "Tumor",          # 高 PTPRZ1, VIM, GFAP
  "10" = "Tumor",         # 高 PTPRZ1, VIM
  "11" = "Tumor",         # 角蛋白、CT83 等
  "12" = "Tumor",         # 高 PTPRZ1, VIM
  "13" = "T_NK",
  "14" = "Fibroblast",
  "15" = "Macrophage",
  "16" = "Macrophage",
  "17" = "B_cell"
)
# 移除旧列
seurat_obj$cell_type <- NULL

# 赋值
cell_type_vec <- manual_mapping[as.character(seurat_obj$seurat_clusters)]
names(cell_type_vec) <- NULL
seurat_obj$cell_type <- cell_type_vec

# 查看细胞类型分布
table(seurat_obj$cell_type)
cat("肿瘤细胞数量：", sum(seurat_obj$cell_type == "Tumor"), "\n")

# 绘制 UMAP 图
pdf(file.path(out_dir, "umap_cell_type_final.pdf"), width = 8, height = 6)
DimPlot(seurat_obj, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE) +
  ggtitle("UMAP by manual cell type (final)")
dev.off()

# 保存最终 Seurat 对象
saveRDS(seurat_obj, file = file.path(out_dir, "cgga_scRNA_processed_final.rds"))

cat("注释完成，对象已保存至：", file.path(out_dir, "cgga_scRNA_processed_final.rds"), "\n")
#------------------
# 9. 保存Seurat对象
# ------------------------------
cat("保存Seurat对象...\n")
saveRDS(seurat_obj, file = file.path(out_dir, "cgga_scRNA_processed.rds"))

# 同时保存UMAP坐标
umap_coords <- Embeddings(seurat_obj, reduction = "umap")
write.csv(umap_coords, file.path(out_dir, "umap_coordinates.csv"))

cat("分析完成！所有结果已保存至：", out_dir, "\n")

