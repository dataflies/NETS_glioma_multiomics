# ------------------------------
# 第二阶段：筛选高 NETS 活性的中性粒细胞群（二分法，按中位数分组）
# ------------------------------

# ==================== 路径设置 ====================
raw_data_dir <- "E:/NETS/单细胞数据/原数据"                # 存放 GO 基因集文件
out_dir <- "E:/NETS/单细胞数据/结果"                       # 结果输出文件夹
# 请确认这里使用的是最终注释好的 Seurat 对象（包含正确的 cell_type 列）
rds_file <- file.path(out_dir, "cgga_scRNA_processed_final.rds")  # 根据实际文件名修改
# ==================================================

# 加载必要的包
library(Seurat)
library(dplyr)
library(ggplot2)
library(UCell)       # 需要先安装：BiocManager::install("UCell")
library(ggpubr)      # 用于箱线图添加统计检验

# 1. 加载 Seurat 对象
if (!file.exists(rds_file)) stop("Seurat 对象文件不存在：", rds_file)
seurat_obj <- readRDS(rds_file)
cat("Seurat 对象加载成功，细胞数：", ncol(seurat_obj), "\n")

# 2. 读取本地 GO 基因集文件（提取第一列基因名）
go_files <- c("GO0140644.txt", "GO0140645.txt")
go_genes_list <- list()
for (f in go_files) {
  file_path <- file.path(raw_data_dir, f)
  if (!file.exists(file_path)) {
    stop("文件不存在：", file_path)
  }
  go_data <- read.table(file_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  genes <- unique(go_data[[1]])   # 第一列是基因符号
  go_genes_list[[f]] <- genes
  cat(f, "中包含", length(genes), "个唯一基因\n")
}

# 合并两个 GO 集的基因
go_genes <- unique(unlist(go_genes_list))
cat("合并后 GO 基因总数：", length(go_genes), "\n")

# 3. 添加文献中的三个关键基因
key_genes <- c("PTGS2", "MME", "SLC2A3")
nets_genes <- unique(c(go_genes, key_genes))
cat("最终 NETS 基因总数：", length(nets_genes), "\n")

# 4. 提取中性粒细胞子集（用于后续分析）
if (!"cell_type" %in% colnames(seurat_obj@meta.data)) {
  stop("seurat_obj 中未找到 cell_type 列，请先完成细胞类型注释")
}
neutro_cells <- WhichCells(seurat_obj, expression = cell_type == "Neutrophil")
if (length(neutro_cells) == 0) {
  stop("没有找到标记为 Neutrophil 的细胞，请检查 cell_type 列中的命名")
}
neutro_obj <- subset(seurat_obj, cells = neutro_cells)
cat("中性粒细胞数量：", ncol(neutro_obj), "\n")

# 5. 使用 UCell 计算 NETS 活性评分（在整个 Seurat 对象上进行，方便绘图）
nets_genes_present <- intersect(nets_genes, rownames(seurat_obj))
cat("数据中存在的 NETS 基因数：", length(nets_genes_present), "\n")
if (length(nets_genes_present) < 5) {
  warning("数据中存在的 NETS 基因较少，评分可能不稳健")
}

gene_set <- list(NETS = nets_genes_present)
seurat_obj <- AddModuleScore_UCell(seurat_obj, features = gene_set, name = "NETS_UCell")

# 查看评分列名
score_col <- grep("NETS_UCell", colnames(seurat_obj@meta.data), value = TRUE)[1]
cat("UCell 评分列名：", score_col, "\n")
summary(seurat_obj[[score_col]])


# 6. 二分法：按排序均分（高、低各一半）
# 计算中性粒细胞子集的评分
score_vals <- seurat_obj[[score_col]][neutro_cells, 1]

# 按评分排序，取中位索引，将细胞分为两组数量相等的组（如果奇数，高组多1）
rank_idx <- order(score_vals)
n_neutro <- length(neutro_cells)
cut_point <- ceiling(n_neutro / 2)  # 高组个数（向上取整）

# 低组：评分较小的前一半；高组：后一半
low_cells <- neutro_cells[rank_idx[1:cut_point]]
high_cells <- neutro_cells[rank_idx[(cut_point+1):n_neutro]]

seurat_obj$NETS_group <- "Medium"  # 临时
seurat_obj$NETS_group[high_cells] <- "High"
seurat_obj$NETS_group[low_cells] <- "Low"

# 将非中性粒细胞设置为 "Other"（用于绘图）
seurat_obj$NETS_group[setdiff(colnames(seurat_obj), neutro_cells)] <- "Other"

# 查看中性粒细胞分组情况
cat("中性粒细胞分组统计（按排序均分）：\n")
print(table(seurat_obj$NETS_group[neutro_cells]))

# 7. 箱线图：比较高低组评分差异（使用中性粒细胞子集）
neutro_obj$NETS_group <- seurat_obj$NETS_group[neutro_cells]
p_box <- ggplot(neutro_obj@meta.data, aes(x = NETS_group, y = .data[[score_col]], fill = NETS_group)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  labs(title = "NETS activity (binary split by rank)", x = "Group", y = "UCell Score") +
  theme_minimal() +
  scale_fill_manual(values = c("High" = "red", "Low" = "blue"))

pdf(file.path(out_dir, "NETS_boxplot_binary_rank.pdf"), width = 5, height = 4)
print(p_box)
dev.off()

# 8. 整体 UMAP 高亮显示中性粒细胞（点大小调小）
pdf(file.path(out_dir, "umap_NETS_binary_rank_overall.pdf"), width = 8, height = 6)
DimPlot(seurat_obj, reduction = "umap", group.by = "NETS_group",
        cols = c("High" = "red", "Low" = "blue", "Other" = "grey90"),
        order = c("High", "Low", "Other"),
        pt.size = 0.2) +
  ggtitle("Neutrophil NETS binary groups (rank split)")
dev.off()

# 仅高亮中性粒细胞
pdf(file.path(out_dir, "umap_highlight_neutrophils_binary_rank.pdf"), width = 6, height = 5)
DimPlot(seurat_obj, reduction = "umap",
        cells.highlight = list(
          "High" = high_cells,
          "Low" = low_cells
        ),
        cols.highlight = c("red", "blue"),
        cols = "grey90",
        pt.size = 0.2) +
  ggtitle("Neutrophil binary groups (highlighted)")
dev.off()

# 9. 保存中性粒细胞子集及分组信息
neutro_obj <- subset(seurat_obj, cells = neutro_cells)
saveRDS(neutro_obj, file = file.path(out_dir, "neutrophils_with_NETS_binary_rank.rds"))
write.csv(neutro_obj@meta.data, file.path(out_dir, "neutrophil_metadata_binary_rank.csv"), row.names = TRUE)

# 10. 保存完整对象
saveRDS(seurat_obj, file = file.path(out_dir, "cgga_scRNA_processed_with_NETS_binary_rank.rds"))

cat("第二阶段完成！\n")
cat("中性粒细胞分组（按排序均分）：高活性", sum(neutro_obj$NETS_group == "High"),
    "，低活性", sum(neutro_obj$NETS_group == "Low"), "\n")

