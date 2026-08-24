# ==================== CellChat 分析（基于正确 cluster 编号） ====================
# 成纤维细胞 = cluster 5 和 14
# 高NETs中性粒细胞 = cluster 8 且 NETS_group == "High"

library(Seurat)
library(CellChat)
library(ggplot2)
library(patchwork)
library(dplyr)

# 路径设置（请按实际情况修改）
seurat_obj_path <- "E:/NETS/单细胞数据/结果/cgga_scRNA_processed_with_NETS_binary_rank.rds"
out_dir <- "E:/NETS/单细胞数据/结果/cell_communication"
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 读取对象
seurat_obj <- readRDS(seurat_obj_path)

# 确保有 seurat_clusters 列，并转为数值型
if(!"seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
  stop("未找到 seurat_clusters 列，请检查 Seurat 对象")
}
seurat_obj$cluster_num <- as.numeric(as.character(seurat_obj$seurat_clusters))

# 提取目标细胞：成纤维细胞 (5,14) + 高NETs中性粒细胞 (8)
target_cells <- WhichCells(seurat_obj,
                           expression = (cluster_num == 8 & NETS_group == "High") |
                             cluster_num %in% c(5, 14))
cat("目标细胞数量:", length(target_cells), "\n")
cat("  高NETs中性粒细胞 (cluster 8):", sum(seurat_obj$cluster_num[target_cells] == 8), "\n")
cat("  成纤维细胞 (cluster 5 & 14):", sum(seurat_obj$cluster_num[target_cells] %in% c(5,14)), "\n")

# 子集
subset_obj <- subset(seurat_obj, cells = target_cells)

# 创建通讯用的细胞类型标签（基于正确的 cluster 编号）
subset_obj$cell_type_comm <- case_when(
  subset_obj$cluster_num %in% c(5,14) ~ "Fibroblast",
  subset_obj$cluster_num == 8 ~ "Neutrophil_HighNETS",
  TRUE ~ "Other"
)
# 剔除可能的 Other（保险）
subset_obj <- subset(subset_obj, cell_type_comm != "Other")
subset_obj$cell_type_comm <- factor(subset_obj$cell_type_comm)
Idents(subset_obj) <- "cell_type_comm"

# 创建 CellChat 对象
data.input <- GetAssayData(subset_obj, assay = "RNA", layer = "data")
cellchat <- createCellChat(object = data.input,
                           meta = subset_obj@meta.data,
                           group.by = "cell_type_comm")

# 使用人类数据库
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB

# 预处理
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# 计算通讯概率
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 5)

# 通路水平推断
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# 可视化：整体网络
pdf(file.path(out_dir, "cellchat_network_aggregated.pdf"), width = 8, height = 6)
netVisual_circle(cellchat@net$count,
                 vertex.weight = as.numeric(table(cellchat@idents)),
                 weight.scale = TRUE,
                 label.edge = TRUE,
                 title.name = "Number of interactions")
dev.off()

# 可视化：主要信号通路（前4个）
pathways <- cellchat@netP$pathways
if(length(pathways) > 0) {
  pdf(file.path(out_dir, "cellchat_signaling_pathways.pdf"), width = 12, height = 10)
  for(i in 1:min(4, length(pathways))) {
    netVisual_aggregate(cellchat, signaling = pathways[i],
                        vertex.receiver = which(levels(cellchat@idents) == "Fibroblast"),
                        layout = "circle", vertex.label.cex = 1.2)
  }
  dev.off()
}

# 保存对象和结果
saveRDS(cellchat, file.path(out_dir, "cellchat_object.rds"))
write.csv(cellchat@net$prob, file.path(out_dir, "cellchat_ligand_receptor_prob.csv"))
write.csv(cellchat@net$pval, file.path(out_dir, "cellchat_ligand_receptor_pval.csv"))

# 提取从 Neutrophil_HighNETS 到 Fibroblast 的显著配体-受体对
lr_pairs <- subsetCommunication(cellchat, sources.use = "Neutrophil_HighNETS",
                                targets.use = "Fibroblast")
write.csv(lr_pairs, file.path(out_dir, "cellchat_LR_Neutrophil_to_Fibroblast.csv"), row.names = FALSE)

# 反向：Fibroblast 到 Neutrophil
lr_reverse <- subsetCommunication(cellchat, sources.use = "Fibroblast",
                                  targets.use = "Neutrophil_HighNETS")
write.csv(lr_reverse, file.path(out_dir, "cellchat_LR_Fibroblast_to_Neutrophil.csv"), row.names = FALSE)

cat("CellChat 分析完成！结果保存在：", out_dir, "\n")

