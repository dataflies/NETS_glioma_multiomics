# ==================== 拟时序分析：高NETs中性粒细胞 vs 成纤维细胞 ====================
# 使用Slingshot推断轨迹 + tradeSeq识别随伪时间变化的基因
# 输入：已处理好的Seurat对象（包含细胞类型和NETS_group注释）
# 输出：轨迹图、伪时间分布、随伪时间变化的差异基因列表

# ------------------------------
# 0. 路径设置（请根据实际情况修改）
# ------------------------------
seurat_obj_path <- "E:/NETS/单细胞数据/结果/cgga_scRNA_processed_with_NETS_binary_rank.rds"  # 最终包含NETS_group的Seurat对象
out_dir <- "E:/NETS/单细胞数据/结果/trajectory_analysis"
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ------------------------------
# 1. 加载所需包（如缺失请先安装）
# ------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(slingshot)   # 轨迹推断
library(tradeSeq)    # 伪时间差异表达
library(RColorBrewer)
library(patchwork)

# ------------------------------
# 2. 读取Seurat对象，提取目标细胞
# ------------------------------
cat("读取Seurat对象...\n")
seurat_obj <- readRDS(seurat_obj_path)

# 检查必要列是否存在
if(!"cell_type" %in% colnames(seurat_obj@meta.data)) stop("未找到 cell_type 列")
if(!"NETS_group" %in% colnames(seurat_obj@meta.data)) stop("未找到 NETS_group 列")

# 提取目标细胞：成纤维细胞 + 高NETs中性粒细胞
target_cells <- WhichCells(seurat_obj, 
                           expression = (seurat_clusters %in% c(5, 14)) | 
                             (seurat_clusters == 8 & NETS_group == "High"))
cat(paste0("目标细胞数量：", length(target_cells), "\n"))
cat("  成纤维细胞：", sum(seurat_obj$cell_type[target_cells] == "Fibroblast"), "\n")
cat("  高NETs中性粒细胞：", sum(seurat_obj$cell_type[target_cells] == "Neutrophil" & seurat_obj$NETS_group[target_cells] == "High"), "\n")

if(length(target_cells) < 20) {
  stop("目标细胞数量过少（<20），不适合进行拟时序分析")
}

# 创建子集对象
subset_obj <- subset(seurat_obj, cells = target_cells)

# ------------------------------
# 3. 数据准备：提取UMAP降维坐标（Slingshot的输入）
# ------------------------------
# 确保对象中有UMAP坐标
if(!"umap" %in% names(subset_obj@reductions)) stop("Seurat对象中找不到UMAP降维结果")
umap_emb <- Embeddings(subset_obj, "umap")
colnames(umap_emb) <- c("UMAP1", "UMAP2")

# 提取聚类标签（这里直接用cell_type作为分组，Slingshot会在两组之间拟合曲线）
# 创建正确的细胞标签（用于Slingshot）
subset_obj$trajectory_group <- ifelse(subset_obj$seurat_clusters %in% c(5,14), "Fibroblast", "Neutrophil_HighNETS")
cell_labels <- subset_obj$trajectory_group
# 可选：将中性粒细胞高NETs组单独命名，便于区分
cell_labels[cell_labels == "Neutrophil"] <- "Neutrophil_HighNETS"
table(cell_labels)

# ------------------------------
# 4. 运行Slingshot推断轨迹
# ------------------------------
# 设置随机种子以保证结果可重复
set.seed(123)

# 使用Slingshot，指定起始簇（起始点为中性粒细胞，假设轨迹从高NETs中性粒细胞指向成纤维细胞）
# 若无先验知识，也可让Slingshot自动推断，但指定起始簇有助于解释
start_cluster <- "Neutrophil_HighNETS"
sds <- slingshot(umap_emb, 
                 clusterLabels = cell_labels,
                 start.clus = start_cluster,
                 reducedDim = 'UMAP')

cat("\n轨迹推断完成。\n")
cat("拟合的曲线数量：", length(slingCurves(sds)), "\n")

# ------------------------------
# 5. 可视化轨迹：UMAP图上叠加曲线，颜色为伪时间
# ------------------------------
# 提取伪时间（每个细胞沿最近曲线的伪时间）
pseudotime <- slingPseudotime(sds)
# 对于每个细胞，取存在的伪时间（非NA）的第一个曲线（通常只有一条曲线拟合）
pseudotime_primary <- apply(pseudotime, 1, function(x) {
  val <- x[!is.na(x)]
  if(length(val) > 0) val[1] else NA
})

# 将伪时间添加到meta.data中
subset_obj$pseudotime <- pseudotime_primary

# 绘制UMAP轨迹图
p1 <- ggplot(as.data.frame(umap_emb), aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = pseudotime_primary), size = 1, alpha = 0.8) +
  scale_color_gradient(low = "blue", high = "red", na.value = "grey80") +
  labs(title = "Trajectory inferred by Slingshot", color = "Pseudotime") +
  theme_minimal()

# 添加Slingshot的曲线
curve_df <- slingCurves(sds)[[1]]$s  # 第一条曲线
p1 <- p1 + geom_path(data = as.data.frame(curve_df), aes(x = UMAP1, y = UMAP2), 
                     linewidth = 1.2, linetype = "dashed", color = "black")

# 保存
ggsave(file.path(out_dir, "trajectory_umap_pseudotime.pdf"), p1, width = 7, height = 6)

# 也可以按细胞类型着色
p2 <- ggplot(as.data.frame(umap_emb), aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = cell_labels), size = 1, alpha = 0.8) +
  scale_color_manual(values = c("Neutrophil_HighNETS" = "red", "Fibroblast" = "steelblue")) +
  geom_path(data = as.data.frame(curve_df), aes(x = UMAP1, y = UMAP2), 
            linewidth = 1.2, linetype = "dashed", color = "black") +
  labs(title = "Trajectory by cell type", color = "Group") +
  theme_minimal()
ggsave(file.path(out_dir, "trajectory_umap_celltype.pdf"), p2, width = 7, height = 6)

# 伪时间分布箱线图
p3 <- ggplot(subset_obj@meta.data, aes(x = cell_labels, y = pseudotime, fill = cell_labels)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  scale_fill_manual(values = c("Neutrophil_HighNETS" = "red", "Fibroblast" = "steelblue")) +
  labs(title = "Pseudotime distribution", x = "", y = "Pseudotime") +
  theme_minimal() +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "pseudotime_distribution.pdf"), p3, width = 5, height = 4)

# ------------------------------
# 6. 使用tradeSeq识别随伪时间变化的基因（沿轨迹的表达模式）
# ------------------------------
# 准备表达矩阵（counts 或 log-normalized count）
# tradeSeq推荐使用log-counts，这里使用Seurat对象的RNA数据（已LogNormalize）
expr_matrix <- as.matrix(GetAssayData(subset_obj, assay = "RNA", layer = "data"))  # log1p normalized

# 确保细胞顺序与sds对象一致
cells_in_sds <- rownames(umap_emb)  # Slingshot使用的细胞
expr_matrix <- expr_matrix[, cells_in_sds]

# 构建tradeSeq的拟合对象（负二项式模型，适用于count数据；但对于log数据，使用gaussian可能更合适）
# 更稳健：使用count数据 + offset，这里简化：因为数据已经是log标准化，使用tradeSeq的"gaussian"家族
# 注意：tradeSeq需对每个基因拟合GAM，可能较慢，可先筛选高变基因或在一定表达水平的基因。

# 筛选表达基因：至少在10%细胞中表达>0（log标准化后>0对应原始count>1）
expressed_genes <- names(which(rowSums(expr_matrix > 0) > 0.1 * ncol(expr_matrix)))
cat(paste0("共有 ", length(expressed_genes), " 个基因表达达标\n"))

# 从子集对象中提取高变基因（如果尚未计算，先计算）
if(length(VariableFeatures(subset_obj)) == 0) {
  cat("正在计算高变基因...\n")
  subset_obj <- FindVariableFeatures(subset_obj, selection.method = "vst", nfeatures = 2000)
}
hv_genes <- VariableFeatures(subset_obj)

# 取高变基因与表达基因的交集
genes_to_fit <- intersect(hv_genes, expressed_genes)
cat(paste0("高变基因数量（与表达基因交集）：", length(genes_to_fit), "\n"))

# 如果高变基因太少，则退回到所有表达基因的前2000个（按现有顺序）
if(length(genes_to_fit) < 100) {
  cat("警告：高变基因数量过少，将使用所有表达基因的前2000个\n")
  genes_to_fit <- expressed_genes[1:min(2000, length(expressed_genes))]
}

# 使用fitGAM函数拟合广义加性模型
set.seed(123)
cat("开始拟合GAM模型（基因数量：", length(genes_to_fit), "，预计耗时数分钟至数十分钟）...\n")

gam_list <- fitGAM(sds, 
                   counts = expr_matrix[genes_to_fit, ], 
                   family = "gaussian",  # 因为expr_matrix已经是log标准化数据
                   nknots = 6, 
                   verbose = TRUE)

# 识别随伪时间显著变化的基因（关联检验）
asso_res <- associationTest(gam_list)
asso_res$gene <- rownames(asso_res)
asso_sig <- asso_res[asso_res$pvalue < 0.05, ]
asso_sig <- asso_sig[order(asso_sig$pvalue), ]
cat(paste0("与伪时间显著相关的基因数量（p<0.05）：", nrow(asso_sig), "\n"))

# 保存结果
write.csv(asso_res, file.path(out_dir, "tradeSeq_association_test_all.csv"), row.names = FALSE)
if(nrow(asso_sig) > 0) {
  write.csv(asso_sig, file.path(out_dir, "tradeSeq_significant_genes.csv"), row.names = FALSE)
}

# 可选：绘制top显著基因随伪时间变化的表达趋势
if(nrow(asso_sig) > 0) {
  top_genes <- head(asso_sig$Gene, 6)
  pdf(file.path(out_dir, "top_genes_pseudotime_expression.pdf"), width = 8, height = 6)
  for(gene in top_genes) {
    p <- plotSmoothers(gam_list, Gene) +
      ggtitle(paste0(gene, " (p = ", format(asso_sig[gene, "pvalue"], digits = 3), ")"))
    print(p)
  }
  dev.off()
}

# ------------------------------
# 7. 保存子集Seurat对象（包含伪时间）
# ------------------------------
saveRDS(subset_obj, file = file.path(out_dir, "trajectory_subset_NeutrophilHigh_Fibroblast.rds"))
write.csv(subset_obj@meta.data, file = file.path(out_dir, "trajectory_cell_metadata.csv"), row.names = TRUE)

cat("\n拟时序分析完成！所有结果保存在：", out_dir, "\n")




# 提取伪时间（每个细胞一个值，取第一列非NA）
pseudotime <- slingPseudotime(sds)[,1]  # 第一曲线
names(pseudotime) <- rownames(slingPseudotime(sds))

# 获取表达矩阵（log标准化数据）
expr_matrix <- as.matrix(GetAssayData(subset_obj, assay = "RNA", layer = "data"))

# 创建数据框用于绘图
plot_df <- data.frame(
  pseudotime = pseudotime,
  S100A8 = expr_matrix["S100A8", names(pseudotime)],
  FN1 = expr_matrix["FN1", names(pseudotime)],
  cell_type = subset_obj$cell_type[names(pseudotime)]
)

# 绘制散点图 + 平滑曲线
library(ggplot2)
p1 <- ggplot(plot_df, aes(x = pseudotime, y = S100A8)) +
  geom_point(aes(color = cell_type), alpha = 0.6) +
  geom_smooth(method = "loess", se = TRUE, color = "red") +
  labs(title = "S100A8 expression along pseudotime", y = "Log-normalized expression") +
  theme_minimal()

p2 <- ggplot(plot_df, aes(x = pseudotime, y = FN1)) +
  geom_point(aes(color = cell_type), alpha = 0.6) +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  labs(title = "FN1 expression along pseudotime", y = "Log-normalized expression") +
  theme_minimal()

# 合并显示
library(patchwork)
p1 + p2

# 保存为 PDF 文件
pdf(file.path(out_dir, "pseudotime_validation_S100A8_FN1.pdf"), width = 10, height = 5)
p1 + p2
dev.off()

# 同时保存为 PNG（便于快速查看）
png(file.path(out_dir, "pseudotime_validation_S100A8_FN1.png"), width = 1000, height = 500, res = 120)
p1 + p2
dev.off()

cat("图片已保存至：", file.path(out_dir, "pseudotime_validation_S100A8_FN1.pdf"), "\n")
# 查看起始簇是否确实是 Neutrophil_HighNETS
startNode(sds)
# 或者查看曲线方向（第一条曲线的起始细胞坐标）
curve <- slingCurves(sds)[[1]]
head(curve$s)  # 查看曲线起始点坐标
# 方法1：查看 sds 对象的 metadata
sds@metadata$startClust

# 方法2：查看 slingParams
slingParams(sds)$start.clus



# ==================== 中性粒细胞按伪时间三分位比较（高 vs 低） ====================
# 前提：subset_obj 已经运行过 slingshot，且包含 pseudotime 列
# 中性粒细胞为 cluster 8

library(Seurat)
library(ggplot2)
library(ggrepel)
library(dplyr)

# ------------------------------
# 1. 提取中性粒细胞子集（基于正确的 cluster 编号 8）
# ------------------------------
if(!"cluster_num" %in% colnames(subset_obj@meta.data)) {
  subset_obj$cluster_num <- as.numeric(as.character(subset_obj$seurat_clusters))
}

neutro_cells <- WhichCells(subset_obj, expression = cluster_num == 8)
neutro_obj <- subset(subset_obj, cells = neutro_cells)

# 检查 pseudotime 是否存在
if(!"pseudotime" %in% colnames(neutro_obj@meta.data)) {
  stop("neutro_obj 中缺少 pseudotime 列，请先运行 slingshot 并赋值 pseudotime。")
}

# ------------------------------
# 2. 根据伪时间三分位数分组（低、中、高）
# ------------------------------
pseudo_vals <- neutro_obj$pseudotime
q33 <- quantile(pseudo_vals, 1/3, na.rm = TRUE)
q67 <- quantile(pseudo_vals, 2/3, na.rm = TRUE)

low_cells  <- names(pseudo_vals[pseudo_vals <= q33])
mid_cells  <- names(pseudo_vals[pseudo_vals > q33 & pseudo_vals < q67])
high_cells <- names(pseudo_vals[pseudo_vals >= q67])

cat("低伪时间组细胞数:", length(low_cells), "\n")
cat("中伪时间组细胞数:", length(mid_cells), "\n")
cat("高伪时间组细胞数:", length(high_cells), "\n")

# 添加分组标签
neutro_obj$pseudotime_group <- "Middle"
neutro_obj$pseudotime_group[low_cells]  <- "Low"
neutro_obj$pseudotime_group[high_cells] <- "High"
Idents(neutro_obj) <- "pseudotime_group"

# ------------------------------
# 3. 差异表达分析：High vs Low
# ------------------------------
deg <- FindMarkers(neutro_obj, 
                   ident.1 = "High", 
                   ident.2 = "Low",
                   logfc.threshold = 0.25, 
                   min.pct = 0.1, 
                   verbose = TRUE)
deg$gene <- rownames(deg)
deg <- deg[order(deg$avg_log2FC, decreasing = TRUE), ]

# 保存结果
out_dir <- "E:/NETS/单细胞数据/结果/trajectory_analysis"
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
write.csv(deg, file.path(out_dir, "neutro_pseudotime_High_vs_Low.csv"), row.names = FALSE)

# 显示前20个上调基因
cat("\nHigh vs Low 上调基因 TOP 20:\n")
print(head(deg[deg$avg_log2FC > 0, ], 20))

# ------------------------------
# 4. 火山图（可自定义阈值和标记基因）
# ------------------------------
# 处理 p 值为0的情况
if(any(deg$p_val == 0, na.rm = TRUE)) {
  min_nonzero <- min(deg$p_val[deg$p_val > 0], na.rm = TRUE)
  deg$p_val[deg$p_val == 0] <- min_nonzero
}
deg$log_p <- -log10(deg$p_val)
deg <- deg[is.finite(deg$log_p), ]

# 设定阈值（可根据实际调整）
log2FC_threshold <- 1.0   # 建议先用1.0，如果上下调基因太少可降低
p_threshold <- 0.05

# 标记您感兴趣的基因（示例，可自行修改）
ligands<- c("HBEGF", "PDGFA", "OSM", "IL16", "SPP1", "NAMPT", "GRN", "LGALS9", "KLK6", 
                         "PTPRS", "FN1", "THBS2", "COL6A1", "ENTPD1", "CD46", "CNTN2", "JAM2", 
                 "MAG", "HLA-DPB1", 
                         "HLA-DMA", "HLA-DMB", "HLA-DRA", "HLA-DRB1", "SELL", "SEMA4D", "CD55")
receptors <- c("TGFBR2", "CXCR1", "CXCR2", "IL1R2", "CSF1R", "ITGA5_ITGB1", "TLR4", "LILRB3",
               "C3AR1", "ITGAM_ITGB2", "ITGAX_ITGB2", "FPR2", "FPR1", "AXL", "MERTK", "CD47", 
               "PTGER4", "CD74", "TMIGD3",
               "PILRA", "CNTN2", "CD4", "PTPRF", "ADGRE5", "PTPRS", "SORL1", "TREM2_TYROBP")
label_genes <- unique(c(ligands, receptors))
# 只保留在 deg 中存在的基因
label_genes <- intersect(label_genes, deg$gene)

deg$significance <- "Not significant"
deg$significance[abs(deg$avg_log2FC) >= log2FC_threshold & deg$p_val < p_threshold] <- "Significant"
deg$significance[deg$avg_log2FC >= log2FC_threshold & deg$p_val < p_threshold] <- "Up"
deg$significance[deg$avg_log2FC <= -log2FC_threshold & deg$p_val < p_threshold] <- "Down"

deg$label <- ifelse(deg$gene %in% label_genes, deg$gene, NA)

# 绘制火山图
p_volcano <- ggplot(deg, aes(x = avg_log2FC, y = log_p, color = significance, size = significance)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("Not significant" = "grey80", 
                                "Significant" = "grey50",
                                "Up" = "red", 
                                "Down" = "blue")) +
  scale_size_manual(values = c("Not significant" = 0.8, 
                               "Significant" = 1,
                               "Up" = 1.2, 
                               "Down" = 1.2)) +
  geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), 
             linetype = "dashed", color = "black", alpha = 0.7) +
  geom_hline(yintercept = -log10(p_threshold), 
             linetype = "dashed", color = "black", alpha = 0.7) +
  geom_text_repel(data = subset(deg, !is.na(label)), 
                  aes(label = label),
                  color = "black", size = 4, box.padding = 0.5, 
                  point.padding = 0.3, force = 10, max.overlaps = 20) +
  labs(x = expression(Log[2]~Fold~Change~(High~vs~Low)),
       y = expression(-Log[10]~P~value),
       title = "Neutrophil: High pseudotime vs Low pseudotime") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  coord_cartesian(ylim = c(0, max(deg$log_p, na.rm = TRUE) * 0.95))

# 保存
ggsave(file.path(out_dir, "neutro_pseudotime_High_vs_Low_volcano.pdf"), 
       p_volcano, width = 8, height = 7)

cat("\n分析完成！结果保存在：", out_dir, "\n")

# 确保 deg 数据框包含 gene, avg_log2FC, p_val（或 p_val_adj）
# 如果使用了 p_val_adj，请将下面的 p_val 改为 p_val_adj

# 定义阈值（与火山图一致）
log2FC_threshold <- 1.0
p_threshold <- 0.05

# 只考虑标记基因中存在于 deg 的
label_genes_present <- intersect(label_genes, deg$gene)

# 创建分类函数
classify_gene <- function(gene) {
  fc <- deg$avg_log2FC[deg$gene == gene]
  p <- deg$p_val[deg$gene == gene]
  if (is.na(fc) | is.na(p)) return(NA)
  
  if (fc <= -log2FC_threshold & p < p_threshold) return("Top-Left (sig down)")
  if (fc >= log2FC_threshold & p < p_threshold) return("Top-Right (sig up)")
  if (fc <= -log2FC_threshold & p >= p_threshold) return("Bottom-Left (trend down, ns)")
  if (fc >= log2FC_threshold & p >= p_threshold) return("Bottom-Right (trend up, ns)")
  if (abs(fc) < log2FC_threshold & p < p_threshold) return("Top-Middle (small FC, sig)")
  if (abs(fc) < log2FC_threshold & p >= p_threshold) return("Bottom-Middle (no diff)")
  return("Other")
}

# 对每个标记基因分类
gene_class <- sapply(label_genes_present, classify_gene)

# 按分区整理列表
classified_list <- split(label_genes_present, gene_class)

# 打印结果（使用 cat 输出）
cat("\n========== 火山图标记基因分区结果 ==========\n\n")
for (region in names(classified_list)) {
  cat(sprintf("【%s】 (%d 个基因):\n", region, length(classified_list[[region]])))
  cat(paste(classified_list[[region]], collapse = ", "), "\n\n")
}

# 可选：输出未在 deg 中找到的标记基因
missing <- setdiff(label_genes, label_genes_present)
if (length(missing) > 0) {
  cat("⚠️ 以下标记基因不在 deg 数据中（无法分类）:\n")
  cat(paste(missing, collapse = ", "), "\n")
}



# ==================== GO BP 富集分析（浅蓝 vs 深蓝 上调基因） ====================
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)

# 1. 筛选上调基因（LightBlue 高表达）
# 根据您的数据，top 基因 log2FC 极高，阈值可放宽
up_genes <- deg %>%
  filter(avg_log2FC > 1, p_val_adj < 0.05) %>%
  pull(gene)

cat(paste("显著上调基因数量:", length(up_genes), "\n"))

# 2. ID 转换（Symbol -> Entrez ID）
entrez_ids <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

# 3. GO 富集分析
ego <- enrichGO(gene = entrez_ids$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont = "BP",          # 生物学过程
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.10,
                readable = TRUE)

# 4. 输出结果
if(!is.null(ego) && nrow(ego) > 0) {
  write.csv(as.data.frame(ego), file.path(out_dir, "GO_BP_upregulated_LightBlue.csv"), row.names = FALSE)
  
  # 绘制 dotplot
  pdf(file.path(out_dir, "GO_BP_dotplot.pdf"), width = 10, height = 8)
  print(dotplot(ego, showCategory = 15, title = "GO BP Enrichment (Up in LightBlue vs DeepBlue)"))
  dev.off()
  
  # 绘制 barplot
  pdf(file.path(out_dir, "GO_BP_barplot.pdf"), width = 10, height = 6)
  print(barplot(ego, showCategory = 15, title = "GO BP Enrichment"))
  dev.off()
  
  cat("GO BP 富集分析完成，结果保存至：", out_dir, "\n")
} else {
  cat("未发现显著富集的 GO BP 条目，可以尝试降低 log2FC 阈值（例如 >0.5）\n")
}


