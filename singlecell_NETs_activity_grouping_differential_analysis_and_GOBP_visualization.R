# ==================== 第三阶段：中性粒细胞子集深入分析 ====================
# 前提：已加载必要的包，已有 seurat_obj 和 neutro_obj 对象，且 neutro_obj 包含 NETS_group 列。

# 1. 单独的中性粒细胞群 UMAP 展示 NETS 高低情况
# 如果中性粒细胞子集还没有 UMAP 坐标，我们使用原始降维结果（整体 UMAP）并高亮显示
# 或者重新对中性粒细胞子集运行 UMAP（更专注）
library(ggplot2)
library(Seurat)

# 方法一：直接在整体 UMAP 上高亮中性粒细胞（已做过）
# 方法二：单独对中性粒细胞子集重新计算 UMAP（推荐，因为可以看到子群结构）
cat("对中性粒细胞子集重新运行 PCA 和 UMAP...\n")
neutro_obj <- RunPCA(neutro_obj, npcs = 30, verbose = FALSE)
neutro_obj <- RunUMAP(neutro_obj, dims = 1:20, reduction = "pca", seed.use = 42)

# 绘制 UMAP，按 NETS_group 着色
pdf(file.path(out_dir, "neutrophil_umap_NETS_groups.pdf"), width = 6, height = 5)
DimPlot(neutro_obj, reduction = "umap", group.by = "NETS_group",
        cols = c("High" = "red", "Low" = "blue"), pt.size = 0.5) +
  ggtitle("Neutrophil subpopulations by NETS activity")
dev.off()

# 2. 两个高低评分的箱线图，带统计检验（已经做过，但可以重新美化或保存）
# 确保 neutro_obj 有评分列，我们使用之前计算的 UCell 评分
score_col <- grep("NETS_UCell", colnames(neutro_obj@meta.data), value = TRUE)[1]
p_box <- ggplot(neutro_obj@meta.data, aes(x = NETS_group, y = .data[[score_col]], fill = NETS_group)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y = max(neutro_obj[[score_col]]) * 0.9) +
  labs(title = "NETS activity (rank split)", x = "Group", y = "UCell Score") +
  theme_minimal() +
  scale_fill_manual(values = c("High" = "red", "Low" = "blue"))

pdf(file.path(out_dir, "NETS_boxplot_binary_rank_v2.pdf"), width = 5, height = 4)
print(p_box)
dev.off()

# 3. 高、低组分别跑 GO 富集分析（基于差异表达基因）
# 首先提取高组和低组的中性粒细胞
high_cells <- WhichCells(neutro_obj, expression = NETS_group == "High")
low_cells <- WhichCells(neutro_obj, expression = NETS_group == "Low")

# 计算高组相对于低组的差异表达基因（上调基因）
Idents(neutro_obj) <- "NETS_group"
markers <- FindMarkers(neutro_obj, ident.1 = "High", ident.2 = "Low",
                       logfc.threshold = 0.25, min.pct = 0.1, only.pos = TRUE)
# 筛选显著差异基因（p_val_adj < 0.05）
sig_up <- rownames(markers[markers$p_val_adj < 0.05, ])
cat("高组相对于低组显著上调的基因数：", length(sig_up), "\n")

# 如果希望也看下调基因，可以设置 only.pos = FALSE，这里我们只取上调基因进行 GO 富集

# 如果没有显著基因，可以放宽阈值，比如 logfc.threshold = 0.1
if (length(sig_up) == 0) {
  warning("没有显著上调基因，尝试放宽阈值...")
  markers <- FindMarkers(neutro_obj, ident.1 = "High", ident.2 = "Low",
                         logfc.threshold = 0.1, min.pct = 0.1, only.pos = TRUE)
  sig_up <- rownames(markers[markers$p_val_adj < 0.05, ])
  cat("放宽后显著上调的基因数：", length(sig_up), "\n")
}

# 准备 GO 富集分析（需要安装 clusterProfiler 和 org.Hs.eg.db）
if (!require("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler")
  library(clusterProfiler)
}
if (!require("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
  library(org.Hs.eg.db)
}

# 将基因符号转换为 Entrez ID
if (length(sig_up) > 0) {
  gene_entrez <- bitr(sig_up, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  # 去除重复和NA
  gene_entrez <- unique(gene_entrez$ENTREZID)
  
  # 背景基因（所有在 neutro_obj 中表达的基因）
  all_genes <- rownames(neutro_obj)
  all_entrez <- bitr(all_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  all_entrez <- unique(all_entrez$ENTREZID)
  
  # GO 富集分析（BP, CC, MF 均可，此处以 BP 为例）
  ego <- enrichGO(gene = gene_entrez,
                  universe = all_entrez,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)
  
  # 保存结果
  if (!is.null(ego) && nrow(ego) > 0) {
    write.csv(as.data.frame(ego), file.path(out_dir, "GO_BP_high_vs_low.csv"), row.names = FALSE)
    # 绘制气泡图
    pdf(file.path(out_dir, "GO_BP_high_vs_low.pdf"), width = 8, height = 6)
    print(dotplot(ego, showCategory = 15, title = "GO BP enrichment (High vs Low)"))
    dev.off()
    cat("GO 富集分析完成，结果已保存。\n")
  } else {
    cat("未发现显著富集的 GO 条目。\n")
  }
} else {
  cat("没有显著差异基因，无法进行 GO 富集分析。\n")
}

# 可选：保存更新后的中性粒细胞子集（含UMAP坐标）
saveRDS(neutro_obj, file = file.path(out_dir, "neutrophils_with_NETS_umap.rds"))

cat("第三阶段完成！结果已保存至：", out_dir, "\n")
