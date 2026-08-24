# ==================== 0. 加载必要的 R 包 ====================
library(Seurat)          # 单细胞数据分析
library(CellChat)        # 细胞间通讯分析
library(dplyr)           # 数据处理
library(ggplot2)         # 绘图
library(clusterProfiler) # 富集分析（本例未使用，但保留）
library(org.Hs.eg.db)    # 人类基因注释（配合 clusterProfiler）

# ==================== 1. 读取数据 ====================
out_dir <- "E:/NETS/单细胞数据/结果"               # 结果输出目录
seurat_obj <- readRDS(file.path(out_dir, "cgga_scRNA_processed_with_NETS_binary_rank.rds"))

# ==================== 2. 提取中性粒细胞和肿瘤细胞 ====================
# 根据细胞类型注释筛选
seurat_sub <- subset(seurat_obj, subset = cell_type %in% c("Neutrophil", "Tumor"))

# 获取两类细胞的细胞名称
tumor_cells <- WhichCells(seurat_sub, expression = cell_type == "Tumor")
neut_cells  <- WhichCells(seurat_sub, expression = cell_type == "Neutrophil")

# ==================== 3. 计算每个细胞的 NETs 评分 ====================
# 读取 NETs 相关基因列表（由之前分析得到）
nets_genes <- scan(file.path(out_dir, "tumor_genes_up_in_high_NETS.txt"), what = "character")

# 获取表达矩阵（使用“data”层，即标准化后的数据）
expr <- GetAssayData(seurat_sub, layer = "data")

# 计算每个细胞中 NETs 基因的平均表达量作为评分
seurat_sub$NETS_score <- colMeans(expr[rownames(expr) %in% nets_genes, ])

# ==================== 4. 根据 NETs 评分将中性粒细胞分为高、低两组 ====================
scores <- seurat_sub$NETS_score[neut_cells]          # 只取中性粒细胞的评分
ord <- order(scores)                                 # 按评分排序
half <- floor(length(ord) / 2)                       # 取一半作为分界

low_neut  <- neut_cells[ord[1:half]]                 # 评分较低的一半
high_neut <- neut_cells[ord[(half + 1):length(ord)]] # 评分较高的一半

# ==================== 5. 为 CellChat 分析定义分组 ====================
# 初始化分组列
seurat_sub$chat_group <- NA

# 中性粒细胞根据评分高低分别标记为 "Neutrophil_Low" 和 "Neutrophil_High"
seurat_sub$chat_group[low_neut]  <- "Neutrophil_Low"
seurat_sub$chat_group[high_neut] <- "Neutrophil_High"

# 肿瘤细胞按其原有的 Seurat 聚类（seurat_clusters）分组，命名为 "Tumor_G0", "Tumor_G1" 等
tumor_clusters <- seurat_sub$seurat_clusters[tumor_cells]
seurat_sub$chat_group[tumor_cells] <- paste0("Tumor_G", tumor_clusters)

# ==================== 6. 构建 CellChat 对象并计算通讯 ====================
# 定义一个函数，给定要分析的细胞列表，自动创建 CellChat 对象并完成基本分析
make_chat <- function(cells_use) {
  # 提取指定细胞的表达矩阵（只保留这些细胞）
  expr <- GetAssayData(seurat_sub, layer = "data")[, cells_use, drop = FALSE]
  
  # 准备元数据：分组信息
  meta <- data.frame(
    group = seurat_sub$chat_group[cells_use],
    row.names = colnames(expr)
  )
  
  # 创建 CellChat 对象
  chat <- createCellChat(expr, meta = meta, group.by = "group")
  
  # 加载人类细胞通讯数据库（确保已加载，如果没有则需提前加载：data("CellChatDB.human")）
  chat@DB <- CellChatDB.human
  
  # 以下为标准 CellChat 分析流程
  chat <- subsetData(chat)                       # 子集化数据，过滤低表达基因
  chat <- identifyOverExpressedGenes(chat)       # 识别过表达的配体/受体基因
  chat <- identifyOverExpressedInteractions(chat) # 识别过表达的相互作用
  chat <- computeCommunProb(chat)                # 计算通讯概率
  chat <- computeCommunProbPathway(chat)         # 计算通路水平通讯概率
  chat <- aggregateNet(chat)                     # 聚合网络，得到 cluster 间总权重
  
  return(chat)
}

# 分别对高 NETs 中性粒细胞 + 肿瘤细胞，以及低 NETs 中性粒细胞 + 肿瘤细胞构建 CellChat
chat_high <- make_chat(c(high_neut, tumor_cells))
chat_low  <- make_chat(c(low_neut,  tumor_cells))

# ==================== 7. 将 CellChat 的 cluster-level 信号分配到单个肿瘤细胞 ====================
# 这一步的目的是获得每个肿瘤细胞从不同中性粒细胞群体接收到的信号强度

# 查看 cluster 名称（确保与后续匹配）
cat("Sender clusters (rows):\n")
print(rownames(chat_high@net$weight))
cat("Receiver clusters (cols):\n")
print(colnames(chat_high@net$weight))

# 确定中性粒细胞群体在 CellChat 对象中的组名
neut_high_cluster <- "Neutrophil_High"
neut_low_cluster  <- "Neutrophil_Low"

# 所有肿瘤组名（以 "Tumor_G" 开头）
tumor_clusters <- grep("^Tumor_G", colnames(chat_high@net$weight), value = TRUE)

# 提取 cluster 间通讯权重（从高/低中性粒到各肿瘤组）
signal_high_cluster <- chat_high@net$weight[neut_high_cluster, tumor_clusters]
signal_low_cluster  <- chat_low@net$weight[neut_low_cluster, tumor_clusters]

# 获取所有肿瘤细胞的 Seurat 聚类 ID（字符形式）
tumor_cells <- WhichCells(seurat_sub, expression = cell_type == "Tumor")
tumor_clusters_char <- as.character(seurat_sub$seurat_clusters[tumor_cells])

# 在 Seurat 对象中初始化两列，用于存储每个肿瘤细胞的信号强度
seurat_sub$neut_signal_high <- NA
seurat_sub$neut_signal_low  <- NA

# 遍历每个肿瘤组，将 cluster 级别的信号值赋给组内所有细胞
for (clu_name in names(signal_high_cluster)) {
  # 去掉 "Tumor_G" 前缀，得到数字编号
  clu_id <- gsub("Tumor_G", "", clu_name)
  
  # 找出属于该肿瘤组的细胞
  cells_in_clu <- tumor_cells[tumor_clusters_char == clu_id]
  
  if (length(cells_in_clu) == 0) {
    warning(clu_name, " 没有对应单细胞，跳过")
    next
  }
  
  # 将高/低组的信号值赋给这些细胞
  seurat_sub$neut_signal_high[cells_in_clu] <- signal_high_cluster[clu_name]
  seurat_sub$neut_signal_low[cells_in_clu]  <- signal_low_cluster[clu_name]
}

# 检查结果
cat("High group signal summary (tumor cells):\n")
print(summary(seurat_sub$neut_signal_high[tumor_cells]))

cat("Low group signal summary (tumor cells):\n")
print(summary(seurat_sub$neut_signal_low[tumor_cells]))

# ==================== 8. 根据信号差重新划分肿瘤细胞群体 ====================
# 计算每个肿瘤细胞在高、低中性粒细胞信号上的差异
signal_diff <- seurat_sub$neut_signal_high[tumor_cells] - 
  seurat_sub$neut_signal_low[tumor_cells]

# 以信号差的中位数作为分界点
median_diff <- median(signal_diff, na.rm = TRUE)

# 初始化新的分组列
seurat_sub$tumor_group <- NA

# 根据信号差是否大于中位数，将肿瘤细胞分为两组
seurat_sub$tumor_group[tumor_cells] <- ifelse(
  signal_diff > median_diff,
  "Tumor_NeutSig_High",   # 对高 NETs 中性粒细胞信号更敏感
  "Tumor_NeutSig_Low"     # 对低 NETs 中性粒细胞信号更敏感
)

# 查看分组结果
table(seurat_sub$tumor_group, useNA = "ifany")

# ==================== 9. 可视化：高低信号散点图 ====================
# 提取肿瘤细胞的高、低信号值
tumor_df <- data.frame(
  high = seurat_sub$neut_signal_high[tumor_cells],
  low  = seurat_sub$neut_signal_low[tumor_cells]
)

# 绘制散点图，添加分位数线（可选，用于辅助划分）
ggplot(tumor_df, aes(x = high, y = low)) +
  geom_point(alpha = 0.5) +
  geom_vline(xintercept = quantile(tumor_df$high, c(1/3, 2/3)), 
             linetype = "dashed", color = "grey") +
  geom_hline(yintercept = quantile(tumor_df$low, c(1/3, 2/3)), 
             linetype = "dashed", color = "grey") +
  labs(
    x = "High-NETS signal (from Neutrophil_High)",
    y = "Low-NETS signal (from Neutrophil_Low)"
  ) +
  theme_classic()

# 保存为文件
ggsave(file.path(out_dir, "signal_scatter.pdf"), width = 6, height = 5)

# 查看每个肿瘤亚群的high和low信号
unique_vals <- unique(data.frame(
  tumor_cluster = seurat_sub$seurat_clusters[tumor_cells],
  high = seurat_sub$neut_signal_high[tumor_cells],
  low = seurat_sub$neut_signal_low[tumor_cells]
))
print(unique_vals)
# ==================== 10. 计算每个肿瘤亚群的平均信号差 ====================
# 信号差：高NETs中性粒细胞信号 - 低NETs中性粒细胞信号
seurat_sub$signal_diff <- seurat_sub$neut_signal_high - seurat_sub$neut_signal_low

# 获取所有肿瘤细胞的亚群信息
tumor_cells <- WhichCells(seurat_sub, expression = cell_type == "Tumor")
tumor_clusters <- seurat_sub$seurat_clusters[tumor_cells]  # 因子或整数

# 计算每个亚群信号差的均值和中位数（用于排序）
cluster_signal_summary <- seurat_sub@meta.data[tumor_cells, ] %>%
  group_by(seurat_clusters) %>%
  summarise(
    mean_diff = mean(signal_diff, na.rm = TRUE),
    median_diff = median(signal_diff, na.rm = TRUE),
    n_cells = n()
  ) %>%
  arrange(desc(mean_diff))  # 按均值降序排列

print(cluster_signal_summary)

# ==================== 11. 选取极端亚群 ====================
# 用户可修改选取的亚群数量，这里示例：信号差最高2个亚群和最低2个亚群
n_extreme <- 2  # 每组选取2个亚群

high_diff_clusters <- cluster_signal_summary$seurat_clusters[1:n_extreme]          # 最高n个
low_diff_clusters <- cluster_signal_summary$seurat_clusters[(nrow(cluster_signal_summary)-n_extreme+1):nrow(cluster_signal_summary)]  # 最低n个

cat("高信号差亚群：", paste(high_diff_clusters, collapse = ", "), "\n")
cat("低信号差亚群：", paste(low_diff_clusters, collapse = ", "), "\n")

# 提取这些极端亚群的细胞
high_sig_cells <- tumor_cells[tumor_clusters %in% high_diff_clusters]
low_sig_cells  <- tumor_cells[tumor_clusters %in% low_diff_clusters]

# 创建新的分组标识（用于后续差异分析）
seurat_sub$extreme_group <- NA
seurat_sub$extreme_group[high_sig_cells] <- "HighSignal_Extreme"
seurat_sub$extreme_group[low_sig_cells]  <- "LowSignal_Extreme"

# 检查各组细胞数量
table(seurat_sub$extreme_group)

# ==================== 12. 差异表达分析 ====================
# 仅使用极端亚群的细胞进行差异分析
Idents(seurat_sub) <- seurat_sub$extreme_group
degs <- FindMarkers(
  seurat_sub,
  ident.1 = "HighSignal_Extreme",
  ident.2 = "LowSignal_Extreme",
  test.use = "wilcox",        # 可选 "wilcox", "bimod", "MAST" 等
  logfc.threshold = 0.25,     # 最小log2FC阈值（过滤）
  min.pct = 0.1,              # 最小表达比例
  only.pos = FALSE            # 同时返回上下调
)

# 添加基因符号列（方便查看）
degs$gene <- rownames(degs)

# 筛选显著差异基因（根据p值调整后阈值）
sig_degs <- degs[degs$p_val_adj < 0.05 & abs(degs$avg_log2FC) > 0.5, ]  # 用户可调整阈值
cat("显著差异基因数量：", nrow(sig_degs), "\n")

# 分别提取上调和下调基因
up_genes <- rownames(sig_degs[sig_degs$avg_log2FC > 0, ])
down_genes <- rownames(sig_degs[sig_degs$avg_log2FC < 0, ])

cat("上调基因数量：", length(up_genes), "\n")
cat("下调基因数量：", length(down_genes), "\n")

# ==================== 12b. 火山图可视化 ====================
library(ggplot2)
library(dplyr)
library(ggrepel)
# 准备绘图数据：添加一列标识显著性/方向
volcano_data <- degs %>%
  mutate(
    # 定义分组：上调（adj.p<0.05 & log2FC > 0.5）、下调（adj.p<0.05 & log2FC < -0.5）、不显著
    regulation = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0.5   ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < -0.5  ~ "Down",
      TRUE                                  ~ "Not significant"
    ),
    # 可选：标记特定基因（例如 top 5 上调/下调）
    gene_label = if_else(
      (regulation == "Up" & rank(-avg_log2FC) <= 5) |
        (regulation == "Down" & rank(avg_log2FC) <= 5),
      gene, NA_character_
    )
  )

# 计算 -log10(adjusted p-value)
volcano_data$log10_padj <- -log10(volcano_data$p_val_adj)

# 绘制火山图
p_volcano <- ggplot(volcano_data, aes(x = avg_log2FC, y = log10_padj, color = regulation)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(
    values = c("Up" = "red", "Down" = "blue", "Not significant" = "grey50"),
    name = "Regulation"
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", alpha = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  labs(
    title = "Volcano Plot: HighSignal_Extreme vs LowSignal_Extreme",
    x = "Log2 Fold Change",
    y = "-Log10(Adjusted P-value)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) +
  # 添加基因标签（避免重叠）
  geom_text_repel(
    aes(label = gene_label),
    max.overlaps = 20,
    size = 3,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50"
  )

# 显示图形
print(p_volcano)

# 保存为 PDF 或 PNG
ggsave(file.path(out_dir, "volcano_DE_High_vs_Low.pdf"), p_volcano, width = 8, height = 6)
ggsave(file.path(out_dir, "volcano_DE_High_vs_Low.png"), p_volcano, width = 8, height = 6, dpi = 300)

cat("火山图已保存至:", out_dir, "\n")

# ==================== 13. 去除假信号 / 免疫噪音 / housekeeping ====================
# 定义需要过滤的基因名称模式（根据实际数据可调整）
bad_prefix <- c("^RPS","^RPL","^MT-","^KRT","^HSP","^EEF","^TUB","^ACT",
                "^IFIT","^ISG","^OAS","^HLA","^S100")

# 1️⃣ 基于差异表达结果筛选出高置信度、高特异性的“真信号”基因
#    这些基因将保留，其余符合模式的基因则被视为噪音
high_conf_genes <- rownames(degs[
  degs$p_val_adj < 1e-10 &
    degs$avg_log2FC > 0.8 &
    degs$pct.1 > 0.2 &
    degs$pct.2 < 0.1, ])

# 2️⃣ 再从中剔除模式匹配的“噪音”基因
genes_to_keep <- high_conf_genes[!grepl(paste(bad_prefix, collapse="|"), high_conf_genes)]

# 3️⃣ 最终的上调基因 = 原上调基因中与“待保留”基因的交集（即保留高置信、非噪音的上调基因）
up_genes_filtered <- intersect(up_genes, genes_to_keep)

cat("高置信度、非噪音上调基因数量:", length(up_genes_filtered), "\n")

# ==================== 14. 保存上调基因及相关信息 ====================
save_file <- file.path(out_dir, "DE_up_genes_filtered.csv")
write.csv(
  degs[rownames(degs) %in% up_genes_filtered, ], 
  file = save_file, 
  row.names = TRUE
)
cat("上调基因及差异信息已保存到:", save_file, "\n")

# ==================== 15. GOBP 富集分析 ====================
library(clusterProfiler)
library(org.Hs.eg.db)

entrez_up <- bitr(up_genes_filtered, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

go_bp_res <- enrichGO(gene=entrez_up$ENTREZID,
                      OrgDb=org.Hs.eg.db,
                      keyType="ENTREZID",
                      ont="BP",
                      pAdjustMethod="BH",
                      pvalueCutoff=0.05,
                      qvalueCutoff=0.2)

go_bp_sig <- as.data.frame(go_bp_res)
go_bp_sig <- go_bp_sig[go_bp_sig$p.adjust < 0.05, ]
cat("显著 GO BP 条目数量:", nrow(go_bp_sig), "\n")
write.csv(go_bp_sig, file.path(out_dir, "GOBP_up_genes_filtered.csv"), row.names=FALSE)

# ==================== 16. GOBP 气泡图（dotplot 版本） ====================
if (nrow(go_bp_sig) > 0) {
  library(stringr)
  p <- dotplot(go_bp_res, 
               showCategory = 20, 
               color = "p.adjust", 
               size = "Count",
               label_format = 45)   # 超过45个字符自动折行
  
  # 保存：增大宽度
  ggsave(file.path(out_dir, "GOBP_Up_dotplot.pdf"), p, width = 10, height = 9)
  ggsave(file.path(out_dir, "GOBP_Up_dotplot.png"), p, width = 10, height = 9, dpi = 300)
  cat("GOBP 气泡图已保存。\n")
} else {
  cat("无显著 GO 条目，跳过绘图。\n")
}
# ==================== 16. 保存 GOBP 气泡图为 PDF ====================
# 确保输出目录存在
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 保存为 PDF 文件
pdf_file <- file.path(out_dir, "GOBP_Up_BubblePlot.pdf")
ggsave(pdf_file, plot = last_plot(), width = 10, height = 9, device = "pdf")

cat("气泡图已保存至：", pdf_file, "\n")


# ==================== 17. CellChat 比较可视化 ====================
library(ggplot2)
library(tidyr)
library(patchwork)   # 用于组合 gg1 和 gg2

# 检查对象存在
if (!exists("chat_high") | !exists("chat_low")) {
  stop("请先运行 CellChat 分析！")
}

# ---------- 图1：总通讯数量与强度组合条形图 ----------
comm_data <- data.frame(
  Group = c("High NETs", "Low NETs"),
  Count = c(sum(chat_high@net$count), sum(chat_low@net$count)),
  Strength = c(sum(chat_high@net$weight), sum(chat_low@net$weight))
)

comm_long <- comm_data %>%
  pivot_longer(cols = c(Count, Strength), names_to = "Metric", values_to = "Value")

p_count_strength <- ggplot(comm_long, aes(x = Group, y = Value, fill = Group)) +
  geom_bar(stat = "identity", width = 0.6) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  labs(title = "Communication quantity and strength: High vs Low NETs") +
  scale_fill_manual(values = c("High NETs" = "#E41A1C", "Low NETs" = "#377EB8")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))

ggsave(file.path(out_dir, "Comm_count_strength_combined.pdf"), 
       p_count_strength, width = 8, height = 4.5)

# ============ 图2：比较每个信号通路的信息流（手动实现） ============
library(ggplot2)
library(dplyr)
library(tidyr)

# 获取两组所有通路的并集
all_pathways <- union(chat_high@netP$pathways, chat_low@netP$pathways)
if (length(all_pathways) == 0) stop("未检测到任何通路！")

# 安全提取通路总强度（第三维为通路）
get_pathway_strength <- function(chat, pathway) {
  prob_array <- chat@netP$prob
  if (is.null(prob_array)) return(NA)
  pathway_names <- dimnames(prob_array)[[3]]   # 通路名在第三维
  if (is.null(pathway_names)) return(NA)
  idx <- which(pathway_names == pathway)
  if (length(idx) == 0) return(0)
  mat <- prob_array[, , idx]   # 发送者×接收者矩阵
  return(sum(mat, na.rm = TRUE))
}

# 计算所有通路在两组中的强度
df_pathway <- data.frame(
  Pathway = all_pathways,
  High = sapply(all_pathways, function(p) get_pathway_strength(chat_high, p)),
  Low  = sapply(all_pathways, function(p) get_pathway_strength(chat_low, p))
)

# 过滤掉两组均为0的通路（无通讯）
df_pathway <- df_pathway[!(df_pathway$High == 0 & df_pathway$Low == 0), ]

# 计算差值（High - Low），并取绝对值排序，选 Top 20（可根据需要调整）
df_pathway$Diff <- df_pathway$High - df_pathway$Low
df_pathway <- df_pathway[order(-abs(df_pathway$Diff)), ]
top20 <- head(df_pathway, 20)

# 将数据转换为长格式，以便绘制堆叠或分面图
top20_long <- top20 %>%
  pivot_longer(cols = c(High, Low), names_to = "Group", values_to = "Strength") %>%
  mutate(Group = factor(Group, levels = c("High", "Low")))

# 按差值排序（使得 High 更高的放在上方）
top20_long$Pathway <- factor(top20_long$Pathway, levels = rev(top20$Pathway))

# 绘制堆叠条形图（类似 rankNet stacked = TRUE 的效果）
p_stacked <- ggplot(top20_long, aes(x = Pathway, y = Strength, fill = Group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("High" = "#E41A1C", "Low" = "#377EB8"),
                    name = "Group") +
  labs(title = "Comparison of signaling pathway information flow",
       x = "", y = "Total communication strength") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 9))

# 绘制分面条形图（类似 rankNet stacked = FALSE，并列显示）
p_dodged <- ggplot(top20_long, aes(x = Pathway, y = Strength, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +
  coord_flip() +
  scale_fill_manual(values = c("High" = "#E41A1C", "Low" = "#377EB8"),
                    name = "Group") +
  labs(title = "Comparison of signaling pathway information flow (dodged)",
       x = "", y = "Total communication strength") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 9))

# 使用 patchwork 组合两个图（可选）
library(patchwork)
p_combined <- p_stacked + p_dodged +
  plot_annotation(title = "Top 20 differential pathways (High vs Low NETs)")

# 保存
ggsave(file.path(out_dir, "Pathway_infoflow_stacked.pdf"), p_stacked, width = 10, height = 8)
ggsave(file.path(out_dir, "Pathway_infoflow_dodged.pdf"), p_dodged, width = 10, height = 8)
ggsave(file.path(out_dir, "Pathway_infoflow_combined.pdf"), p_combined, width = 14, height = 8)

cat("图2：信号通路信息流比较已保存。\n")


# ==================== 18. 肿瘤亚群通讯偏好分组柱状图 ====================
library(ggplot2)
library(dplyr)
library(tidyr)

# 1. 计算每个肿瘤亚群的平均通讯强度（高和低）
tumor_cells <- WhichCells(seurat_sub, expression = cell_type == "Tumor")
tumor_clusters <- seurat_sub$seurat_clusters[tumor_cells]
high_vals <- seurat_sub$neut_signal_high[tumor_cells]
low_vals  <- seurat_sub$neut_signal_low[tumor_cells]

cluster_comm <- data.frame(
  cluster = tumor_clusters,
  high = high_vals,
  low = low_vals
) %>%
  group_by(cluster) %>%
  summarise(
    mean_high = mean(high, na.rm = TRUE),
    mean_low  = mean(low, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    diff = mean_high - mean_low,
    # 标记极端亚群（需要确保 high_diff_clusters 和 low_diff_clusters 已定义）
    extreme_group = case_when(
      cluster %in% high_diff_clusters ~ "High extreme",
      cluster %in% low_diff_clusters  ~ "Low extreme",
      TRUE ~ "Other"
    )
  ) %>%
  arrange(desc(diff))   # 按信号差降序排列

# 保存表格（供论文补充材料）
write.csv(cluster_comm, file.path(out_dir, "Tumor_subcluster_communication_means.csv"), row.names = FALSE)

# 2. 转换为长格式用于绘图
cluster_comm_long <- cluster_comm %>%
  pivot_longer(
    cols = c(mean_high, mean_low),
    names_to = "Group",
    values_to = "Strength"
  ) %>%
  mutate(
    Group = factor(Group, levels = c("mean_high", "mean_low"),
                   labels = c("High NETs", "Low NETs")),
    # 为边框颜色标记极端组
    border_group = ifelse(cluster %in% c(high_diff_clusters, low_diff_clusters),
                          "Extreme", "Other")
  )

# 3. 绘制分组柱状图（每个亚群两根柱）
p_preference <- ggplot(cluster_comm_long,
                       aes(x = factor(cluster, levels = cluster_comm$cluster),
                           y = Strength,
                           fill = Group,
                           color = border_group)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8),
           width = 0.7,
           size = 0.8) +   # 边框粗细
  scale_fill_manual(
    values = c("High NETs" = "#E41A1C", "Low NETs" = "#377EB8"),
    name = "Communication from"
  ) +
  scale_color_manual(
    values = c("Extreme" = "black", "Other" = NA),
    guide = "none"          # 不显示边框图例
  ) +
  labs(
    title = "Differential communication preference of tumor subclusters",
    subtitle = "High vs Low NETs neutrophils",
    x = "Tumor subcluster",
    y = "Mean communication strength (aggregated probability)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45, hjust = 1,
      face = ifelse(levels(factor(cluster_comm$cluster)) %in% c(high_diff_clusters, low_diff_clusters),
                    "bold", "plain")
    ),
    legend.position = "bottom",
    panel.grid.major.x = element_blank()
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")

# 4. 保存图形
ggsave(file.path(out_dir, "Tumor_subcluster_communication_preference_grouped.pdf"),
       p_preference, width = 8, height = 5)
ggsave(file.path(out_dir, "Tumor_subcluster_communication_preference_grouped.png"),
       p_preference, width = 8, height = 5, dpi = 300)

cat("Figure 3e (differential communication preference) 已保存。\n")

