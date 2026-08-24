# ============================================================
#  二分法 CellChat 分析（Inner vs Outer）
#  动态划分：按层级数量均分，Center 归入 Inner
#  输出：四个方向的配体通讯表格、条形图、netVisual_circle 网络图
#  修正版：消除 samples 警告、统一条形宽度、使用 netVisual_circle
# ============================================================

rm(list = ls())
library(Seurat)
library(CellChat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(writexl)
library(readxl)

# ----- 1. 配置 -----
CONFIG <- list(
  base_dir      = "E:/NETS/空转组数据/结果/Lesion_Analysis/Oscillation_Field",
  output_dir    = "E:/NETS/空转组数据/结果/Lesion_Analysis/CellChat_Analysis/Binary_Results",
  layer_order   = c("Center", "L6", "L5", "L4", "L3", "L2", "L1"),
  p_threshold   = 0.05,
  top_ligands   = 10,
  min_cells     = 10,
  raw_use       = TRUE,
  compute_prob_type = "truncatedMean",
  trim          = 0.1,
  dpi           = 300
)

if (!dir.exists(CONFIG$output_dir)) dir.create(CONFIG$output_dir, recursive = TRUE)

# ----- 2. 辅助函数 -----

#' 为 Seurat 对象添加层级标签（layer_detailed）
define_layers <- function(obj, layer_order = CONFIG$layer_order) {
  meta <- obj@meta.data
  grade_vec <- meta$oscillation_grade
  lesion_vec <- meta$lesion_id
  
  valid_idx <- !is.na(grade_vec) & grade_vec != "Zero"
  if (sum(valid_idx) == 0) {
    meta$layer_detailed <- NA
    obj@meta.data <- meta
    return(obj)
  }
  
  exist_levels <- intersect(layer_order, unique(grade_vec[valid_idx]))
  if (length(exist_levels) == 0) {
    meta$layer_detailed <- NA
    obj@meta.data <- meta
    return(obj)
  }
  
  meta$layer_detailed <- NA
  meta$layer_detailed[!is.na(lesion_vec)] <- "Center"
  for (lv in setdiff(exist_levels, "Center")) {
    meta$layer_detailed[grade_vec == lv & is.na(lesion_vec)] <- lv
  }
  obj@meta.data <- meta
  return(obj)
}

#' 动态二分分组：按层级顺序均分，Center 必属 Inner
get_binary_groups <- function(exist_levels) {
  n <- length(exist_levels)
  if (n < 2) return(list(Inner = exist_levels, Outer = character(0)))
  center_pos <- which(exist_levels == "Center")
  if (length(center_pos) == 0) {
    split_idx <- ceiling(n / 2)
  } else {
    split_idx <- ceiling(n / 2)
  }
  Inner <- exist_levels[1:split_idx]
  Outer <- exist_levels[(split_idx + 1):n]
  list(Inner = Inner, Outer = Outer)
}

#' 运行 CellChat 并提取四个方向的通讯（返回 dirs 和 cellchat 对象）
run_binary_cellchat <- function(counts, meta_sub, config) {
  # 添加 samples 列，消除警告
  meta_sub$samples <- "sample1"
  
  # 创建 CellChat 对象，按 layer_detailed 分组
  cellchat <- createCellChat(counts, meta = meta_sub, group.by = "layer_detailed")
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  if (config$compute_prob_type == "triMean") {
    cellchat <- computeCommunProb(cellchat, raw.use = config$raw_use)
  } else {
    cellchat <- computeCommunProb(cellchat, 
                                  type = "truncatedMean", 
                                  trim = config$trim,
                                  raw.use = config$raw_use)
  }
  cellchat <- filterCommunication(cellchat, min.cells = config$min_cells)
  cellchat <- computeCommunProbPathway(cellchat)
  
  # 提取详细层级的通讯矩阵
  net_prob <- cellchat@net$prob
  net_pval <- cellchat@net$pval
  if (is.null(net_prob) || is.null(net_pval)) return(NULL)
  
  senders <- dimnames(net_prob)[[1]]
  receivers <- dimnames(net_prob)[[2]]
  pathways <- dimnames(net_prob)[[3]]
  
  # 获取二分分组
  exist_levels <- intersect(CONFIG$layer_order, unique(meta_sub$layer_detailed))
  groups <- get_binary_groups(exist_levels)
  Inner <- groups$Inner
  Outer <- groups$Outer
  if (length(Outer) == 0) {
    warning("没有 Outer 层级，跳过")
    return(NULL)
  }
  
  # 构建发送者、接收者的二分标签
  sender_binary <- ifelse(senders %in% Inner, "Inner", ifelse(senders %in% Outer, "Outer", NA))
  receiver_binary <- ifelse(receivers %in% Inner, "Inner", ifelse(receivers %in% Outer, "Outer", NA))
  
  # 初始化四个方向的列表
  dirs <- list()
  for (s in c("Inner", "Outer")) {
    for (r in c("Inner", "Outer")) {
      dir_name <- paste0(s, "_to_", r)
      dirs[[dir_name]] <- data.frame()
    }
  }
  
  # 遍历所有配体-受体对
  for (i in seq_along(senders)) {
    for (j in seq_along(receivers)) {
      s_bin <- sender_binary[i]
      r_bin <- receiver_binary[j]
      if (is.na(s_bin) || is.na(r_bin)) next
      dir_name <- paste0(s_bin, "_to_", r_bin)
      for (p in seq_along(pathways)) {
        prob <- net_prob[i, j, p]
        pval <- net_pval[i, j, p]
        if (!is.na(prob) && prob > 0 && pval < config$p_threshold) {
          dirs[[dir_name]] <- rbind(dirs[[dir_name]], data.frame(
            Sender = senders[i],
            Receiver = receivers[j],
            Pathway = pathways[p],
            Probability = prob,
            P_val = pval,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  # 汇总各方向配体强度（按接收层细胞数归一化）
  cellcounts <- table(meta_sub$layer_detailed)
  for (dir_name in names(dirs)) {
    df <- dirs[[dir_name]]
    if (nrow(df) == 0) next
    df$Ligand <- gsub("_.*", "", df$Pathway)
    df$ReceiverLayer <- df$Receiver
    df$CellCount <- cellcounts[df$ReceiverLayer]
    df_sum <- df %>%
      group_by(ReceiverLayer, Ligand) %>%
      summarise(
        TotalProb = sum(Probability, na.rm = TRUE),
        MinPval = min(P_val, na.rm = TRUE),
        CellCount = first(CellCount),
        .groups = "drop"
      ) %>%
      mutate(AvgProb = TotalProb / CellCount) %>%
      ungroup()
    dirs[[dir_name]] <- df_sum %>%
      group_by(Ligand) %>%
      summarise(
        TotalProb = sum(TotalProb, na.rm = TRUE),
        AvgProb = sum(AvgProb, na.rm = TRUE),
        MinPval = min(MinPval, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(Direction = dir_name)
  }
  
  # 返回 dirs 和 cellchat 对象
  return(list(dirs = dirs, cellchat = cellchat))
}

#' 绘制条形图（分面展示四个方向的 Top 配体）- 统一 x 轴范围，固定条形宽度
plot_binary_bar <- function(df_merged, sample_name, config) {
  plot_data <- df_merged %>%
    group_by(Direction) %>%
    arrange(desc(AvgProb), .by_group = TRUE) %>%
    slice_head(n = config$top_ligands) %>%
    ungroup()
  
  if (nrow(plot_data) == 0) return(NULL)
  
  # 计算 -log10(AvgProb)
  plot_data$NegLog10 <- -log10(plot_data$AvgProb + 1e-12)
  plot_data$Label <- format(plot_data$AvgProb, scientific = TRUE, digits = 2)
  
  # 按方向排序配体
  plot_data <- plot_data %>%
    group_by(Direction) %>%
    mutate(Ligand = factor(Ligand, levels = Ligand[order(NegLog10, decreasing = FALSE)])) %>%
    ungroup()
  
  # 计算全局最大 x 值，统一坐标轴范围
  max_x <- max(plot_data$NegLog10, na.rm = TRUE) * 1.1
  
  p <- ggplot(plot_data, aes(x = NegLog10, y = Ligand, fill = Direction)) +
    geom_col(width = 0.6, alpha = 0.85) +
    geom_text(aes(label = Label), hjust = -0.1, size = 2.8, colour = "black") +
    scale_fill_manual(values = c("Inner_to_Inner" = "#E41A1C", 
                                 "Inner_to_Outer" = "#F78104",
                                 "Outer_to_Inner" = "#4DAF4A",
                                 "Outer_to_Outer" = "#377EB8")) +
    scale_x_continuous(
      name = "-log10(Avg Communication Probability)",
      expand = expansion(mult = c(0, 0.15)),
      limits = c(0, max_x)
    ) +
    facet_wrap(~ Direction, scales = "free_y", ncol = 2) +
    theme_bw(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(face = "bold", size = 11),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.3, "lines"),
      axis.text.y = element_text(size = 8, face = "bold"),
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
      axis.title.x = element_text(size = 10),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13)
    ) +
    labs(
      title = paste0(sample_name, " - Binary Communication (Top ", config$top_ligands, " ligands per direction)"),
      x = "-log10(Avg Communication Probability)",
      y = NULL
    )
  
  return(p)
}

# ----- 3. 主流程 -----

main <- function() {
  # 加载数据
  seurat_list <- readRDS(file.path(CONFIG$base_dir, "seurat_list_with_weights.rds"))
  sample_names <- names(seurat_list)
  
  # 为每个样本定义详细层级
  for (s in sample_names) {
    seurat_list[[s]] <- define_layers(seurat_list[[s]], layer_order = CONFIG$layer_order)
  }
  
  all_results <- list()  # 保存各样本的方向汇总
  
  for (s in sample_names) {
    cat("\n========================================\n")
    cat("🔄 处理样本:", s, "\n")
    
    obj <- seurat_list[[s]]
    meta <- obj@meta.data
    valid <- !is.na(meta$layer_detailed)
    if (sum(valid) < 10) {
      message("   ⚠️ 有效 spots 太少，跳过")
      next
    }
    
    counts <- GetAssayData(obj, assay = "Spatial", layer = "counts")[, valid, drop = FALSE]
    meta_sub <- meta[valid, ]
    
    exist_levels <- intersect(CONFIG$layer_order, unique(meta_sub$layer_detailed))
    if (length(exist_levels) < 2) {
      message("   ⚠️ 层级数 < 2，跳过")
      next
    }
    groups <- get_binary_groups(exist_levels)
    cat("   Inner:", paste(groups$Inner, collapse = ", "), "\n")
    cat("   Outer:", paste(groups$Outer, collapse = ", "), "\n")
    
    # 运行 CellChat 并提取四个方向（同时返回 cellchat 对象）
    res <- run_binary_cellchat(counts, meta_sub, CONFIG)
    if (is.null(res)) {
      message("   ⚠️ 无显著通讯，跳过")
      next
    }
    dirs <- res$dirs
    cellchat <- res$cellchat
    
    # 合并所有方向到一个数据框
    df_merged <- bind_rows(dirs, .id = "Direction")
    if (nrow(df_merged) == 0) {
      message("   ⚠️ 无通讯数据")
      next
    }
    
    # 保存详细表格
    write_xlsx(df_merged, path = file.path(CONFIG$output_dir, paste0("Binary_Summary_", s, ".xlsx")))
    
    # 生成条形图（统一 x 轴）
    p_bar <- plot_binary_bar(df_merged, s, CONFIG)
    if (!is.null(p_bar)) {
      ggsave(file.path(CONFIG$output_dir, paste0("Binary_Bar_", s, ".pdf")), 
             p_bar, width = 12, height = 10, dpi = CONFIG$dpi)
    }
    
    # ---- 新增：使用 netVisual_circle 绘制二分网络图 ----
    # 从 cellchat 对象中构建 Inner/Outer 聚合矩阵（总概率）
    # ---- 新增：使用 netVisual_circle 绘制二分网络图 ----
    net_prob <- cellchat@net$prob
    senders <- dimnames(net_prob)[[1]]
    receivers <- dimnames(net_prob)[[2]]
    Inner <- groups$Inner
    Outer <- groups$Outer
    
    sender_bin <- ifelse(senders %in% Inner, "Inner", ifelse(senders %in% Outer, "Outer", NA))
    receiver_bin <- ifelse(receivers %in% Inner, "Inner", ifelse(receivers %in% Outer, "Outer", NA))
    
    agg_mat <- matrix(0, nrow = 2, ncol = 2, 
                      dimnames = list(c("Inner","Outer"), c("Inner","Outer")))
    for (i in seq_along(senders)) {
      for (j in seq_along(receivers)) {
        sender_label <- sender_bin[i]
        receiver_label <- receiver_bin[j]
        if (is.na(sender_label) || is.na(receiver_label)) next
        total_prob <- sum(net_prob[i, j, ], na.rm = TRUE)
        agg_mat[sender_label, receiver_label] <- agg_mat[sender_label, receiver_label] + total_prob
      }
    }
    
    weight_inner <- sum(meta_sub$layer_detailed %in% Inner)
    weight_outer <- sum(meta_sub$layer_detailed %in% Outer)
    vertex_weight <- c(weight_inner, weight_outer)
    names(vertex_weight) <- c("Inner", "Outer")
    
    # 绘制并保存（此时 s 仍是正确的样本名）
    pdf(file.path(CONFIG$output_dir, paste0("Binary_Ring_", s, ".pdf")), 
        width = 6, height = 6)
    netVisual_circle(agg_mat, 
                     vertex.weight = vertex_weight,
                     weight.scale = TRUE, 
                     label.edge = TRUE,
                     title.name = paste0(s, " - Binary Communication (Total Probability)"))
    dev.off()
    
    # （可选）如果您仍想保留原来自定义的环形图，可以取消下面注释
    # p_ring <- plot_binary_ring(df_merged, s, CONFIG)
    # if (!is.null(p_ring)) {
    #   ggsave(file.path(CONFIG$output_dir, paste0("Binary_Ring_", s, ".pdf")), 
    #          p_ring, width = 10, height = 10, dpi = CONFIG$dpi)
    # }
    
    all_results[[s]] <- df_merged
    cat("   ✅ 样本处理完成\n")
  }
  
  # 汇总所有样本的结果
  if (length(all_results) > 0) {
    combined <- bind_rows(all_results, .id = "Sample")
    write_xlsx(combined, path = file.path(CONFIG$output_dir, "AllSamples_Binary_Summary.xlsx"))
    cat("\n✅ 汇总文件已保存\n")
  }
  
  cat("\n🎉 所有样本处理完成！\n")
  cat("📁 输出目录:", CONFIG$output_dir, "\n")
}

# 执行
main()

