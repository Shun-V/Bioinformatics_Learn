# ==========================================
# 0. 环境准备与包加载 (如果报错找不到包，请先 install.packages 或 BiocManager::install)
# ==========================================
# BiocManager::install(c("clusterProfiler", "org.Hs.eg.db")) # 首次运行需安装富集分析包

library(scTenifoldKnk)
library(Matrix)
library(igraph)
# 制图与富集所需包
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db) # 人类基因注释数据库

# ==========================================
# 1. 读取 Python 导出的三个 10x 标准文件
# ==========================================
cat("正在读入大型稀疏矩阵...\n")
raw_matrix <- readMM("matrix.mtx")

genes <- read.delim("features.tsv", header = FALSE, stringsAsFactors = FALSE)
barcodes <- read.delim("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)

# ==========================================
# 2. 组装矩阵：赋予行名和列名
# ==========================================
rownames(raw_matrix) <- genes[, 1]
colnames(raw_matrix) <- barcodes[, 1]

cat("原始矩阵维度：", dim(raw_matrix)[1], "个基因，", dim(raw_matrix)[2], "个细胞\n")

# ==========================================
# 3. 基因过滤与防报错机制
# ==========================================
cells_expressing_gene <- rowSums(raw_matrix > 0)
min_cells <- ncol(raw_matrix) * 0.05
filtered_matrix <- raw_matrix[cells_expressing_gene >= min_cells, ]

cat("过滤后剩余高置信度基因数：", nrow(filtered_matrix), "\n")

# ==========================================
# 4. 运行 scTenifoldKnk 虚拟敲除
# ==========================================
# TGFB1 的 Ensembl ID 是 "ENSG00000105329"
target_gene <- "ENSG00000105329" 

if (!(target_gene %in% rownames(filtered_matrix))) {
  stop("警告：目标基因在过滤后被剔除了，或者本来就不在矩阵里！请降低过滤阈值或换一个基因。")
}

cat(paste("开始对", target_gene, "(TGFB1) 进行虚拟敲除...\n"))
# 注意：这一步在 M4 Pro 上预计耗时 30分钟 ~ 1.5小时
lung_result <- scTenifoldKnk(countMatrix = filtered_matrix, gKO = target_gene)

# ==========================================
# 5. 结果提取与 Gene ID 翻译 (ENSG -> 基因名)
# ==========================================
diff_regulation_df <- lung_result$diffRegulation
diff_regulation_df <- diff_regulation_df[order(diff_regulation_df$p.adj), ]

# 使用人类数据库翻译基因名，方便作图展示
cat("正在进行基因 ID 翻译...\n")
mapped_symbols <- mapIds(org.Hs.eg.db,
                         keys = diff_regulation_df$gene,
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")

# 将翻译成功的基因名覆盖原有的 ENSG ID，翻译失败的保持原样
diff_regulation_df$gene <- ifelse(is.na(mapped_symbols), diff_regulation_df$gene, mapped_symbols)

cat("敲除分析完成！受影响最严重的前 20 个基因如下：\n")
print(head(diff_regulation_df, 20))


# ==========================================
# 6. 数据可视化展示 (统计学图表)
# ==========================================
cat("开始生成可视化图表...\n")

# ----------------- a. 柱状图 (Top 20) -----------------
top20_diff_genes <- head(diff_regulation_df[order(-diff_regulation_df$FC), ], 20)
p1 <- ggplot(top20_diff_genes, aes(x=reorder(gene, FC), y=FC)) +
  geom_bar(stat='identity', fill='#5A9BD4') +
  coord_flip() + 
  labs(title="Top 20 Differentially Regulated Genes", x="Gene", y="FC") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

pdf(file="barplot.pdf", width=6, height=5)
print(p1)
dev.off()

# ----------------- b. 火山图 (Volcano Plot) -----------------
diff_regulation_df$log_p.adj <- -log10(diff_regulation_df$p.adj + 1e-300) # 加极小值防 log10(0) 报错
diff_regulation_df$significant <- ifelse(diff_regulation_df$p.adj < 0.05, "Significant", "Not significant")
labeled_significant_genes <- subset(diff_regulation_df, p.adj < 0.05)

y_axis_upper_limit <- quantile(diff_regulation_df$log_p.adj, 0.999, na.rm = TRUE)
p2 <- ggplot(diff_regulation_df, aes(x=Z, y=log_p.adj, color=significant)) +
  geom_point(alpha=0.7, size=1) +  
  scale_color_manual(values = c("Significant" = "red", "Not significant" = "gray50")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="red") +
  geom_text_repel(data=labeled_significant_genes, aes(label=gene), size=3, max.overlaps=30) +
  labs(title="Volcano Plot of Virtual Knockout", x="Z-score", y="-log10(FDR)") +
  theme_classic() + 
  coord_cartesian(ylim = c(0, y_axis_upper_limit)) + theme(legend.position = "none")

pdf(file="volcano.pdf", width=6, height=5)
print(p2)
dev.off()

# ----------------- c. 饼图 (Pie Chart) -----------------
significant_gene_count <- table(diff_regulation_df$significant)
significant_gene_df <- as.data.frame(significant_gene_count)
colnames(significant_gene_df) <- c("category", "count")
significant_gene_df$percentage <- paste0(round(significant_gene_df$count / sum(significant_gene_df$count) * 100, 1), "%")

p3 <- ggplot(significant_gene_df, aes(x="", y=count, fill=category)) +
  geom_bar(stat="identity", width=1) +
  geom_text(aes(label=percentage), position=position_stack(vjust=0.5), size=4) + 
  coord_polar("y", start=0) + 
  scale_fill_manual(values=c("Significant"="red", "Not significant"="lightgray")) + 
  labs(title="Proportion of Significant Genes", fill="") +
  theme_minimal() + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text = element_blank(),
        panel.grid = element_blank(), legend.position = "right", plot.title = element_text(hjust = 0.5))

pdf(file="pie.pdf", width=6, height=5)
print(p3)
dev.off()


# ==========================================
# 7. 进阶生物学意义分析 (GO 富集气泡图)
# ==========================================
cat("开始执行 GO 生物学通路富集分析...\n")

# 获取显著变化的基因 ID (使用原始的 ENSG ID，以便数据库能精准识别)
# 我们这里重新从原始提取，因为上面为了作图把 diff_regulation_df$gene 替换成了 Symbol
sig_ensembl_ids <- rownames(lung_result$diffRegulation)[lung_result$diffRegulation$p.adj < 0.05]

# 如果显著基因太少可能无法富集，做一个保护判断
if (length(sig_ensembl_ids) > 5) {
  ego <- enrichGO(gene          = sig_ensembl_ids,
                  keyType       = "ENSEMBL",
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",        # BP: Biological Process
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05)
  
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    pdf("GO_enrichment_dotplot.pdf", width=8, height=6)
    # 画气泡图，展示最受影响的前 15 条生物学通路
    print(dotplot(ego, showCategory=15, title="Top 15 Affected Biological Processes"))
    dev.off()
    cat("富集气泡图已生成！请查看 GO_enrichment_dotplot.pdf\n")
  } else {
    cat("提示：富集分析未发现显著聚集的通路 (可能网络震荡过于分散)。\n")
  }
} else {
  cat("提示：显著受影响的基因数量太少 (<5)，不足以进行富集分析。\n")
}

cat("运行完毕！\n")