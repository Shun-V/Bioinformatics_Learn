# =====================================================================
# scTenifoldKnk 基因敲除分析与可视化测试
# =====================================================================

# 1. 导入必要的 R 包
library(Seurat)
library(Matrix)
library(ggplot2)
library(ggrepel)

# ==========================================
# 2. 基础参数与读取数据
# ==========================================
data_dir <- "./A1.Matrix_0"
target_gene_name <- "Rhbdd1"

cat("正在读取 10X 数据...\n")
counts <- Read10X(data.dir = data_dir)

# 创建 Seurat 对象
seurat_obj <- CreateSeuratObject(counts = counts)

# ==========================================
# 3. 提取高变基因与构建输入矩阵
# ==========================================
cat("正在计算高变基因...\n")
# 经测试，2000基因的计算耗时10分钟；前4分钟在构建计算网络，显示为0%
# 4000基因的网络构建就已经耗时过长，效率起见中断计算
seurat_obj <- FindVariableFeatures(object = seurat_obj, selection.method = "vst", nfeatures = 2000)


# 提取高变基因，并将目标敲除基因强行加入名单并去重
high_variable_genes <- VariableFeatures(seurat_obj)
selected_genes <- unique(c(target_gene_name, high_variable_genes))

# 提取输入矩阵
count_matrix <- GetAssayData(seurat_obj, layer = "counts")[selected_genes, ]

# 清理内存，防止后续多线程计算时 RAM 爆满
rm(counts, seurat_obj)
gc()

# ==========================================
# 4. 运行 scTenifoldKnk 虚拟敲除
# ==========================================
cat(paste("开始对", target_gene_name, "进行虚拟敲除...\n"))
results <- scTenifoldKnk(
    countMatrix = count_matrix,
    gKO = target_gene_name,  
    qc = FALSE,              
    nCores = 4               
)
cat("敲除分析完成！\n")

# ==========================================
# 5. 结果提取与排序
# ==========================================
diff_regulation_df <- results$diffRegulation
diff_regulation_df <- diff_regulation_df[order(diff_regulation_df$p.adj), ]

cat("受影响最严重的前 20 个基因如下：\n")
print(head(diff_regulation_df, 20))

# ==========================================
# 6. 数据可视化展示 (统计学图表)
# ==========================================
cat("开始生成可视化图表...\n")

# ----------------- a. 柱状图 (Top 20) -----------------
# 提取 Fold Change 变化最大的前 20 个基因
top20_diff_genes <- head(diff_regulation_df[order(-abs(diff_regulation_df$FC)), ], 20)
p1 <- ggplot(top20_diff_genes, aes(x=reorder(gene, FC), y=FC)) +
  geom_bar(stat='identity', fill='#5A9BD4') +
  coord_flip() + 
  labs(title="Top 20 Differentially Regulated Genes", x="Gene", y="FC") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

pdf(file="barplot.pdf", width=6, height=5)
print(p1)
dev.off()
cat("柱状图 (barplot.pdf) 已生成！\n")

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
cat("火山图 (volcano.pdf) 已生成！\n")

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
cat("饼图 (pie.pdf) 已生成！\n")