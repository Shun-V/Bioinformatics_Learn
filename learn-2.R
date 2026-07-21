# 部分读取dataset并转换到R的依赖
# install.packages("Seurat")
# install.packages("devtools")
# devtools::install_github("mojaveazure/seurat-disk")

# 加载包
library(scTenifoldKnk)
library(Seurat)
library(SeuratDisk)
library(Matrix)
library(igraph)

# ==========================
# 该方式只适用于旧版的h5ad文件，新版结构有所变化，无法读取

# # 读取以下文件
# h5ad_file <- "af6e81be-e65c-4821-987e-e0eb6c8acd59.h5ad"
# h5seurat_file <- "lung.h5seurat"

# # 1. 格式转换：从 Python 的 h5ad 转成 R 的 h5seurat
# # 这个过程可能会打印很多进度条，如果提示是否覆盖，输入 y 或设置 overwrite = TRUE
# Convert(h5ad_file, dest = "h5seurat", overwrite = TRUE)

# # 2. 读取到内存中，成为一个 Seurat 对象
# lung_seurat <- LoadH5Seurat(h5seurat_file, meta.data = FALSE, misc = FALSE)

# # 看看导入了什么
# print(lung_seurat)

# ==========================
# 更换为自己的转换脚本

# ==========================================
# 1. 读取 Python 导出的三个 10x 标准文件
# ==========================================
cat("正在读入大型稀疏矩阵...\n")
# readMM 专门用来瞬间读取极大的 .mtx 稀疏矩阵文件
raw_matrix <- readMM("matrix.mtx")

# 读取基因名和细胞名 (header = FALSE 因为我们用 Python 导出时没有表头)
genes <- read.delim("features.tsv", header = FALSE, stringsAsFactors = FALSE)
barcodes <- read.delim("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)

# ==========================================
# 2. 组装矩阵：赋予行名和列名
# ==========================================
# 将 genes 的第一列设为行名，barcodes 的第一列设为列名
rownames(raw_matrix) <- genes[, 1]
colnames(raw_matrix) <- barcodes[, 1]

# ==========================================
# 3. 基本维度查看 (行是基因，列是细胞)
# ==========================================
cat("原始矩阵维度：", dim(raw_matrix)[1], "个基因，", dim(raw_matrix)[2], "个细胞\n")

# 查看一下组装好的矩阵前 5 行 5 列，确认行名列名都已经贴上去了
print(raw_matrix[1:5, 1:5])

# ==========================================
# 4. 开始基因敲除测试
# ==========================================

# 1. 极度重要的防报错机制：基因过滤
# 计算每个基因在多少个细胞中表达量大于 0
cells_expressing_gene <- rowSums(raw_matrix > 0)

# 设定阈值：至少在 5% 的细胞中表达
min_cells <- ncol(raw_matrix) * 0.05

# 保留满足条件的基因
filtered_matrix <- raw_matrix[cells_expressing_gene >= min_cells, ]
cat("过滤后剩余高置信度基因数：", nrow(filtered_matrix), "\n")


# 2. 检查 TGFB1 的 Ensembl ID 是否在矩阵中
# TGFB1 的 Ensembl ID 是 "ENSG00000105329"
target_gene <- "ENSG00000105329" 

if (!(target_gene %in% rownames(filtered_matrix))) {
  stop("警告：目标基因在过滤后被剔除了，或者本来就不在矩阵里！请降低过滤阈值或换一个基因。")
}


# 3. 开始执行虚拟敲除分析
cat(paste("开始对", target_gene, "(TGFB1) 进行虚拟敲除，这可能需要很长时间...\n"))

# 运行核心算法
lung_result <- scTenifoldKnk(countMatrix = filtered_matrix, gKO = target_gene)


# 4. 提取并查看受波及最严重的前 20 个基因
final_table <- lung_result$diffRegulation
final_table <- final_table[order(final_table$p.adj), ]

print("敲除分析完成！受影响最严重的前 20 个基因如下：")
print(head(final_table, 20))