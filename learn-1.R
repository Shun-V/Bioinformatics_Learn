# 加载包
library(scTenifoldKnk)
library(Matrix)
library(igraph)

# 设定随机种子保证每次运行结果一致 (类似 np.random.seed)
set.seed(42) 

# 设定矩阵维度：100个基因，300个细胞
n_genes <- 100
n_cells <- 300

# 生成负二项分布数据
raw_counts <- matrix(rnbinom(n_genes * n_cells, size = 1, prob = 0.1), nrow = n_genes)

# 命名行 (基因) 和列 (细胞)
rownames(raw_counts) <- paste0("Gene_", 1:n_genes)
colnames(raw_counts) <- paste0("Cell_", 1:n_cells)

# 转换为稀疏矩阵格式 (dgCMatrix，类似 scipy.sparse.csc_matrix)
sparse_counts <- as(raw_counts, "dgCMatrix")

# 查看矩阵维度和前 5行5列
print(dim(sparse_counts))
print(sparse_counts[1:5, 1:5])

# 使用稀疏矩阵，敲除目标是 "Gene_10"
cat("开始运行 scTenifoldKnk，请稍候...\n")
result <- scTenifoldKnk(countMatrix = sparse_counts, gKO = "Gene_10")
cat("运行完毕！\n")

# 提取差异调控结果
affected_genes <- result$diffRegulation

# 按照统计学显著性 (p.adj) 从小到大排序
affected_genes <- affected_genes[order(affected_genes$p.adj), ]

# 查看受敲除冲击最严重的前 10 个基因
head(affected_genes, 10)

# 绘制 Gene_10 敲除后的核心响应网络
plotKO(result, gKO = "Gene_10", annotate = FALSE)

# 保存图片
# 第一步：打开一个名叫 "Gene10_KO_network.png" 的画布文件
# 你可以设置宽度(width)、高度(height)和分辨率(res)
png("Gene10_KO_network.png", width = 800, height = 800, res = 120)
# pdf("Gene10_KO_network.pdf", width = 7, height = 7)

# 第二步：执行画图命令（将图片在画布上渲染）
plotKO(result, gKO = "Gene_10", annotate = FALSE)

# 第三步：关闭画布并保存（非常重要！不跑这句图片文件会损坏）
dev.off()