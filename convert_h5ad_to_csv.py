import scanpy as sc
import scipy.io as sio
import pandas as pd

print("1. 读取 h5ad 文件...")
adata = sc.read_h5ad("af6e81be-e65c-4821-987e-e0eb6c8acd59.h5ad")

# CELLxGENE 的数据通常把原始 counts 藏在 adata.raw.X 里
if adata.raw is not None:
    counts = adata.raw.X
    genes = adata.raw.var_names
else:
    counts = adata.X
    genes = adata.var_names
    
cells = adata.obs_names

print("2. 开始导出为 10x 格式 (这可能需要一两分钟)...")
# 核心避坑：Python 的矩阵是 [细胞 x 基因]，而 R 是 [基因 x 细胞]
# 所以导出前必须用 .T 进行矩阵转置！
sio.mmwrite("matrix.mtx", counts.T)

# 导出基因名和细胞名
pd.Series(genes).to_csv("features.tsv", index=False, header=False)
pd.Series(cells).to_csv("barcodes.tsv", index=False, header=False)

print("3. 转换完成，输出保存在features.tsv和barcodes.tsv！")