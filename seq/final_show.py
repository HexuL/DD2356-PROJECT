import numpy as np
import matplotlib.pyplot as plt
import os

# 文件路径
file_path = "finitedifference.txt"

# 检查文件是否存在
if os.path.exists(file_path):
    print("Loading data from:", file_path)

    # 加载数据
    data = np.loadtxt(file_path)

    # 创建图像
    fig, ax = plt.subplots()
    im = ax.imshow(data, cmap='bwr', origin='lower')
    plt.colorbar(im)

    # 设置标题
    plt.title("finitedifference.txt")

    # 显示图像
    plt.show()
else:
    print("File finitedifference.txt not found in the current directory.")
