import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os

# 获取当前目录下所有txt文件
folder_path = "./outputdata"
files = [os.path.join(folder_path, file) for file in os.listdir(folder_path) if file.endswith(".txt")]
files.sort()  # 排序文件以确保按顺序显示
print("show ",len(files),"data")

# 创建一个函数来更新图像数据
def update_data(file):
    # 加载数据
    data = np.loadtxt(file)

    # 更新图像数据
    im.set_array(data)

    # 更新标题
    plt.title(os.path.basename(file))

    return im,

# 创建初始图像
fig, ax = plt.subplots()
im = ax.imshow(np.zeros((10, 10)), cmap='bwr', origin='lower')
plt.colorbar(im)

# 创建动画
ani = FuncAnimation(fig, update_data, frames=files, blit=True, interval=50)
plt.show()
