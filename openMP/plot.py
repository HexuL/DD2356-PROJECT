import numpy as np
import matplotlib.pyplot as plt
import os

## @file plot.py
#  @brief This script loads data from a file and visualizes it using a heatmap.
# 
#  This script checks if the specified file exists, loads the data from the file,
#  and then creates a heatmap using matplotlib. It also adds a color bar and sets
#  the title of the plot to the filename.

file_path = "output.txt"

## @brief Check if the specified file exists.
#
#  If the file does not exist, print an error message and exit the program.
if not os.path.exists(file_path):
    print(f"File {file_path} does not exist.")
    exit()

## @brief Load data from the specified file.
#
#  Use numpy to load the data from the file into a 2D array.
data = np.loadtxt(file_path)

## @brief Create a new figure and axis object for plotting.
fig, ax = plt.subplots()

## @brief Display the data as a heatmap.
#
#  Use the 'bwr' colormap for the heatmap and set the origin to the lower left corner.
im = ax.imshow(data, cmap='bwr', origin='lower')

## @brief Add a color bar to the plot.
plt.colorbar(im)

## @brief Set the title of the plot.
#
#  The title is set to the basename of the file path.
plt.title(os.path.basename(file_path))

## @brief Display the plot.
plt.show()

