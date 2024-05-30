import numpy as np
import matplotlib.pyplot as plt
import os

file_path = "output.txt"

# Check if file exists
if not os.path.exists(file_path):
    print(f"File {file_path} does not exist.")
    exit()

# Load data from file
data = np.loadtxt(file_path)

# Create a new figure and axis object
fig, ax = plt.subplots()

# Display data using 'bwr' colormap, with origin at the lower left corner
im = ax.imshow(data, cmap='bwr', origin='lower')

# Add a color bar
plt.colorbar(im)

# Set the title of the figure to the filename
plt.title(os.path.basename(file_path))

# Show the plot
plt.show()

