import sys
import os
import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV file
# Determine the folder path
folder_path = sys.argv[1] if len(sys.argv) > 1 else '../../build/sempro/'

# Construct the full file path
file_path = os.path.join(folder_path, 'curl.csv')

# Load the CSV file
data = pd.read_csv(file_path)

# Column to plot
column_to_plot = sys.argv[2] if len(sys.argv) > 2 else 'VorticityMagnitude'

# Generate the histogram
plt.figure(figsize=(8, 6))
plt.hist(data[column_to_plot], bins=30, color='blue', edgecolor='black', alpha=0.7)
plt.title('Histogram of Vorticity Magnitude')
plt.xlabel('Vorticity Magnitude')
plt.ylabel('Frequency')
plt.grid(True, linestyle='--', alpha=0.6)

# Show the plot
plt.show()