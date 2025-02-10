import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV file
data = pd.read_csv('curl.csv')

# Column to plot
column_to_plot = 'VorticityMagnitude'

# Generate the histogram
plt.figure(figsize=(8, 6))
plt.hist(data[column_to_plot], bins=30, color='blue', edgecolor='black', alpha=0.7)
plt.title('Histogram of Vorticity Magnitude')
plt.xlabel('Vorticity Magnitude')
plt.ylabel('Frequency')
plt.grid(True, linestyle='--', alpha=0.6)

# Show the plot
plt.show()
