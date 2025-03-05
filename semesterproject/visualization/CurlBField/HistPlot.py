import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the CSV files
curl_new_path = 'curlNew.csv'  # Replace with the correct path if needed
curl_old_path = 'curlOld.csv'

curl_new_df = pd.read_csv(curl_new_path)
curl_old_df = pd.read_csv(curl_old_path)

# Columns to plot histograms for
columns_to_plot = ['Vorticity_X', 'Vorticity_Y', 'Vorticity_Z', 'VorticityMagnitude']
# Rename columns
curl_new_df.rename(columns={'Vorticity_X': 'Magnetic curl X', 'Vorticity_Y': 'Magnetic curl Y', 'Vorticity_Z': 'Magnetic curl Z', 'VorticityMagnitude': 'Magnetic curl Magnitude'}, inplace=True)
curl_old_df.rename(columns={'Vorticity_X': 'Magnetic curl X', 'Vorticity_Y': 'Magnetic curl Y', 'Vorticity_Z': 'Magnetic curl Z', 'VorticityMagnitude': 'Magnetic curl Magnitude'}, inplace=True)

# Update columns to plot
columns_to_plot = ['Magnetic curl X', 'Magnetic curl Y', 'Magnetic curl Z', 'Magnetic curl Magnitude']

# Function to plot individual histograms
def plot_individual_histograms(dataframe, title_prefix, output_file=None, normalize=False):
    plt.figure(figsize=(16, 10))
    for i, column in enumerate(columns_to_plot, 1):
        plt.subplot(2, 2, i)
        weights = np.ones_like(dataframe[column]) / len(dataframe) if normalize else None
        plt.hist(dataframe[column], bins=100, color='skyblue', edgecolor='black', alpha=0.7, weights= weights, log=True)
        plt.title(f'{title_prefix} {column}')
        plt.xlabel(f'{column} in T/m')
        plt.ylabel('Frequency')
    plt.tight_layout()
    if output_file:
        plt.savefig(output_file)
    else:
        plt.show()

# Function to plot comparison histograms
def plot_comparison_histograms(dataframe1, dataframe2, title1, title2, output_file=None, normalize=False):
    plt.figure(figsize=(16, 10))
    for i, column in enumerate(columns_to_plot, 1):
        plt.subplot(2, 2, i)
        weights1 = np.ones_like(dataframe1[column]) / len(dataframe1) if normalize else None
        weights2 = np.ones_like(dataframe2[column]) / len(dataframe2) if normalize else None
        min_val = min(dataframe1[column].min(), dataframe2[column].min())
        max_val = max(dataframe1[column].max(), dataframe2[column].max())
        bins = np.linspace(min_val, max_val, 100)
        plt.hist(dataframe2[column], bins=bins, color='orange', alpha=0.6, label=title2, edgecolor='black', weights=weights2, log=True)
        plt.hist(dataframe1[column], bins=bins, color='skyblue', alpha=0.6, label=title1, edgecolor='black', weights=weights1, log=True)
        plt.title(f'Comparison: {column}')
        plt.xlabel(f'{column} in T/m')
        plt.ylabel('Frequency')
        plt.legend()
    plt.tight_layout()
    if output_file:
        plt.savefig(output_file)
    else:
        plt.show()


norm = True

# Generate individual histograms
plot_individual_histograms(curl_new_df, "New Mesh -", output_file="histograms_new_data.png", normalize=norm)
print("1")
plot_individual_histograms(curl_old_df, "Old Mesh -", output_file="histograms_old_data.png", normalize=norm)
print("1")

# Generate comparison histograms
plot_comparison_histograms(curl_new_df, curl_old_df, "New Mesh", "Old Mesh", output_file="comparison_histograms.png", normalize=norm)
print("1")