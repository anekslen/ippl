import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot_histogram(dataframe, xlabel, output_file=None, normalize=False, new = True):
    for column in columns:
        plt.figure(figsize=(16, 10))
        weights = np.ones_like(dataframe[column]) / len(dataframe) if normalize else None
        if new:
            plt.hist(dataframe[column], bins=100, color=[128/255, 0/255, 32/255], edgecolor='black', alpha=0.7, weights=weights, log=True)
            plt.legend(["New Mesh"], fontsize=20)
        else:
            plt.hist(dataframe[column], bins=100, color=[0/255, 73/255, 104/255], edgecolor='black', alpha=0.7, weights=weights, log=True)
            plt.legend(["Old Mesh"], fontsize=20)
        plt.xlabel(xlabel, fontsize=20)
        plt.ylabel('Relative Frequency [%]', fontsize=20)
        plt.tight_layout()
        if output_file:
            plt.savefig(f"{output_file}_{column}.png")
        else:
            plt.show()

# Function to plot comparison histograms
def plot_comparison_histograms(dataframe1, dataframe2, xlabel, output_file=None, single_output_file=None, normalize=False):
    for column in columns:
        plt.figure(figsize=(16, 10))
        weights1 = np.ones_like(dataframe1[column]) / len(dataframe1) if normalize else None
        weights2 = np.ones_like(dataframe2[column]) / len(dataframe2) if normalize else None
        min_val = min(dataframe1[column].min(), dataframe2[column].min())
        max_val = max(dataframe1[column].max(), dataframe2[column].max())
        bins = np.linspace(min_val, max_val, 100)
        plt.hist(dataframe2[column], bins=bins, color=[0/255, 73/255, 104/255], alpha=0.7, label='Old Mesh', edgecolor='black', weights=weights2, log=True)
        plt.hist(dataframe1[column], bins=bins, color=[128/255, 0/255, 32/255], alpha=0.7, label='New Mesh', edgecolor='black', weights=weights1, log=True)
        plt.xlabel(xlabel, fontsize=25)
        plt.ylabel('Relative  Frequency [%]', fontsize=25)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        plt.legend(fontsize=25)
        if single_output_file:
            plt.savefig(f"{single_output_file}_{column}.png")
        else:
            plt.show()



# Load the CSV files
curl_new_path = '/home/annah/semesterproject/code/ippl/semesterproject/mesh/reactormesh/vtkInput/NewMesh/curl.csv'
grad_new_path = '/home/annah/semesterproject/code/ippl/semesterproject/mesh/reactormesh/vtkInput/NewMesh/grad.csv'
curl_old_path = '/home/annah/semesterproject/code/ippl/semesterproject/mesh/reactormesh/vtkInput/OldMesh/curl.csv'
grad_old_path = '/home/annah/semesterproject/code/ippl/semesterproject/mesh/reactormesh/vtkInput/OldMesh/grad.csv'

curl_new_df = pd.read_csv(curl_new_path)
curl_old_df = pd.read_csv(curl_old_path)
grad_new_df = pd.read_csv(grad_new_path)
grad_old_df = pd.read_csv(grad_old_path)

# Columns to plot histograms for
columns = ["Magnitude"]

norm = True

# Generate individual histograms for curl
plot_histogram(curl_new_df, xlabel = r'$|\nabla \times \mathbf{B}|$ [T/m]', output_file='curl_new', normalize=norm, new=True)
plot_histogram(curl_old_df, xlabel = r'$|\nabla \times \mathbf{B}|$ [T/m]', output_file='curl_old', normalize=norm, new=False)
# Generate comparison histograms for curl
plot_comparison_histograms(curl_new_df, curl_old_df, xlabel=r'$|\nabla \times \mathbf{B}|$ [T/m]', output_file='curl_comparison', single_output_file='curl_comparison', normalize=norm)

# Generate individual histograms for grad
plot_histogram(grad_new_df, xlabel=r'$|\nabla \cdot \mathbf{B}|$ [T/m]', output_file='grad_new', normalize=norm, new=True)
plot_histogram(grad_old_df, xlabel=r'$|\nabla \cdot \mathbf{B}|$ [T/m]', output_file='grad_old', normalize=norm, new=False)
# Generate comparison histograms for grad
plot_comparison_histograms(grad_new_df, grad_old_df, xlabel=r'$|\nabla \cdot \mathbf{B}|$ [T/m]', output_file='grad_comparison', single_output_file='grad_comparison', normalize=norm)