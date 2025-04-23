import os
import pandas as pd
import re
import shutil

# Path to the data folder
data_folder = "../runs/100EV/"

# Dictionary to store DataFrames with unique names
filetypes = {}
filetypes["ExitedParticles"] = {}
filetypes["LostParticles"] = {}
filetypes["MissedCells"] = {}
filetypes["MissedWeights"] = {}
filetypes["Particles"] = {}
filetypes["InitialParticles"] = {}

nparticles = 0

# Walk through all subdirectories and files
for root, _, files in os.walk(data_folder):

    # Skip the root folder
    if root == data_folder:
        continue

    # Only load folders with data from a single run, not combinations of runs, which always end with "Np..."
    if not re.search(r"Np\d+$", root):
        continue

    for file in files:
        if file.endswith("manager.csv"):
            file_path = os.path.join(root, file)
            # Create a unique name for the DataFrame based on the file path
            df_name = os.path.relpath(file_path, data_folder).replace("/", "_").replace("\\", "_").replace(".csv", "_newIds")

            # Get the file type
            filetype = file.split("_")[0]

            # Load the file into a DataFrame
            df = pd.read_csv(file_path)

            # Adjust particle ids to be unique across all files
            df['Particle_id'] += nparticles

            # Store the DataFrame in the dictionary
            filetypes[filetype][df_name] = df

            # Write the DataFrame with unique particle ids
            new_file_path = file_path.replace(".csv", "_newIds.csv")
            df.to_csv(new_file_path, index=False)

        if re.search(r"Np\d+\.csv$", file):
            file_path = os.path.join(root, file)
            print(file_path)

            # Create a unique name for the DataFrame based on the file path
            df_name = os.path.join(root, "InitialParticles")
            df_name = os.path.relpath(df_name, data_folder).replace("/", "_").replace("\\", "_")

            # Load the file into a DataFrame
            df = pd.read_csv(file_path)
            
            # Create Particle_id and adjust it to be unique across all files
            df['Particle_id'] = range(len(df))
            df['Particle_id'] += nparticles

            # Store the DataFrame in the dictionary
            filetypes["InitialParticles"][df_name] = df

            # Write the DataFrame with unique particle ids
            new_file_path = file_path.replace(".csv", "_newIds.csv")
            df.to_csv(new_file_path, index=False)
            
    nparticles += int(root.split("Np")[1])


# Combine DataFrames
# List of all combinations of DataFrames
comb = [
    "Rand",
    "Uniform",
    "Grid",
    "R0.005_R0.14",
    "R0.318_R0.332",
    "Vspile101",
    "Vspile102",
    "Vspile201",
    "Rand_R0.005_R0.14",
    "Rand_R0.318_R0.332",
    "Uniform_R0.005_R0.14",
    "Uniform_R0.318_R0.332",
    "Grid_R0.005_R0.14",
    "Grid_R0.318_R0.332",
    "Rand.*Vspile101",
    "Rand.*Vspile102",
    "Rand.*Vspile201",
    "Uniform.*Vspile101",
    "Uniform.*Vspile102",
    "Uniform.*Vspile201",
    "Grid.*Vspile101",
    "Grid.*Vspile102",
    "Grid.*Vspile201"
]


# Dictionary to store combined DataFrames
for name in comb:
    output_dir = os.path.join(data_folder, name.replace('.*', '_'))
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    print(f"Combining DataFrames for {name}")
    print(output_dir)
    for filetype, dataframes in filetypes.items():
        if filetype == "InitialParticles":
            file_path = os.path.join(output_dir, f"{name.replace('.*', '_')}_newIds.csv")
        else:
            file_path = os.path.join(output_dir, f"{filetype}_1_manager_newIds.csv")

        # Combine DataFrames and save to a new file
        filetypes[filetype][file_path] = pd.concat([df for df_name, df in dataframes.items() if re.search(name, df_name)])
        filetypes[filetype][file_path].to_csv(file_path, index=False)
