import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
import numpy as np

# Load output data
# Load the CSV files
particle_path = sys.argv[1] if len(sys.argv) > 1 else '/home/annah/semesterproject/code/ippl/build/semesterproject/data/LostParticles_1_manager_newIds.csv'
folder_path = os.path.dirname(particle_path)

# Load the particle data
particles_df = pd.read_csv(particle_path)
wrong_df = particles_df[particles_df['Boundary_type'] == 'Interior']
particles_df = particles_df[particles_df['Boundary_type'] != 'Interior'] # Discard particles that were kicked out of the simulation due to negative interpolation weights

print(f"Percentage of lost particles: {len(wrong_df) / len(particles_df) * 100:.2f}%")

# Prepare data for the bar plot
# Count the percentage of particles exiting at each boundary region
boundary_counts = particles_df['Boundary_type'].value_counts()
total_particles = boundary_counts.sum()

# Calculate standard error for the bar plot
error_bars = np.sqrt((boundary_counts / total_particles) * (1 - (boundary_counts / total_particles)) / total_particles) * 100

# Normalize the counts to get percentages
boundary_counts = boundary_counts / total_particles * 100



# Calculate coordinates for the boundary lines and circle segments for inner and outer boundary regions
# Define inner and outer boundary regions (inner [0], outer [1])
RMin = [0.005, 0.318]
RMax = [0.14, 0.332]

ThetaMin = [0, 20]
ThetaMin = [value / 360 * 2 * np.pi for value in ThetaMin]

ThetaMax = [45, 45]
ThetaMax = [value / 360 * 2 * np.pi for value in ThetaMax]

# Calculate Points for circle segments and lines
X_RMin_ThetaMin = [RMin[0] * np.cos(ThetaMin[0]), RMin[1] * np.cos(ThetaMin[1])]
Y_RMin_ThetaMin = [RMin[0] * np.sin(ThetaMin[0]), RMin[1] * np.sin(ThetaMin[1])]

X_RMin_ThetaMax = [RMin[0] * np.cos(ThetaMax[0]), RMin[1] * np.cos(ThetaMax[1])]
Y_RMin_ThetaMax = [RMin[0] * np.sin(ThetaMax[0]), RMin[1] * np.sin(ThetaMax[1])]

X_RMax_ThetaMin = [RMax[0] * np.cos(ThetaMin[0]), RMax[1] * np.cos(ThetaMin[1])]
Y_RMax_ThetaMin = [RMax[0] * np.sin(ThetaMin[0]), RMax[1] * np.sin(ThetaMin[1])]

X_RMax_ThetaMax = [RMax[0] * np.cos(ThetaMax[0]), RMax[1] * np.cos(ThetaMax[1])]
Y_RMax_ThetaMax = [RMax[0] * np.sin(ThetaMax[0]), RMax[1] * np.sin(ThetaMax[1])]

# Calculate boundary lines of input particle area
lines = [
    [((X_RMin_ThetaMin[0], Y_RMin_ThetaMin[0]), (X_RMax_ThetaMin[0], Y_RMax_ThetaMin[0])),
    ((X_RMin_ThetaMax[0], Y_RMin_ThetaMax[0]), (X_RMax_ThetaMax[0], Y_RMax_ThetaMax[0]))],
    [((X_RMin_ThetaMin[1], Y_RMin_ThetaMin[1]), (X_RMax_ThetaMin[1], Y_RMax_ThetaMin[1])),
    ((X_RMin_ThetaMax[1], Y_RMin_ThetaMax[1]), (X_RMax_ThetaMax[1], Y_RMax_ThetaMax[1]))]
]

# Calculate circle segments of input particle area
theta = [np.linspace(ThetaMin[0], ThetaMax[0], 100), np.linspace(ThetaMin[1], ThetaMax[1])]  # Circle segments

X_RMin_circle = [RMin[0] * np.cos(theta[0]), RMin[1] * np.cos(theta[1])]
y_RMin_circle = [RMin[0] * np.sin(theta[0]), RMin[1] * np.sin(theta[1])]

X_RMax_circle = [RMax[0] * np.cos(theta[0]), RMax[1] * np.cos(theta[1])]
y_RMax_circle = [RMax[0] * np.sin(theta[0]), RMax[1] * np.sin(theta[1])]

# Check which regions need to be plotted for the file (if not specified in folder name both inner and outer regions are included)
inner = True if 'R0.005_R0.14' in folder_path else False
outer = True if 'R0.318_R0.332' in folder_path else False

if not inner and not outer:
    inner = True
    outer = True



# Load coordinates of the start points of the particles
# Load the CSV file
input_path = sys.argv[2] if len(sys.argv) > 2 else '/home/annah/semesterproject/code/ippl/build/semesterproject/data/Rand_newIds.csv'
input_df = pd.read_csv(input_path)

# Convert polar coordinates (R, Theta, Z) of initial positions to Cartesian coordinates (x, y, z)
if 'R' in input_df.columns and 'Theta' in input_df.columns and 'Z' in input_df.columns:
    input_df['x'] = input_df['R'] * np.cos(input_df['Theta'])
    input_df['y'] = input_df['R'] * np.sin(input_df['Theta'])
    input_df['z'] = input_df['Z']

if 'X' in input_df.columns and 'Y' in input_df.columns and 'Z' in input_df.columns:
    input_df.rename(columns={'X': 'x', 'Y': 'y', 'Z': 'z'}, inplace=True)

# Join the two DataFrames on the 'Particle_id' column
merged_df = pd.merge(input_df, particles_df[['Particle_id', 'Time', 'Boundary_type']], on='Particle_id', how='left')
merged_df = merged_df.dropna(subset=['Boundary_type'])



# Create color map for the bar plot
# Read all the filenames in the boundary folder and add them to a list
boundary_folder = sys.argv[3] if len(sys.argv) > 3 else '/home/annah/semesterproject/code/ippl/semesterproject/mesh/reactormesh/vtkInput/NewMesh/BoundaryClassification/'
boundary_files = [os.path.splitext(f)[0] for f in os.listdir(boundary_folder) if os.path.isfile(os.path.join(boundary_folder, f))]

# Color map for the bar plot
dark_cud_colors = [
    (102/255, 51/255, 0/255),
    (0/255, 73/255, 104/255),
    (176/255, 32/255, 36/255),
    (204/255, 85/255, 0/255),
    (103/255, 35/255, 103/255),
    (51/255, 102/255, 0/255),
    (0/255, 109/255, 119/255),
    (70/255, 130/255, 180/255),
    (153/255, 102/255, 0/255),
    (47/255, 79/255, 79/255),
    (128/255, 0/255, 32/255)
]

# Map boundary files to colors
boundary_color_map = {boundary: dark_cud_colors[i % len(dark_cud_colors)] for i, boundary in enumerate(boundary_files)}
ax_colors = [boundary_color_map.get(boundary, (0.5, 0.5, 0.5)) for boundary in boundary_counts.index]

# Define dotsizes and dpis for the scatter plots (each scatter plot will be saved with different dot sizes and dpis)
dpis = [500]
dotsize = [1, 2, 5]



# Plot 1: Bar plot of the exit positions
plt.figure(figsize=(16, 10))
ax = boundary_counts.plot(kind='bar', color=ax_colors, edgecolor='black', linewidth=2, yerr=error_bars, capsize=12, alpha=0.9)

# Add percentage labels to the bars
for i, value in enumerate(boundary_counts):
    ax.text(i, value + (error_bars.iloc[i] + 0.5), f'{value:.2f}%', ha='center', fontsize=20, color='black', fontweight='bold')

# Increase y-axis limit to create space for labels
y_max = ax.get_ylim()[1]
ax.set_ylim(0, y_max + 1)

plt.xlabel('Boundary Region', fontsize=25)
plt.ylabel('Fraction of Particles [%]', fontsize=25)
plt.xticks(rotation=45, ha='right', fontsize=25)
plt.yticks(fontsize=25)
plt.tight_layout()  # Adjust layout to prevent clipping
# Add a legend with boundary names and their corresponding colors
legend_labels = [f'{boundary}' for boundary in boundary_color_map.keys()]
legend_patches = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=boundary_color_map[boundary], markersize=23) for boundary in boundary_color_map.keys()]
plt.legend(legend_patches, legend_labels, loc='upper right', fontsize=23)
plt.savefig(f'{folder_path}/ExitPositions.png')



# Plot 2: Histogram plot of the exit times
plt.figure(figsize=(16, 10))

# Calculate histogram data
hist_data, bin_edges = np.histogram(
    particles_df['Time'], 
    bins=100, 
    weights=(100 * np.ones(len(particles_df)) / len(particles_df))
)

plt.bar(
    bin_edges[:-1], 
    hist_data, 
    width=np.diff(bin_edges), 
    log=True, 
    color=[0/255, 109/255, 119/255], 
    edgecolor='black', 
    linewidth=1.2
)


plt.xlabel('Time [s]', fontsize=25)
plt.ylabel('Fraction of Particles [%]', fontsize=25)
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.tight_layout()  # Adjust layout to prevent clipping
plt.savefig(f'{folder_path}/ExitTimes.png')

# Calculate the percentage of particles that survive longer than 0.00014s
surviving_particles = particles_df[particles_df['Time'] > 0.00014]
percentage_surviving = (len(surviving_particles) / len(particles_df)) * 100
print(f"Percentage of particles surviving longer than 0.00014s: {percentage_surviving:.2f}%")



# Plot 3: Scatter plot of particle positions colored by boundary region
for sz in dotsize:
    plt.figure(figsize=(10, 8))
    # Ensure equal aspect ratio for x and y axes
    plt.gca().set_aspect('equal', adjustable='box')

    for boundary, color in boundary_color_map.items():
        subset = merged_df[merged_df['Boundary_type'] == boundary]
        plt.scatter(subset['x'], subset['y'], label=boundary, color=color, alpha=0.7, s=sz)
    legend_labels = [f'{boundary}' for boundary in boundary_color_map.keys()]
    legend_patches = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=boundary_color_map[boundary], markersize=14) for boundary in boundary_color_map.keys()]
    plt.legend(legend_patches, legend_labels, fontsize=14)
    plt.xlabel('x [m]', fontsize=25)
    plt.ylabel('y [m]', fontsize=25)

    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)

    # Plot lines and circle segments of input particle area
    if inner:
        plt.plot(X_RMin_circle[0], y_RMin_circle[0], color='black', linewidth=1)
        plt.plot(X_RMax_circle[0], y_RMax_circle[0], color='black', linewidth=1)

        for line in lines[0]:
            x_values, y_values = zip(*line)
            plt.plot(x_values, y_values, color='black', linewidth=1)

    if outer:
        plt.plot(X_RMin_circle[1], y_RMin_circle[1], color='black', linewidth=1)
        plt.plot(X_RMax_circle[1], y_RMax_circle[1], color='black', linewidth=1)

        for line in lines[1]:
            x_values, y_values = zip(*line)
            plt.plot(x_values, y_values, color='black', linewidth=1)

    for dpi in dpis:
        plt.savefig(f'{folder_path}/ParticlePositions_{sz}_{dpi}.png', dpi=dpi)



# Plot 4: Scatter plot of particle positions colored by exit time
for sz in dotsize:
    plt.figure(figsize=(10, 8))
    # Ensure equal aspect ratio for x and y axes
    plt.gca().set_aspect('equal', adjustable='box')

    plt.scatter(merged_df['x'], merged_df['y'], c=merged_df['Time'], cmap='viridis', s=sz)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=20)  # Adjust the size of the text in the color bar
    cbar.ax.set_ylabel('Exit Time [s]', fontsize=25)
    plt.xlabel('x', fontsize=25)
    plt.ylabel('y', fontsize=25)

    plt.xlabel('x [m]', fontsize=25)
    plt.ylabel('y [m]', fontsize=25)

    # Plot lines and circle segments of input particle area
    if inner:
        plt.plot(X_RMin_circle[0], y_RMin_circle[0], color='black', linewidth=1)
        plt.plot(X_RMax_circle[0], y_RMax_circle[0], color='black', linewidth=1)
        
        for line in lines[0]:
            x_values, y_values = zip(*line)
            plt.plot(x_values, y_values, color='black', linewidth=1)

    if outer:
        plt.plot(X_RMin_circle[1], y_RMin_circle[1], color='black', linewidth=1)
        plt.plot(X_RMax_circle[1], y_RMax_circle[1], color='black', linewidth=1)
        
        for line in lines[1]:
            x_values, y_values = zip(*line)

            plt.plot(x_values, y_values, color='black', linewidth=1)

    for dpi in dpis:
        plt.savefig(f'{folder_path}/ParticlePositionsTime_{sz}_{dpi}.png', dpi=dpi)



# Plot 5: Scatter plot of particles colored by boundary region for particles having times greater than 0.002s
longConf_df = merged_df[merged_df['Time'] > 0.002]
for sz in dotsize:
    plt.figure(figsize=(10, 8))

    # Ensure equal aspect ratio for x and y axes
    plt.gca().set_aspect('equal', adjustable='box')

    for boundary, color in boundary_color_map.items():
        subset = longConf_df[longConf_df['Boundary_type'] == boundary]
        plt.scatter(subset['x'], subset['y'], label=boundary, color=color, alpha=0.7, s=sz)

    legend_labels = [f'{boundary}' for boundary in boundary_color_map.keys()]
    legend_patches = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=boundary_color_map[boundary], markersize=15) for boundary in boundary_color_map.keys()]
    plt.legend(legend_patches, legend_labels, fontsize=15)
    plt.xlabel('x [m]', fontsize=25)
    plt.ylabel('y [m]', fontsize=25)

    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)

    # Plot lines and circle segments of input particle area
    if inner:
        plt.plot(X_RMin_circle[0], y_RMin_circle[0], color='black', linewidth=1)
        plt.plot(X_RMax_circle[0], y_RMax_circle[0], color='black', linewidth=1)

        for line in lines[0]:
            x_values, y_values = zip(*line)
            plt.plot(x_values, y_values, color='black', linewidth=1)

    if outer:
        plt.plot(X_RMin_circle[1], y_RMin_circle[1], color='black', linewidth=1)
        plt.plot(X_RMax_circle[1], y_RMax_circle[1], color='black', linewidth=1)
        
        for line in lines[1]:
            x_values, y_values = zip(*line)
            plt.plot(x_values, y_values, color='black', linewidth=1)

    for dpi in dpis:
        plt.savefig(f'{folder_path}/ParticlePositionsLongConf_{sz}_{dpi}.png', dpi=dpi)



# Plot 6: Scatter plot of particles exiting at the tunel
for sz in dotsize:
    tunel_df = merged_df[merged_df['Boundary_type'].isin(['Tun_Slant_E', 'Tun_Roof_E', 'Tun_Fillet_E', 'Tun_Side_E'])]

    plt.figure(figsize=(10, 8))
    # Ensure equal aspect ratio for x and y axes
    plt.gca().set_aspect('equal', adjustable='box')

    for boundary, color in boundary_color_map.items():
        subset = tunel_df[tunel_df['Boundary_type'] == boundary]
        plt.scatter(subset['x'], subset['y'], label=boundary, color=color, alpha=0.7, s=sz)

    legend_labels = [f'{boundary}' for boundary in boundary_color_map.keys() if boundary in ['Tun_Slant_E', 'Tun_Roof_E', 'Tun_Fillet_E', 'Tun_Side_E']]
    legend_patches = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=boundary_color_map[boundary], markersize=15) for boundary in legend_labels]
    plt.legend(legend_patches, legend_labels, fontsize=15)
    plt.xlabel('x [m]', fontsize=25)
    plt.ylabel('y [m]', fontsize=25)

    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)

    # Plot lines and circle segments of input particle area
    if inner:
        plt.plot(X_RMin_circle[0], y_RMin_circle[0], color='black', linewidth=1)
        plt.plot(X_RMax_circle[0], y_RMax_circle[0], color='black', linewidth=1)
        
        for line in lines[0]:
            x_values, y_values = zip(*line)
            plt.plot(x_values, y_values, color='black', linewidth=1)

    if outer:
        plt.plot(X_RMin_circle[1], y_RMin_circle[1], color='black', linewidth=1)
        plt.plot(X_RMax_circle[1], y_RMax_circle[1], color='black', linewidth=1)
        
        for line in lines[1]:
            x_values, y_values = zip(*line)
            plt.plot(x_values, y_values, color='black', linewidth=1)
        
    for dpi in dpis:
        plt.savefig(f'{folder_path}/ParticlePositionsTunel_{sz}_{dpi}.png', dpi=dpi)