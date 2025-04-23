import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import ScalarFormatter

# Set the style for the plots
sns.set(style="whitegrid", palette="muted")

# Load the CSV files
particle_path = '/home/annah/semesterproject/code/ippl/build/semesterproject/data/Particles_1_manager.csv'  # Replace with the correct path if needed

particles_df = pd.read_csv(particle_path)

# Columns to plot histograms for
columns_to_plot = ['Time', 'E_kin']

# Calculate E_kin in eV
particles_df['E_kin'] = particles_df['E_kin'] / 1.602e-19  # Convert from Joules to eV

# Filter the dataframe by 'particle' column without changing the order
unique_particles = particles_df['Particle_id'].unique()

# Calculate Energy difference to initial energy for each particle
for particle in unique_particles:
    # Filter the dataframe by particle
    particle_df = particles_df[particles_df['Particle_id'] == particle]
    # Calculate initial energy
    initial_energy = particle_df['E_kin'].iloc[10]
    # Create a new column for energy difference
    particles_df.loc[particles_df['Particle_id'] == particle, 'E_kin_diff'] = (particle_df['E_kin'] - initial_energy) / initial_energy



# Set up the figures
plt.figure(figsize=(14, 4))

# Plot E_kin_diff for each particle
for particle in unique_particles:
    # Filter the dataframe by particle
    particle_df = particles_df[particles_df['Particle_id'] == particle]
    # Filter out rows where velocities are zero
    particle_df = particle_df[(particle_df['Velocity_x'] != 0) | (particle_df['Velocity_y'] != 0) | (particle_df['Velocity_z'] != 0)]
    plt.plot(particle_df['Time'][1:], particle_df['E_kin_diff'][1:], alpha=0.8, color=[0/255, 73/255, 104/255], linewidth=2)

plt.xlabel(r'Time[s]', fontsize=16)
plt.ylabel(r'$(E_{\mathrm{kin}} - E_{\mathrm{kin, t=0}}) / E_{\mathrm{kin, t=0}} [eV]$', fontsize=16)
plt.grid(True)
plt.gca().yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
plt.gca().ticklabel_format(axis='y', style='sci', scilimits=(-3, 3))  # Use scientific notation for y-axis ticks
plt.gca().ticklabel_format(axis='x', style='sci', scilimits=(-1, 3))  # Use scientific notation for x-axis ticks

# Adjust the y-axis ticks to make them more readable
plt.gca().tick_params(axis='y', labelsize=12)

# Adjust layout to make sure everything fits
plt.tight_layout()

# Save the plot as an image file
plt.savefig('E_kin_diff.png', dpi=300)


# Set up the figures
plt.figure(figsize=(14, 4))

# Plot E_kin for each particle
for particle in unique_particles:
    # Filter the dataframe by particle
    particle_df = particles_df[particles_df['Particle_id'] == particle]
    # Filter out rows where velocities are zero
    particle_df = particle_df[(particle_df['Velocity_x'] != 0) | (particle_df['Velocity_y'] != 0) | (particle_df['Velocity_z'] != 0)]
    plt.plot(particle_df['Time'][1:], particle_df['E_kin'][1:] - particle_df['E_kin'][0], alpha=0.8, linewidth=2, color=[0/255, 73/255, 104/255])
plt.xlabel(r'Time[s]', fontsize=16)
plt.ylabel(r'E_kin[eV]', fontsize=16)
plt.grid(True)
plt.gca().yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
plt.gca().ticklabel_format(axis='y', style='sci', scilimits=(-3, 3))  # Use scientific notation for y-axis ticks
plt.gca().ticklabel_format(axis='x', style='sci', scilimits=(-1, 3))  # Use scientific notation for x-axis ticks

# Adjust the y-axis ticks to make them more readable
plt.gca().tick_params(axis='y', labelsize=12)

# Adjust layout to make sure everything fits
plt.tight_layout()

# Save the plot as an image file
plt.savefig('E_kin.png', dpi=300)