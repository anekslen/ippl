import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV files
particle_path = '/home/anna/semesterproject/anekslen_ippl/build_new/semesterproject/data/Particles_1_manager.csv'  # Replace with the correct path if needed

particles_df = pd.read_csv(particle_path)

# Columns to plot histograms for
columns_to_plot = ['Time', 'E_kin']

# Filter the dataframe by 'particle' column without changing the order
unique_particles = particles_df['Particle_id'].unique()

# Calculate Energy difference to initial energy for each particle
for particle in unique_particles:
    # Filter the dataframe by particle
    particle_df = particles_df[particles_df['Particle_id'] == particle]
    # Calculate initial energy
    initial_energy = particle_df['E_kin'].iloc[10]
    print(initial_energy)
    # Create a new column for energy difference
    particles_df.loc[particles_df['Particle_id'] == particle, 'E_kin_diff'] = (particle_df['E_kin'] - initial_energy) / initial_energy

# Plot E_kin_diff for each particle
plt.figure(figsize=(16, 10))
plt.subplot(1, 2, 1)
for particle in unique_particles:
    # Filter the dataframe by particle
    particle_df = particles_df[particles_df['Particle_id'] == particle]
    # Filter out rows where velocities are zero
    particle_df = particle_df[(particle_df['Velocity_x'] != 0) | (particle_df['Velocity_y'] != 0) | (particle_df['Velocity_z'] != 0)]
    plt.plot(particle_df['Time'][1:], particle_df['E_kin_diff'][1:], alpha=0.7)
plt.title('E_kin_diff')
plt.xlabel('Time')
plt.ylabel('E_kin_diff')
plt.legend(unique_particles)

plt.subplot(1, 2, 2)
for particle in unique_particles:
    # Filter the dataframe by particle
    particle_df = particles_df[particles_df['Particle_id'] == particle]
    # Filter out rows where velocities are zero
    particle_df = particle_df[(particle_df['Velocity_x'] != 0) | (particle_df['Velocity_y'] != 0) | (particle_df['Velocity_z'] != 0)]
    plt.plot(particle_df['Time'][1:], particle_df['E_kin'][1:], alpha=0.7)
plt.title('E_kin')
plt.xlabel('Time')
plt.ylabel('E_kin')
plt.legend(unique_particles)

plt.savefig('E_kin.png')

# Plot E_kin for each particle in single plots
plt.figure(figsize=(16, 10))
for i, particle in enumerate(unique_particles, 1):
    plt.subplot(2, 2, i)
    # Filter the dataframe by particle
    particle_df = particles_df[particles_df['Particle_id'] == particle]
    # Filter out rows where velocities are zero
    particle_df = particle_df[(particle_df['Velocity_x'] != 0) | (particle_df['Velocity_y'] != 0) | (particle_df['Velocity_z'] != 0)]
    plt.plot(particle_df['Time'][1:], particle_df['E_kin_diff'][1:], alpha=0.7)
    plt.title('E_kin')
    plt.xlabel('Time')
    plt.ylabel('E_kin')
    plt.legend()
plt.savefig('E_kin_single.png')