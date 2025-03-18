import pandas as pd
import matplotlib.pyplot as plt
import sys

# Load the CSV files
particle_path = sys.argv[1] if len(sys.argv) > 1 else '/home/annah/semesterproject/code/ippl/build/semesterproject/data/LostParticles_1_manager.csv'

particles_df = pd.read_csv(particle_path)
# Count the percentage of particles exiting at each boundary region
boundary_counts = particles_df['Boundary_region'].value_counts()
boundary_counts = boundary_counts / boundary_counts.sum() * 100


# Bar plot of the exit positions
plt.figure(figsize=(16, 10))
boundary_counts.plot(kind='bar')
plt.xlabel('Boundary Region')
plt.ylabel('Percentage of Particles')
plt.title('Particles exiting at each boundary region')
plt.savefig('ExitPositions.png')