import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV files
particle_path = '/home/anna/semesterproject/anekslen_ippl/build_parallel/semesterproject/data/Particles_1_manager.csv'  # Replace with the correct path if needed

particles_df = pd.read_csv(particle_path)

def cut_dataframe(df, time, particle_ids = []):
    if len(particle_ids) != 0:
        return df[(df['Particle_id'].isin(particle_ids)) & (df['Time'] <= time)]
    return df[(df['Time'] <= time)]

# Function to subtract values from Particle_id 0
def subtract_particle_0(df, pmin_id, maxtime):
    result_df = pd.DataFrame()
    for i, (time, group) in enumerate(df.groupby('Time')):
        if time <= maxtime:
            if i % 1000 == 0:
                print(f"Subtracting values from Particle_id min at time {time}")
            particle_0 = group[group['Particle_id'] == pmin_id].iloc[0]
            group['Position_x'] -= particle_0['Position_x']
            group['Position_x'] /= particle_0['Position_x']
            group['Position_y'] -= particle_0['Position_y']
            group['Position_y'] /= particle_0['Position_y']
            group['Position_z'] -= particle_0['Position_z']
            group['Position_z'] /= particle_0['Position_z']
            group['Velocity_x'] -= particle_0['Velocity_x']
            group['Velocity_x'] /= particle_0['Velocity_x']
            group['Velocity_y'] -= particle_0['Velocity_y']
            group['Velocity_y'] /= particle_0['Velocity_y']
            group['Velocity_z'] -= particle_0['Velocity_z']
            group['Velocity_z'] /= particle_0['Velocity_z']
            group['B_x'] -= particle_0['B_x']
            group['B_x'] /= particle_0['B_x']
            group['B_y'] -= particle_0['B_y']
            group['B_y'] /= particle_0['B_y']
            group['B_z'] -= particle_0['B_z']
            group['B_z'] /= particle_0['B_z']
            result_df = pd.concat([result_df, group])
    print("Subtracted values from Particle_id min")
    return result_df

def plot_columns(df, columns, title, xlabel, ylabel):
    plt.figure(figsize=(16, 10))
    for particle in df['Particle_id'].unique():
        particle_df = df[df['Particle_id'] == particle]
        plt.plot(particle_df['Time'], particle_df[columns], alpha=0.7)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(df['Particle_id'].unique())
    plt.savefig(f'{title}.png')

# Find the particle_id of the particle whose velocity changes to 0 at the latest time step
def find_particle_with_latest_zero_velocity(df):
    latest_zero_velocity_time = -1
    particle_id_with_latest_zero_velocity = None
    
    for particle_id in df['Particle_id'].unique():
        particle_df = df[df['Particle_id'] == particle_id]
        zero_velocity_times = particle_df[(particle_df['Velocity_x'] == 0) & 
                                            (particle_df['Velocity_y'] == 0) & 
                                            (particle_df['Velocity_z'] == 0)]['Time']
        if not zero_velocity_times.empty:
            latest_zero_velocity_time_for_particle = zero_velocity_times.min()
            if latest_zero_velocity_time_for_particle > latest_zero_velocity_time:
                latest_zero_velocity_time = latest_zero_velocity_time_for_particle
                particle_id_with_latest_zero_velocity = particle_id

    if latest_zero_velocity_time == -1:
        latest_zero_velocity_time = df['Time'].max()
        particle_id_with_latest_zero_velocity = df['Particle_id'].max()
    
    return particle_id_with_latest_zero_velocity, latest_zero_velocity_time

def plot_multiple_particles(df, columns, particle_ids = []):
    if len(particle_ids) == 0:
        particle_ids = df['Particle_id'].unique()

    for column in columns:
        fig, axs = plt.subplots(len(particle_ids), 1, figsize=(16, 10 * len(particle_ids)), sharex=True)
        if len(particle_ids) == 1:
            axs = [axs]
        for ax, particle_id in zip(axs, particle_ids):
            particle_df = df[df['Particle_id'] == particle_id]
            ax.plot(particle_df['Time'], particle_df[column], alpha=0.7)
            ax.set_title(f'{column} for Particle {particle_id}')
            ax.set_xlabel('Time')
            ax.set_ylabel(column)
        plt.tight_layout()
        plt.savefig(f'{column}_multiple.png')

particle_id_with_latest_zero_velocity, latest_time = find_particle_with_latest_zero_velocity(particles_df)
print(f'Particle ID with the latest zero velocity: {particle_id_with_latest_zero_velocity} at time {latest_time}')

# Apply the function to the dataframe
# particles_df = subtract_particle_0(particles_df, particle_id_with_latest_zero_velocity, latest_time)

# Get data frame with only 2 particles and time <= latest time
particles_df_cut = cut_dataframe(particles_df, latest_time, [particle_id_with_latest_zero_velocity-1, particle_id_with_latest_zero_velocity])
particles_df_timecut = cut_dataframe(particles_df, latest_time)

plot_columns(particles_df_cut, 'Position_x', 'Position_x', 'Time', 'Position_x')
plot_columns(particles_df_cut, 'Position_y', 'Position_y', 'Time', 'Position_y')
plot_columns(particles_df_cut, 'Position_z', 'Position_z', 'Time', 'Position_z')
plot_columns(particles_df_cut, 'Velocity_x', 'Velocity_x', 'Time', 'Velocity_x')
plot_columns(particles_df_cut, 'Velocity_y', 'Velocity_y', 'Time', 'Velocity_y')
plot_columns(particles_df_cut, 'Velocity_z', 'Velocity_z', 'Time', 'Velocity_z')
plot_columns(particles_df_cut, 'B_x', 'B_x', 'Time', 'B_x')
plot_columns(particles_df_cut, 'B_y', 'B_y', 'Time', 'B_y')
plot_columns(particles_df_cut, 'B_z', 'B_z', 'Time', 'B_z')

plot_multiple_particles(particles_df_timecut, ['Position_x','B_x', 'B_y', 'B_z'], [particle_id_with_latest_zero_velocity-1, particle_id_with_latest_zero_velocity])
plot_multiple_particles(particles_df_timecut, ['Position_z'])