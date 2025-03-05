import numpy as np
import csv

# Constants
dt = 1.0e-9  # Time step

# Magnetic field (Tesla)
B_field = np.array([0.0, 0.0, -0.14])
B_array = np.array([[0, B_field[2], -B_field[1]],
                   [-B_field[2], 0, B_field[0]],
                   [B_field[1], -B_field[0], 0]])

# Particle properties
mass = 1.673e-27  # kg
charge = 1.602e-19  # Coulombs

E_0 = 100 * 1.602e-19 # energy in Joules
absv_0 = np.sqrt(2.0 * E_0 / mass)  # absolut value of speed in m/s
print(f"absv_0: {absv_0}")

position = np.array([0.0, 0.0, 0.0])  # meters
velocity = np.array([absv_0 / np.sqrt(2), absv_0 / np.sqrt(2), -10000])  # meters/second
print(f"velocity: {velocity[0], velocity[1], velocity[2]}")

# Function to compute the Lorentz force
def lorentz_force(charge, velocity, B_field):
    return charge * np.cross(velocity, B_field)

# Leapfrog integration scheme
def leapfrog_step(position, velocity, mass, charge, B_field, dt):    
    # Full-step position update
    new_position = position + dt * velocity
    
    # Full-step velocity update
    new_velocity = velocity + dt * lorentz_force(charge, velocity, B_field) / mass
    
    return new_position, new_velocity

def tajimasImplicitMethod_step(position, velocity, mass, charge, B_matrix, dt):
    # Full-step position update
    new_position = position + dt * velocity
    
    # Full-step velocity update
    new_velocity = np.linalg.solve(np.eye(3) - dt * charge / (2 * mass) * B_matrix, (np.eye(3) + dt * charge / (2 * mass) * B_matrix).dot(velocity))
    
    return new_position, new_velocity

def boris_step(position, velocity, mass, charge, B_field, dt):
    
    new_position = position + dt * velocity

    t = 0.5 * charge * B_field / mass * dt
    s = 2.0 * t / (1.0 + np.dot(t, t))
    
    v_prime = velocity + np.cross(velocity, t)
    new_velocity = velocity + np.cross(v_prime, s)

    
    return new_position, new_velocity

# Simulation loop
num_steps = 10000
positions_leapfrog = [position]
positions_tajima = [np.array([1.0, 1.0, 1.0])]
position_boris = [np.array([0.5, 0.5, 0.5])]

# Caculate initial half-step velocity
velocity =  velocity + 0.5 * dt * lorentz_force(charge, velocity, B_field) / mass
velocities_leapfrog = [velocity]
velocities_tajima = [velocity]
velocities_boris = [velocity]

# Calculate initial kinetic energy
E_kins_leapfrog = [0.5 * mass * np.linalg.norm(velocity)**2]
E_kins_tajima = [0.5 * mass * np.linalg.norm(velocity)**2]
E_kins_boris = [0.5 * mass * np.linalg.norm(velocity)**2]

# Leapfrog simulation
position = positions_leapfrog[0]
velocity = velocities_leapfrog[0]
for _ in range(num_steps):
    position, velocity = leapfrog_step(position, velocity, mass, charge, B_field, dt)
    positions_leapfrog.append(position)
    velocities_leapfrog.append(velocity)
    E_kins_leapfrog.append(0.5 * mass * np.linalg.norm(velocity)**2)

# Tajima's implicit method simulation
position = positions_tajima[0]
velocity = velocities_tajima[0]
for _ in range(num_steps):
    position, velocity = tajimasImplicitMethod_step(position, velocity, mass, charge, B_array, dt)
    positions_tajima.append(position)
    velocities_tajima.append(velocity)
    E_kins_tajima.append(0.5 * mass * np.linalg.norm(velocity)**2)

# Boris simulation
position = position_boris[0]
velocity = velocities_boris[0]
for _ in range(num_steps):
    position, velocity = boris_step(position, velocity, mass, charge, B_field, dt)
    position_boris.append(position)
    velocities_boris.append(velocity)
    E_kins_boris.append(0.5 * mass * np.linalg.norm(velocity)**2)

# Convert positions and velocities to numpy arrays for easier manipulation
positions_leapfrog = np.array(positions_leapfrog)
velocities_leapfrog = np.array(velocities_leapfrog)

positions_tajima = np.array(positions_tajima)
velocities_tajima = np.array(velocities_tajima)

position_boris = np.array(position_boris)
velocities_boris = np.array(velocities_boris)

# Save the data to a CSV file
fname = f"/home/anna/semesterproject/anekslen_ippl/build/semesterproject/data/Particles_python_manager.csv"

# Open the file in append mode
with open(fname, mode='a', newline='') as csvfile:
    csvout = csv.writer(csvfile)
    
    # Set precision for floating-point numbers
    csvout.writerows([["Time", "Particle_id", "Position_x", "Position_y", "Position_z", "Cell_id", "Velocity_x", "Velocity_y", "Velocity_z", "E_kin", "E_diff"]])

    # Writing particle data leapfrog
    for i in range(positions_leapfrog.shape[0]):
        position = positions_leapfrog[i]
        velocity = velocities_leapfrog[i]
        particle_id = 0
        cell_id = 0
        ekin = E_kins_leapfrog[i]
        ediff = E_kins_leapfrog[i] - E_kins_leapfrog[0]

        # Write the data row for each particle
        csvout.writerow([i*dt, particle_id, position[0], position[1], position[2], cell_id, velocity[0], velocity[1], velocity[2], ekin, ediff])

    # Writing particle data tajima
    for i in range(positions_tajima.shape[0]):
        position = positions_tajima[i]
        velocity = velocities_tajima[i]
        particle_id = 1
        cell_id = 0
        ekin = E_kins_tajima[i]
        ediff = E_kins_tajima[i] - E_kins_tajima[0]

        # Write the data row for each particle
        csvout.writerow([i*dt, particle_id, position[0], position[1], position[2], cell_id, velocity[0], velocity[1], velocity[2], ekin, ediff])
    
    # Writing particle data boris
    for i in range(position_boris.shape[0]):
        position = position_boris[i]
        velocity = velocities_boris[i]
        particle_id = 2
        cell_id = 0
        ekin = E_kins_boris[i]
        ediff = E_kins_boris[i] - E_kins_boris[0]

        # Write the data row for each particle
        csvout.writerow([i*dt, particle_id, position[0], position[1], position[2], cell_id, velocity[0], velocity[1], velocity[2], ekin, ediff])
