import numpy as np
import matplotlib.pyplot as plt

# Calculate coordinates for the boundary lines and circle segments for inner and outer boundary regions
# Define inner and outer boundary regions (inner [0], outer [1])
RMin = [0.005, 0.318]
RMax = [0.14, 0.332]

ThetaMin = [0, 20]
ThetaMin = [value / 360 * 2 * np.pi for value in ThetaMin]

ThetaMax = [45, 45]
ThetaMax = [value / 360 * 2 * np.pi for value in ThetaMax]

# Calculate Points for circle segments and lines of input particle area
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
theta = [np.linspace(ThetaMin[0], ThetaMax[0], 100), np.linspace(ThetaMin[1], ThetaMax[1], 100)]  # Circle segments

X_RMin_circle = [RMin[0] * np.cos(theta[0]), RMin[1] * np.cos(theta[1])]
y_RMin_circle = [RMin[0] * np.sin(theta[0]), RMin[1] * np.sin(theta[1])]

X_RMax_circle = [RMax[0] * np.cos(theta[0]), RMax[1] * np.cos(theta[1])]
y_RMax_circle = [RMax[0] * np.sin(theta[0]), RMax[1] * np.sin(theta[1])]

# Calculate circle segments for the mesh region
R = [0, 0.158, 0.292, 0.36]
Theta = [0, 45]
Theta = [value / 360 * 2 * np.pi for value in Theta]

# Calculate Points for circle segments
Mesh_X_ThetaMin = R[-1] * np.cos(Theta[0])
Mesh_Y_ThetaMin = R[-1] * np.sin(Theta[0])

Mesh_X_ThetaMax = R[-1] * np.cos(Theta[1])
Mesh_Y_ThetaMax = R[-1] * np.sin(Theta[1])

# Calculate boundary lines of mesh region
lines_mesh = [
    [(0,0), (Mesh_X_ThetaMin, Mesh_Y_ThetaMin)],
    [(0,0), (Mesh_X_ThetaMax, Mesh_Y_ThetaMax)]
]

# Calculate circle segments of mesh region
theta_mesh = np.linspace(Theta[0], Theta[1], 100)  # Circle segments

Mesh_X_circles = [R[i] * np.cos(theta_mesh) for i in range(len(R))]
Mesh_Y_circles = [R[i] * np.sin(theta_mesh) for i in range(len(R))]

# Calculate the rectangle areas for the mesh region
# Calculate the rectangle areas for the mesh region
squares = []
for i in range(len(X_RMin_circle[0]) - 1):  # Loop through the circle segments
    square_inner = [
        (X_RMin_circle[0][i], y_RMin_circle[0][i]),
        (X_RMin_circle[0][i + 1], y_RMin_circle[0][i + 1]),
        (X_RMax_circle[0][i + 1], y_RMax_circle[0][i + 1]),
        (X_RMax_circle[0][i], y_RMax_circle[0][i])
    ]
    squares.append(square_inner)

for i in range(len(X_RMin_circle[1]) - 1):  # Loop through the circle segments
    square_outer = [
        (X_RMin_circle[1][i], y_RMin_circle[1][i]),
        (X_RMin_circle[1][i + 1], y_RMin_circle[1][i + 1]),
        (X_RMax_circle[1][i + 1], y_RMax_circle[1][i + 1]),
        (X_RMax_circle[1][i], y_RMax_circle[1][i])
    ]
    squares.append(square_outer)


# Plot the boundary lines and circle segments
fig, ax = plt.subplots()

for line in lines_mesh:
    ax.plot([line[0][0], line[1][0]], [line[0][1], line[1][1]], color='black', linewidth=1)

for i in range(len(Mesh_X_circles)):
    ax.plot(Mesh_X_circles[i], Mesh_Y_circles[i], color='black', linewidth=1)

for square in squares:
    square = np.array(square)
    ax.fill(square[:, 0], square[:, 1], color=[0/255, 73/255, 104/255], alpha=1.0)

# Remove x and y axes
ax.axis('off')


ax.set_aspect('equal', 'box')
plt.savefig("ParticleRegions.png")




# Plot a circle and lines splitting it into octants with numbers in each octant
fig, ax = plt.subplots()

# Circle
circle = plt.Circle((0, 0), 1, color='black', fill=False, linewidth=1)
ax.add_artist(circle)

# Lines splitting the circle into octants
angles = np.linspace(0, 2 * np.pi, 9)[:-1]  # 8 lines (octants)
for angle in angles:
    ax.plot([0, np.cos(angle)], [0, np.sin(angle)], color='black', linewidth=0.5)

# Add numbers in each octant and fill octant 2 with color
radius = 0.7  # Distance from the center for the numbers
for i, angle in enumerate(angles):
    if i == 0:  # Octant 1 (index 0)
        # Fill the entire circle segment for octant 2
        theta = np.linspace(angle, angles[i + 1], 100)
        x_segment = np.append([0], np.cos(theta))
        y_segment = np.append([0], np.sin(theta))
        ax.fill(x_segment, y_segment, color=[0/255, 73/255, 104/255], alpha=0.7)
    # Add the number
    x = radius * np.cos(angle + np.pi / 8)  # Offset to center the number in the octant
    y = radius * np.sin(angle + np.pi / 8)
    ax.text(x, y, f"{i + 1}", fontsize=12, ha='center', va='center')
    circle = plt.Circle((x, y), 0.1, color='black', fill=False, linewidth=0.5)
    ax.add_artist(circle)

# Add x and y axes with arrows on both sides
ax.arrow(-1.2, 0, 2.4, 0, head_width=0.03, head_length=0.07, fc='black', ec='black', linewidth=0.5, length_includes_head=True)
ax.arrow(1.2, 0, -2.4, 0, head_width=0.03, head_length=0.07, fc='black', ec='black', linewidth=0.5, length_includes_head=True)
ax.arrow(0, -1.2, 0, 2.4, head_width=0.03, head_length=0.07, fc='black', ec='black', linewidth=0.25, length_includes_head=True)
ax.arrow(0, 1.2, 0, -2.4, head_width=0.03, head_length=0.07, fc='black', ec='black', linewidth=0.25, length_includes_head=True)

# Add (0,0) label and 0s to the sides
ax.text(1.3, 0, "0.0", fontsize=10, ha='center', va='center')  # Move left 0 slightly down
ax.text(0, 1.3, "0.0", fontsize=10, ha='center', va='center')  # Move bottom 0 slightly left

# Add mini x y axis symbol lower right
ax.arrow(1.3, -1, 0.3, 0, head_width=0.03, head_length=0.07, fc='black', ec='black', linewidth=0.5, length_includes_head=True)
ax.arrow(1.3, -1, 0, 0.3, head_width=0.03, head_length=0.07, fc='black', ec='black', linewidth=0.5, length_includes_head=True)
ax.text(1.7, -1, "x", fontsize=10, ha='center', va='center')
ax.text(1.3, -0.6, "y", fontsize=10, ha='center', va='center')

# Remove x and y ticks and the square around the plot
ax.axis('off')

ax.set_aspect('equal', 'box')
plt.savefig("CircleOctantsWithNumbers.png")
