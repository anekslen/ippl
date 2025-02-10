import random
import numpy as np

# Create the data for the nodes
def create_nodedata(nodenumbers, filename, delimiter=' ', origin = [0, 0, 0], domain_length = [1, 1, 1]):
    # nodenumbers is the number of nodes in each dimension

    n = 0  # Initialize the number counter
    with open(filename, 'w') as file:
        # Write the header
        #file.write(f"node{delimiter}x{delimiter}y{delimiter}z\n")

        h = [domain_length[i] / (nodenumbers[i] - 1) for i in range(3)]
        
        # Write the data lines
        for z in [origin[2] + i * h[2] for i in range(nodenumbers[2])]:
            for y in [origin[1] + i * h[1] for i in range(nodenumbers[1])]:
                for x in [origin[0] + i * h[0] for i in range(nodenumbers[0])]:
                    file.write(f"{n}{delimiter}{x}{delimiter}{y}{delimiter}{z}\n")
                    n += 1  # Increment the counter

    assert(n == nodenumbers[0] * nodenumbers[1] * nodenumbers[2])
    print("Node data created successfully")

# Create the data for the mesh elements
def create_meshdata(nodenumbers, filename, delimiter=' '):
    # nodenumbers is the number of nodes in each dimension

    with open(filename, 'w') as file:
        # Write the header
        #file.write(f"element{delimiter}n1{delimiter}n2{delimiter}n3{delimiter}n4{delimiter}n5{delimiter}n6{delimiter}n7{delimiter}n8{delimiter}meshtype{delimiter}rand1{delimiter}rand2{delimiter}rand3{delimiter}rand4{delimiter}elementnumber\n")
        
        # number of elements in each dimension
        elemNum = [num - 1 for num in nodenumbers]

        # Write the data lines
        for z in range(elemNum[2]):
            for y in range(elemNum[1]):
                for x in range(elemNum[0]):
                    elem = x + y * elemNum[1] + z * elemNum[1] * elemNum[2]

                    n1 = x + y * nodenumbers[1] + z * nodenumbers[1] * nodenumbers[2]
                    n2 = n1 + 1
                    n3 = n2 + nodenumbers[1]
                    n4 = n3 - 1
                    n5 = n1 + nodenumbers[1] * nodenumbers[2] 
                    n6 = n2 + nodenumbers[1] * nodenumbers[2] 
                    n7 = n3 + nodenumbers[1] * nodenumbers[2] 
                    n8 = n4 + nodenumbers[1] * nodenumbers[2] 

                    meshtype = 0

                    rand1 = 1
                    rand2 = 2
                    rand3 = 3
                    rand4 = 4

                    file.write(f"{n1}{delimiter}{n2}{delimiter}{n3}{delimiter}{n4}{delimiter}{n5}{delimiter}{n6}{delimiter}{n7}{delimiter}{n8}{delimiter}{meshtype}{delimiter}{rand1}{delimiter}{rand2}{delimiter}{rand3}{delimiter}{rand4}{delimiter}{elem}\n")
    print("Mesh data created successfully")

# Create random data for the nodes
def create_nodevalues(nodenumbers, B_field, filename, delimiter=' '):
    # nodenumbers is the number of nodes in each dimension

    with open(filename, 'w') as file:
        # Write the header
        #file.write(f"node{delimiter}Bx{delimiter}By{delimiter}Bz\n")
        
        # Write the data lines
        for i in range(nodenumbers[0] * nodenumbers[1] * nodenumbers[2]):
            file.write(f"{i}{delimiter}{B_field[0]}{delimiter}{B_field[1]}{delimiter}{B_field[2]}\n")
    print("Node values created successfully")

# Create the data for points inside the elements for checking the B-field interpolation
def create_pointdata(nodenumbers, B_field, filename, delimiter=' '):
    # nodenumbers is the number of nodes in each dimension

    meshnodes = nodenumbers[0] * nodenumbers[1] * nodenumbers[2]

    n = 0  # Initialize the number counter
    with open(filename, 'w') as file:
        # Write the header
        #file.write(f"node{delimiter}x{delimiter}y{delimiter}z\n")

        h = [domain_length[i] / (nodenumbers[i] - 1) for i in range(3)]
        
        # Write the data lines, creates a point at the center of each element
        for i in range(nodenumbers[2] - 1):
            for j in range(nodenumbers[1] - 1):
                for k in range(nodenumbers[0] - 1):
                    x = origin[0] + h[0] / 2 + k * h[0]
                    y = origin[1] + h[1] / 2 + j * h[1]
                    z = origin[2] + h[2] / 2 + i * h[2]

                    elem = k + j * (nodenumbers[1] - 1) + i * (nodenumbers[1] - 1) * (nodenumbers[2] - 1)

                    file.write(f"{meshnodes + n}{delimiter}{x}{delimiter}{y}{delimiter}{z}{delimiter}{elem}{delimiter}{B_field[0]}{delimiter}{B_field[1]}{delimiter}{B_field[2]}\n")
                    n += 1  # Increment the counter
    print("Point data created successfully")

# Read the mesh parameters from a file
def read_mesh_parameters(filename):
    data = np.loadtxt(filename, delimiter=',', dtype=str)
    origin = [float(num) for num in (data[1][k] for k in range(1, 4))]
    num = [int(num) for num in (data[2][k] for k in range(1, 4))]
    domain_length = [float(num) for num in (data[3][k] for k in range(1, 4))]
    B_field = [float(num) for num in (data[4][k] for k in range(1, 4))]

    print("B_field: ", B_field)

    print("Mesh parameters read successfully")

    return origin, num, domain_length, B_field

origin, num, domain_length, B_field = read_mesh_parameters('input/meshdata.txt')

create_nodedata(num, 'vtkInput/nodes.txt', origin = origin, domain_length = domain_length)
create_meshdata(num, 'vtkInput/mesh.txt')
create_nodevalues(num, B_field, 'vtkInput/data.txt')
create_pointdata(num, B_field, 'Pointdata/point_data.txt')