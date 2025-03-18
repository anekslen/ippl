import numpy as np
import vtk
import os

version = input("Test mesh? (t) or Reactor mesh old? (rold) or Reactor mesh new? (rnew): ")

while(version != 't' and version != 'rold' and version != 'rnew'):
    print("Invalid input. Please enter 't' for the Test mesh or 'r' for the Reactor mesh.")
    version = input("Test mesh? (t) or Reactor mesh? (r): ")

boundary_files = {}

if version == 't':
    nodefile = 'testmesh/vtkInput/nodes.txt'
    elementfile = 'testmesh/vtkInput/mesh.txt'
    fieldfile = 'testmesh/vtkInput/data.txt'
    resultfile = 'testmesh/mesh.vtk'
    mapingfile = 'testmesh/ElementMapping.csv'
elif version == 'rold':
    nodefile = 'reactormesh/vtkInput/OldMesh/Node22642_I8_3G20_Rev.txt'
    elementfile = 'reactormesh/vtkInput/OldMesh/Elem19482_8I6_I6_24X_I6.txt'
    fieldfile = 'reactormesh/vtkInput/OldMesh/B22642_XYZmod_I8_1X_4E12.txt'
    resultfile = 'reactormesh/mesh.vtk'
    mapingfile = 'reactormesh/ElementMapping.csv'

elif version == 'rnew':
    nodefile = 'reactormesh/vtkInput/NewMesh/Plasm18592N.txt'
    elementfile = 'reactormesh/vtkInput/NewMesh/PlasmE15641.txt'
    fieldfile = 'reactormesh/vtkInput/NewMesh/Plasma_B.txt' 

    boundary_folder = 'reactormesh/vtkInput/NewMesh/BoundaryClassification/'
    for filename in os.listdir(boundary_folder):
        if filename.endswith('.txt'):
            filepath = os.path.join(boundary_folder, filename)
            # Determine the number of columns in the file
            with open(filepath, 'r') as f:
                first_line = f.readline()
                num_columns = len(first_line.split())
            # Create the converters dictionary
            converters = {i: int for i in range(num_columns)}
            boundary_files[filename] = np.loadtxt(filepath, converters=converters)
    resultfile = 'reactormesh/mesh.vtk'
    mapingfile = 'reactormesh/ElementMapping.csv'


# Create a mapping from original IDs to VTK IDs
point_id_map = {}
element_id_map = {}

# Load the node coordinates
nodes = np.loadtxt(nodefile, dtype=float)

# Load the element connectivity
elements = np.loadtxt(elementfile, dtype=int)

# Load the B field data
b_field_data = np.loadtxt(fieldfile, dtype=float)

boundary_data = {}
boundary_sets = {}

if version == 'rnew':
    for key, value in boundary_files.items():
        boundarytype = key.split('.')[0]
        boundary_data[boundarytype] = np.loadtxt(os.path.join(boundary_folder, key), dtype=int)
        boundary_sets[boundarytype] = set(tuple(row[:]) for row in boundary_data[boundarytype])

    print("Files loades")

# Create a VTK Points object to store the coordinates
points = vtk.vtkPoints()

# Create an array to store custom point IDs (Node IDs)
point_ids = vtk.vtkIntArray()
point_ids.SetName("Node_ID")  # Name the array

# Create an array to store the B field
b_field = vtk.vtkFloatArray()
b_field.SetName("B_Field")  # Name the array
b_field.SetNumberOfComponents(3)  # 3 components (Bx, By, Bz)

for i in range(nodes.shape[0]):
    points.InsertNextPoint(nodes[i][1], nodes[i][2], nodes[i][3])
    point_ids.InsertNextValue(int(nodes[i][0]))  # Insert custom point ID (Node_ID)
    point_id_map[int(nodes[i][0])] = i

    assert(int(nodes[i][0]) == int(b_field_data[i][0]))  # Ensure the node IDs match
    b_field.InsertNextTuple3(b_field_data[i][1], b_field_data[i][2], b_field_data[i][3])  # Insert B field values

# Create a VTK Unstructured Grid
ugrid = vtk.vtkUnstructuredGrid()
ugrid.SetPoints(points)

# Add the custom point IDs to the grid's point data
ugrid.GetPointData().AddArray(point_ids)

# Add the B field data to the grid's point data
ugrid.GetPointData().AddArray(b_field)

# Initialize five scalar arrays for element values
cell_values = []
for i in range(1, 6):
    array = vtk.vtkFloatArray()
    array.SetName(f"Value{i}")
    cell_values.append(array)

array = vtk.vtkFloatArray()
array.SetName(f"Element_Type")
cell_values.append(array)

array = vtk.vtkFloatArray()
array.SetName(f"Element_Number")
cell_values.append(array)

array = vtk.vtkStringArray()
array.SetName(f"Boundary_Type")
cell_values.append(array)

element_id = 0

num_boundary_types = len(boundary_sets)
boundary_number = np.zeros(num_boundary_types + 1, dtype=int)
counter = np.zeros(num_boundary_types + 1, dtype=int)

for i, (key, value) in enumerate(boundary_sets.items()):
    boundary_number[i+1] = len(value)

boundary_number[0] = elements.shape[0] - sum(boundary_number[1:])

# Iterate over elements and insert them into the unstructured grid
for element in elements:
    if element[8] == 1:  # If the element is in the vaccuum mesh (VTK_LINE)
        type = 0
    else:
        type = 1

    if element[2] == element[3] and element[6] == element[7]:  # If the element has 6 nodes (VTK_WEDGE)
        wedge = vtk.vtkWedge()
        for i in range(6):
            wedge.GetPointIds().SetId(i, point_id_map[element[i+int(i/3)]])
        ugrid.InsertNextCell(wedge.GetCellType(), wedge.GetPointIds())

    else:  # If the element has 8 nodes (VTK_HEXAHEDRON)
        hexahedron = vtk.vtkHexahedron()
        for i in range(8):
            hexahedron.GetPointIds().SetId(i, point_id_map[element[i]])
        ugrid.InsertNextCell(hexahedron.GetCellType(), hexahedron.GetPointIds())

    # Insert the values for the current element
    for i in range(5):
        cell_values[i].InsertNextValue(element[i+8])  # Value1 to Value5
    cell_values[5].InsertNextValue(type)  # Element_Type
    cell_values[6].InsertNextValue(element[13]) #Element_Number

    # Check if element is in any boundary set and insert the additional value
    for i, (boundarytype, value) in enumerate(boundary_sets.items()):
        if tuple(element[:]) in value:
            cell_values[7].InsertNextValue(boundarytype)
            counter[i+1] += 1
            break
    else:
        cell_values[7].InsertNextValue("Interior")
        counter[0] += 1

    element_id_map[element[13]] = element_id
    element_id += 1

# Add the cell values array to the unstructured grid's cell data
for vals in cell_values:
    ugrid.GetCellData().AddArray(vals)

# Write the grid to a VTK file
writer = vtk.vtkUnstructuredGridWriter()
writer.SetFileName(resultfile)
writer.SetInputData(ugrid)
writer.Write()

print("Mesh written to 'mesh.vtk'")
print(f"Number of elements in the boundary types: {counter}")
print(f"Number of elements should boundary types: {boundary_number}")
for(i, (key, value)) in enumerate(boundary_sets.items()):
    print(f"Boundary type {i + 1} has name {key}")

# Create file that matches vtk node numbers to input file node number
with open(mapingfile, 'w') as f:
    for key in element_id_map:
        f.write(f"{key}, {element_id_map[key]}\n")
    f.close()