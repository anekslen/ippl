import numpy as np
import vtk

version = input("Test mesh? (t) or Reactor mesh? (r): ")

while(version != 't' and version != 'r'):
    print("Invalid input. Please enter 't' for the Test mesh or 'r' for the Reactor mesh.")
    version = input("Test mesh? (t) or Reactor mesh? (r): ")

if version == 't':
    nodefile = 'testmesh/vtkInput/nodes.txt'
    elementfile = 'testmesh/vtkInput/mesh.txt'
    fieldfile = 'testmesh/vtkInput/data.txt'
    resultfile = 'mesh.vtk'
elif version == 'r':
    nodefile = 'reactormesh/vtkInput/Node22642_I8_3G20_Rev.txt'
    elementfile = 'reactormesh/vtkInput/Elem19482_8I6_I6_24X_I6.txt'
    fieldfile = 'reactormesh/vtkInput/B22642_XYZmod_I8_1X_4E12.txt'
    resultfile = 'mesh.vtk'


# Create a mapping from original IDs to VTK IDs
point_id_map = {}
element_id_map = {}

# Load the node coordinates
nodes = np.loadtxt(nodefile, dtype=float)

# Load the element connectivity
elements = np.loadtxt(elementfile, dtype=int)

# Load the B field data
b_field_data = np.loadtxt(fieldfile, dtype=float)

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

element_id = 0

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
        type = 1
        hexahedron = vtk.vtkHexahedron()
        for i in range(8):
            hexahedron.GetPointIds().SetId(i, point_id_map[element[i]])
        ugrid.InsertNextCell(hexahedron.GetCellType(), hexahedron.GetPointIds())

    # Insert the values for the current element
    for i in range(5):
        cell_values[i].InsertNextValue(element[i+8])  # Value1 to Value5
    cell_values[5].InsertNextValue(type)  # Element_Type
    cell_values[6].InsertNextValue(element[13]) #Element_Number

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