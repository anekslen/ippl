import vtk
from collections import defaultdict
import csv
import sys
import os

# Read data from CSV file
csv_file = sys.argv[1] if len(sys.argv) > 1 else '/home/annah/semesterproject/code/ippl/build/semesterproject/data/Particles_1_manager.csv'
particles = defaultdict(list)

with open(csv_file, 'r') as file:
    reader = csv.DictReader(file)
    headers = next(reader)  # Skip the header row
    for row in reader:
        point_id = int(row['Particle_id'])
        position = (float(row['Position_x'].strip()), float(row['Position_y'].strip()), float(row['Position_z'].strip()))
        velocity = (float(row['Velocity_x'].strip()), float(row['Velocity_y'].strip()), float(row['Velocity_z'].strip()))
        particles[point_id].append((position, velocity))

print(f"Read positions for {len(particles)} particles from CSV file")

# Create a vtkMultiBlockDataSet to store each particle's path as a separate block
multi_block = vtk.vtkMultiBlockDataSet()
block_index = 0

for point_id, data in particles.items():
    polydata = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    velocities = vtk.vtkDoubleArray()
    velocities.SetNumberOfComponents(3)
    velocities.SetName('Velocity')
    cells = vtk.vtkCellArray()

    # Store positions in vtkPoints
    start_index = points.GetNumberOfPoints()  # Keep track of the starting index for this particle's points
    for pos, vel in data:
        points.InsertNextPoint(pos)
        velocities.InsertNextTuple(vel)

    # Create vtkPolyLine for this particle's path
    polyline = vtk.vtkPolyLine()
    polyline.GetPointIds().SetNumberOfIds(len(data))
    for i in range(len(data)):
        polyline.GetPointIds().SetId(i, start_index + i)

    # Add polyline to cells
    cells.InsertNextCell(polyline)

    # Store all points and lines in vtkPolyData
    polydata.SetPoints(points)
    polydata.SetLines(cells)
    polydata.GetPointData().SetVectors(velocities)

    # Add the polydata to the multi-block dataset
    multi_block.SetBlock(block_index, polydata)
    block_index += 1

# Write the multi-block dataset to a VTK file
writer = vtk.vtkXMLMultiBlockDataWriter()
output_file = os.path.join(os.path.dirname(csv_file), 'particles.vtm')
writer.SetFileName(output_file)
writer.SetInputData(multi_block)
writer.Write()

print("Particle paths written to a single VTM file")
