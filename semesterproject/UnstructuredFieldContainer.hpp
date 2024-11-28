#ifndef IPPL_FIELD_CONTAINER_H
#define IPPL_FIELD_CONTAINER_H

// Define the FieldsContainer class
template <typename T, unsigned Dim = 3>
class UnstructuredFieldContainer{
public:
    UnstructuredFieldContainer(const char* filename, const char* B_field_name = "B_Field") : B_field_name_m(B_field_name) {
        // Read the grid from the file
        
        // Create a reader .vtk file
        vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
        reader->SetFileName(filename);
        reader->Update();
        
        // Get the unstructured grid from the reader
        grid = reader->GetOutput();
        assert(grid);
        
        // Check that the grid contains points, cells, and a B vector field
        assert(grid->GetPoints());
        assert(grid->GetPointData()->HasArray(field_name));
        assert(grid->GetPointData()->GetArray(field_name)->GetNumberOfComponents() == Dim);
        assert(grid->GetCells());
    }

    ~UnstructuredFieldContainer(){}

private:
    const char* B_field_name_m;
    vtkSmartPointer<vtkUnstructuredGrid> grid;

public:

    // Write the grid to a file
    void writeGrid(const char* filename) {
        // Check that the magnitude field exists in the grid and is a scalar field, if required
        vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
        writer->SetFileName(filename);
        writer->SetInputData(grid);
        writer->Write();
    }

    // Write a vector field to a file
    void writeField(const char* filename, const char* field_name, bool Magnitude = false, const char* magnitude_field_name = "Magnitude") {
        
        // Check that the field exists in the grid and is a vector field
        assert(grid->GetPointData()->HasArray(field_name));
        assert(grid->GetPointData()->GetArray(field_name)->GetNumberOfComponents() == Dim);

        assert(Magnitude && grid->GetPointData()->HasArray(magnitude_field_name));
        assert(Magnitude && grid->GetPointData()->GetArray(magnitude_field_name)->GetNumberOfComponents() == 1);

        // Write the field values to csv file
        std::ofstream file(filename);

        if (Magnitude) {
            file << "Point_Index," << field_name << "_X," << field_name << "_Y," << field_name << "_Z," << magnitude_field_name << std::endl;
        }
        else {
            file << "Point_Index," << field_name << "_X," << field_name << "_Y," << field_name << "_Z" << std::endl;
        }

        // Write the field values to csv file
        for (int i = 0; i < grid->GetPoints()->GetNumberOfPoints(); i++) {
            double val[Dim];
            double magnitude;
            grid->GetPointData()->GetArray(field_name)->GetTuple(i, val);
            if (Magnitude) {
                magnitude = grid->GetPointData()->GetArray(magnitude_field_name)->GetTuple1(i);
            }

            // Set the format for floating-point values
            file << std::scientific << std::setprecision(4) << std::setw(10);

            if (Magnitude) {
                file << i << ","
                    << std::setw(10) << val[0] << ","
                    << std::setw(10) << val[1] << ","
                    << std::setw(10) << val[2] << ","
                    << std::setw(10) << magnitude << std::endl;
            }
            else {
                file << i << ","
                    << std::setw(10) << val[0] << ","
                    << std::setw(10) << val[1] << ","
                    << std::setw(10) << val[2] << std::endl;
            }
        }
        file.close();
    }

    // Calculate the curl of the vector field field_name and save it in the output_field_name
    void calculateCurl(const char* field_name, const char* output_field_name = "Vorticity") {
        // Check that the field exists in the grid and is a vector field
        assert(grid->GetPointData()->HasArray(field_name));
        assert(grid->GetPointData()->GetArray(field_name)->GetNumberOfComponents() == Dim);

        // Check that the output field does not exist in the grid
        assert(!grid->GetPointData()->HasArray(output_field_name));

        // Calculate the curl of the vector field at each point
        vtkSmartPointer<vtkGradientFilter> gradientFilter = vtkSmartPointer<vtkGradientFilter>::New();
        gradientFilter->SetInputData(grid);
        gradientFilter->SetInputScalars(vtkDataObject::FIELD_ASSOCIATION_POINTS, field_name);
        gradientFilter->SetComputeGradient(false);   // Disable the gradient computation
        gradientFilter->SetComputeVorticity(true);   // Enable the curl computation
        gradientFilter->Update();

        // Add the curl field to the grid and rename it to the output_field_name
        grid->GetPointData()->AddArray(gradientFilter->GetUnstructuredGridOutput()->GetPointData()->GetArray("Vorticity"));
        grid->GetPointData()->GetArray("Vorticity")->SetName(output_field_name);
    }

    // Define functions to calculate the magnitude of a vector field
    void calculateMagnitude(const char* field_name, const char* output_field_name = "Magnitude") {
        // Check that the field exists in the grid and is a vector field
        assert(grid->GetPointData()->HasArray(field_name));
        assert(grid->GetPointData()->GetArray(field_name)->GetNumberOfComponents() == Dim);

        // Check that the output field does not exist in the grid
        assert(!grid->GetPointData()->HasArray(output_field_name));

        // Calculate the magnitude of the curl at each point
        std::string function = "mag(" + std::string(field_name) + ")";

        vtkSmartPointer<vtkArrayCalculator> magnitudeCalculator = vtkSmartPointer<vtkArrayCalculator>::New();
        magnitudeCalculator->SetInputData(grid);
        magnitudeCalculator->SetAttributeTypeToPointData();
        magnitudeCalculator->AddVectorArrayName(field_name);
        magnitudeCalculator->SetFunction(function.c_str());    // Calculate the magnitude of the curl
        magnitudeCalculator->SetResultArrayName(output_field_name);
        magnitudeCalculator->Update();

        // Add the magnitude of the curl to the curl grid
        grid->GetPointData()->AddArray(magnitudeCalculator->GetUnstructuredGridOutput()->GetPointData()->GetArray(output_field_name));
    }

    void getGridBounds(Vector<double, Dim> &min, Vector<double, Dim> &max) {
        double bounds[6];
        grid->GetBounds(bounds);
        for (int i = 0; i < Dim; i++) {
            min[i] = bounds[i * 2];
            max[i] = bounds[i * 2 + 1];
        }
    }


    // Return the Id of the grid cell in which a point is located
    int GetGridCell(double x, double y, double z) {
        vtkIdType cellId = -1;
        double p[Dim] = { x, y, z };
        double pcoords[Dim];
        double weights[8];
        int subId;
        double dist2 = -1.0;

        // Optimization if previous cell is known (set cellId and cell to previous cell)
        // TODO: Not thread save, check this !!!!!!
        cellId = grid->FindCell(p, NULL, 0, 0, dist2, subId, pcoords, weights);
        return cellId;
    }

    // Interpolate a field at a given point TODO: Finish this function
    void interpolateField(const char* field_name, double x, double y, double z, std::vector<int> &interpolatedField, int cellID = -1, double weights[8] = NULL)
    {
        // Check that the field exists in the grid and is a vector field
        assert(grid->GetPointData()->HasArray(field_name));
        assert(grid->GetPointData()->GetArray(field_name)->GetNumberOfComponents() == Dim); return;

        if (cellID == -1 || weights == NULL) {
            // Find the cell that contains the point (x, y, z)
            double p[Dim] = { x, y, z };
            double pcoords[Dim];
            double weights[8];
            int subId;
            double dist2 = -1.0;
            // Optimization if previous cell is known (set cellId and cell to previous cell)
            // TODO: Not thread save, check this !!!!!!
            cellID = grid->FindCell(p, NULL, 0, 0, dist2, subId, pcoords, weights);

            // Access the field
            vtkDataArray* fieldArray = grid->GetPointData()->GetArray(field_name);
            if (!fieldArray) {
                std::cerr << "Error: Failed to get the " << field_name << " array from the mesh!" << std::endl;
                return;
            }

            // Get the point indices of the cell in which the point (x, y, z) is located
            vtkIdList* pointIds = grid->GetCell(cellID)->GetPointIds();
            
            // Interpolate the field at the point (x, y, z)
            vtkIdType interpolationValues[Dim];

            fieldArray->InterpolateTuple(interpolationValues, pointIds, fieldArray, weights);

            for(int i = 0; i < Dim; i++) {
                interpolatedField[i] = interpolationValues[i];
            }
            
            return;
        }
    }
};

#endif