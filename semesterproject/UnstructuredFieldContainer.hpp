#ifndef IPPL_FIELD_CONTAINER_HPP
#define IPPL_FIELD_CONTAINER_HPP

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkGradientFilter.h>
#include <vtkArrayCalculator.h>
#include <vtkCellLocator.h>
#include <vtkPointData.h>
#include <vtkGeometryFilter.h>
#include <vtkPolyDataNormals.h>

#include "datatypes.h"

// TODO: Fix fieldname issue

// Define the FieldsContainer class
template <typename T, unsigned Dim>
class UnstructuredFieldContainer{
public:
    UnstructuredFieldContainer(const char* grid_filename, const char* B_field_name = "B_Field") : B_field_name_m(B_field_name) {
        // Read the grid from the file
        readGrid(grid_filename);
        
        // Print all the arrays in the elements
        vtkSmartPointer<vtkCellData> cellData = grid->GetCellData();
        for (int i = 0; i < cellData->GetNumberOfArrays(); ++i) {
            vtkSmartPointer<vtkDataArray> dataArray = cellData->GetArray(i);
            std::cout << "Array " << i << ": " << dataArray->GetName() << std::endl;
        }
    }

    ~UnstructuredFieldContainer(){}

private:
    const char* B_field_name_m;
    vtkSmartPointer<vtkUnstructuredGrid> grid;

    void readGrid(const char* grid_filename) {
        // Create a reader .vtk file
        vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
        reader->SetFileName(grid_filename);
        reader->Update();
        
        // Get the unstructured grid from the reader
        grid = reader->GetOutput();
        assert(grid);
        
        // Check that the grid contains points, cells, and a B vector field
        assert(grid->GetPoints());
        assert(grid->GetPointData()->HasArray(B_field_name_m));
        assert(grid->GetPointData()->GetArray(B_field_name_m)->GetNumberOfComponents() == Dim);
        assert(grid->GetCells());
    }

public:

    // Write the grid to a file
    void writeGrid(const char* grid_filename) {
        // Check that the magnitude field exists in the grid and is a scalar field, if required
        vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
        writer->SetFileName(grid_filename);
        writer->SetInputData(grid);
        writer->Write();
    }

    // Write a vector field to a file
    void writeField(const char* field_filename, const char* field_name = "B_Field", bool Magnitude = true, const char* magnitude_field_name = "Magnitude") {
        
        // Check that the field exists in the grid and is a vector field
        assert(grid->GetPointData()->HasArray(field_name));
        assert(grid->GetPointData()->GetArray(field_name)->GetNumberOfComponents() == Dim);

        assert(Magnitude && grid->GetPointData()->HasArray(magnitude_field_name));
        assert(Magnitude && grid->GetPointData()->GetArray(magnitude_field_name)->GetNumberOfComponents() == 1);

        // Write the field values to csv file
        std::ofstream file(field_filename);

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
    void calculateCurl(const char* field_name = "B_Field", const char* output_field_name = "Vorticity") {
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
    void calculateMagnitude(const char* field_name = "B_Field", const char* output_field_name = "Magnitude") {
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

    void getGridBounds(Vector_t<double, Dim> &min, Vector_t<double, Dim> &max) {
        double bounds[6];
        grid->GetBounds(bounds);
        for (unsigned i = 0; i < Dim; i++) {
            min[i] = bounds[i * 2];
            max[i] = bounds[i * 2 + 1];
        }
    }

    // TODO: this is now done for a grid in cylindric coordinates where we have 1/8 of the xy cyrcle and the upper z part of the grid, generalize this possibly
    Vector_t<double, Dim> getReferencePosition(Vector_t<double, Dim> R) {
        Vector_t<double, Dim> R_ref = ippl::fabs(R);

        if(R_ref[0] <= R_ref[1]) {
            R_ref[0] = std::abs(R[1]);
            R_ref[1] = std::abs(R[0]);
        }
        
        return R_ref;

    }

    Vector_t<double, Dim> getRealB_field(Vector_t<double, Dim> B, Vector_t<double, Dim> R) {
        
        Vector_t<double, Dim> B_real;

        if(std::abs(R[0]) <= std::abs(R[1])) {
            B_real[0] = B[1];
            B_real[1] = B[0];
            B_real[2] = B[2];
        }
        else {
            B_real = B;
        }

        if(R[0] < 0) {
            B_real[0] = -B_real[0];
        }

        if(R[1] < 0) {
            B_real[1] = -B_real[1];
        }

        if(R[2] < 0) {
            B_real[0] = -B_real[0];
            B_real[1] = -B_real[1];
        }

        return B_real;
    }


    // Return the Id of the grid cell in which a point is located
    int GetGridCell(Vector_t<double, Dim> R, double (&weights)[8]) {
        // Get reference position
        Vector_t<double, Dim> R_ref = getReferencePosition(R);

        // Create and build a cell locator
        vtkSmartPointer<vtkCellLocator> locator = vtkSmartPointer<vtkCellLocator>::New();
        locator->SetDataSet(grid);
        locator->BuildLocator();

        // Find the cell that contains the point (x, y, z)
        double point[Dim] = { R_ref[0], R_ref[1], R_ref[2] };
        double tol2 = 0.0;
        vtkSmartPointer<vtkGenericCell> GenCell = vtkSmartPointer<vtkGenericCell>::New();
        double pcoords[Dim];

        double cellId = locator->FindCell(point, tol2, GenCell, pcoords, weights);
        
        return cellId;
    }

    int GetGridCell(Vector_t<double, Dim> R) {
        double weights[8];
        return GetGridCell(R, weights);
    }

    vtkIdType FindCellAndInterpolateField(Vector_t<double, Dim> R, Vector_t<double, Dim> &interpolatedField, const char* field_name = "B_Field") {
        // Check that the field exists in the grid and is a vector field
        assert(grid->GetPointData()->HasArray(field_name));
        assert(grid->GetPointData()->GetArray(field_name)->GetNumberOfComponents() == Dim);

        // Set the interpolated field to zero
        Vector_t<double, Dim> B_field_ref = Vector_t<double, Dim>(0.0);

        // Get gridcell and interpolation weights
        double weights[8];
        vtkIdType cellId = GetGridCell(R, weights);
        // Check that the point is inside the grid
        if(cellId == -1) {
            return cellId;
        }

        // Get the cell that contains the point
        vtkSmartPointer<vtkCell> cell = grid->GetCell(cellId);

        // Access the field
        vtkSmartPointer<vtkDataArray> fieldArray = grid->GetPointData()->GetArray(B_field_name_m);

        double B_val[Dim];        
        for(unsigned i = 0; i < cell->GetNumberOfPoints(); i++) {
            vtkIdType pointId = cell->GetPointId(i);
            fieldArray->GetTuple(pointId, B_val);
            for(unsigned j = 0; j < Dim; j++) {
                B_field_ref[j] += weights[i] * B_val[j];
            }
        }

        // Get the real field
        interpolatedField = getRealB_field(B_field_ref, R);

        return cellId;
    }

    void interpolateField(Vector_t<double, Dim> R, Vector_t<double, Dim> &interpolatedField, const char* field_name = "B_Field")
    {
        FindCellAndInterpolateField(R, interpolatedField, field_name);
    }

    void getBoundaryInformation(vtkIdType cellID, std::string &boundaryType) {
        vtkSmartPointer<vtkDataArray> boundaryArray = grid->GetCellData()->GetArray("Boundary_Type");
        if (!boundaryArray) {
            std::cerr << "Error: Boundary_Type array not found!" << std::endl;
            boundaryType = "Unknown";
            return;
        }
    
        if (cellID < 0 || cellID >= boundaryArray->GetNumberOfTuples()) {
            std::cerr << "Error: " << cellID << "is an invalid cellID!" << std::endl;
            boundaryType = "Unknown";
            return;
        }
    
        double boundaryValue = boundaryArray->GetTuple1(cellID);
        boundaryType = std::to_string(static_cast<int>(boundaryValue));
    }
};

#endif