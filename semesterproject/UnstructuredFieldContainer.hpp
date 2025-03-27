#ifndef IPPL_FIELD_CONTAINER_HPP
#define IPPL_FIELD_CONTAINER_HPP

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkGradientFilter.h>
#include <vtkArrayCalculator.h>
#include <vtkStaticCellLocator.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkGeometryFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkAbstractArray.h>
#include <vtkStringArray.h>

#include "datatypes.h"

// TODO: Fix fieldname issue

// Define the FieldsContainer class
template <typename T, unsigned Dim>
class UnstructuredFieldContainer{
public:
    UnstructuredFieldContainer(const char* grid_filename, const char* B_field_name = "B_Field") : B_field_name_m(B_field_name) {
        // Get the number of threads from Kokkos
        int numThreads = Kokkos::DefaultExecutionSpace::concurrency();
        std::cout << "Number of threads available: " << numThreads << std::endl;
        
        // Read the grid and grid copies for each thread from the file
        grids.resize(numThreads);
        readGrid(grid_filename);
        
        // Print all the arrays in the elements
        vtkSmartPointer<vtkCellData> cellData = grid->GetCellData();
        for (int i = 0; i < cellData->GetNumberOfArrays(); ++i) {
            std::cout << "Array " << i << ": " << cellData->GetArrayName(i) << std::endl;
        }

	    std::cout << "Grid read from file: " << grid_filename << std::endl;

        locators.resize(numThreads);
        for (int i = 0; i < numThreads; ++i) {
            locators[i] = vtkSmartPointer<vtkStaticCellLocator>::New();
            locators[i]->SetDataSet(grid);
            locators[i]->BuildLocator();
        }

	    std::cout << "locators build successfully" << std::endl;

        // Call grid->FindCell to make it thread safe (it is not thread safe if called for the first time in a parallel region)
        double point[3] = {0.1, 0.1, 0.1};
        vtkSmartPointer<vtkGenericCell> GenCell = vtkSmartPointer<vtkGenericCell>::New();
        vtkIdType cellId = 0;
        double tol2 = 1e-6;
        int subId;
        double pcoords[3];
        double weights[8];
        grid->FindCell(point, nullptr, GenCell, cellId, tol2, subId, pcoords, weights);

        // Call cell->GetCell() to make it thread safe (it is not thread safe if called for the first time in a parallel region)
        grid->GetCell(cellId, GenCell);

        // Call grid->GetPointCells() to make it thread safe (it is not thread safe if called for the first time in a parallel region)
        vtkIdType pointId = 0;
        vtkSmartPointer<vtkIdList> cellList = vtkSmartPointer<vtkIdList>::New();
        grid->GetPointCells(pointId, cellList);

        // Get B-field array
        fieldArray = grid->GetPointData()->GetArray(B_field_name_m);
    }

    ~UnstructuredFieldContainer(){}

private:
    const char* B_field_name_m;
    vtkSmartPointer<vtkUnstructuredGrid> grid;
    std::vector<vtkSmartPointer<vtkUnstructuredGrid>> grids;
    std::vector<vtkSmartPointer<vtkStaticCellLocator>> locators;
    vtkSmartPointer<vtkDataArray> fieldArray;

    void readGrid(const char* grid_filename) {
        // Create a reader .vtk file
        vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
        reader->SetFileName(grid_filename);
        reader->Update();
        
        // Get the unstructured grid and copies for each thread from the reader
        int numThreads = Kokkos::DefaultExecutionSpace::concurrency();
        
        grid = reader->GetOutput();
        for (int i = 0; i < numThreads; ++i) {
            grids[i] = reader->GetOutput();
        }

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
    vtkIdType GetCellId(const Vector_t<double, Dim> R, double (&weights)[8]) {
        // Get reference position TODO: change this back for test grid
        Vector_t<double, Dim> R_ref = getReferencePosition(R);
        // Vector_t<double, Dim> R_ref = R;

        // Find the cell that contains the point (x, y, z)
        double point[Dim] = { R_ref[0], R_ref[1], R_ref[2] };
        double tol2 = 1e-6;
        vtkSmartPointer<vtkGenericCell> GenCell = vtkSmartPointer<vtkGenericCell>::New();
        double pcoords[Dim];

        int threadId = Kokkos::DefaultExecutionSpace::impl_hardware_thread_id();
        vtkIdType cellId = locators[threadId]->FindCell(point, tol2, GenCell, pcoords, weights);

        return cellId;
    }

    vtkIdType GetCellId(Vector_t<double, Dim> R) {        
        // Find the cell that contains the point (x, y, z)
        double weights[8];
        return GetCellId(R, weights);
    }

    vtkIdType CheckNeighboringCells(const Vector_t<double, Dim> R, double (&weights)[8], vtkIdType wrongId) {
        // Get reference position TODO: change this back for test grid
        Vector_t<double, Dim> R_ref = getReferencePosition(R);
        // Vector_t<double, Dim> R_ref = R;

        // Get the cell that does not contain the point, but was found by the locator
        vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();

        int threadId = Kokkos::DefaultExecutionSpace::impl_hardware_thread_id();
        grids[threadId]->GetCell(wrongId, cell);

        // Get the points of the cell
        vtkSmartPointer<vtkIdList> pointIds = cell->GetPointIds();

        // Over all points of the cell find cells that contain the point and check if point is inside the cell with correct interpolation weights
        for (int i = 0; i < pointIds->GetNumberOfIds(); i++) {
            // Get cells containing the point
            vtkSmartPointer<vtkIdList> cellList = vtkSmartPointer<vtkIdList>::New();
            grids[threadId]->GetPointCells(pointIds->GetId(i), cellList);

            for(int j = 0; j < cellList->GetNumberOfIds(); j++) {
                // Get the cell that contains the point
                vtkSmartPointer<vtkGenericCell> GenCell = vtkSmartPointer<vtkGenericCell>::New();
                vtkIdType cellId = cellList->GetId(j);
                grids[threadId]->GetCell(cellId, GenCell);

                // Check if the point is inside the cell
                double point[3] = { R_ref[0], R_ref[1], R_ref[2] };
                double closestPoint[3];
                int subId;
                double pcoords[Dim];
                double dist2;
                double w[8] = {0, 0, 0, 0, 0, 0, 0, 0};
                int insideTest = GenCell->EvaluatePosition(point, closestPoint, subId, pcoords, dist2, w);

                if(insideTest == 1 && cellId != wrongId) {
                    double weightSum = 0;
                    for (int k = 0; k < cell->GetNumberOfPoints(); k++) {
                        weightSum += w[k];
                    }
                    if(std::fabs(weightSum - 1) < 1e-6) {
                        for (int k = 0; k < 8; k++) {
                            weights[k] = w[k];
                        }
                        return cellId;
                    }
                }
            }
        }

        return -1;
    }

    vtkIdType FindCellAndInterpolateField(Vector_t<double, Dim> R, Vector_t<double, Dim> &interpolatedField, Vector_t<double, 8> &w) {
        // Get thread Id
        int threadId = Kokkos::DefaultExecutionSpace::impl_hardware_thread_id();

        // Set the interpolated field to zero
        Vector_t<double, Dim> B_field_ref = Vector_t<double, Dim>(0.0);
        // TODO: change this back after not plotting the weights
        w = Vector_t<double, 8>(0.0);

        // Get gridcell and interpolation weights
        double weights[8] = {0.0};
        vtkIdType cellId = GetCellId(R, weights);
        
        // Check that the point is inside the grid
        if(cellId == -1) {
            return cellId;
        }

        // get the cell that contains the point
        vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
        grids[threadId]->GetCell(cellId, cell);

        // Check if the interpolation weights are normalized
        double weightSum = 0.0;
        bool positiveWeights = true;
        assert(cell->GetNumberOfPoints() <= 8);

        for(unsigned i = 0; i < cell->GetNumberOfPoints(); i++) {
            weightSum += weights[i];
            if(weights[i] < -1e-10) {
                positiveWeights = false;
            }
        }

        if(std::fabs(weightSum - 1) > 1e-6 || !positiveWeights) {
            for(unsigned i = 0; i < cell->GetNumberOfPoints(); i++) {
                // TODO: change this back after not plotting the weights
                w[i] = weights[i];
            }
            return -3;
        }
        
        double B_val[Dim];     
        for(unsigned i = 0; i < cell->GetNumberOfPoints(); i++) {
            // TODO: change this back after not plotting the weights
            w[i] = weights[i];
            
            vtkIdType pointId = cell->GetPointId(i);
            this->fieldArray->GetTuple(pointId, B_val);
            for(unsigned j = 0; j < Dim; j++) {
                B_field_ref[j] += weights[i] * B_val[j];
            }
        }

        // Get the real field
        interpolatedField = getRealB_field(B_field_ref, R);

        return cellId;
    }

    void interpolateField(Vector_t<double, Dim> R, Vector_t<double, Dim> &interpolatedField)
    {
        FindCellAndInterpolateField(R, interpolatedField);
    }

    void getBoundaryInformation(vtkIdType cellID, std::string &boundaryType) {
        vtkSmartPointer<vtkStringArray> boundaryArray = vtkStringArray::SafeDownCast(grid->GetCellData()->GetAbstractArray("Boundary_Type"));
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
    
        boundaryType = boundaryArray->GetValue(cellID);
    }
};

#endif