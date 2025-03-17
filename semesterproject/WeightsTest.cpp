#include "Ippl.h"
#include "Utility/IpplTimings.h"
#include "Utility/Inform.h"

#include <iostream>
#include <string>
#include <fstream>

#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkCellLocator.h>
#include <vtkAbstractArray.h>
#include <vtkStringArray.h>

constexpr unsigned Dim = 3;
using T                = double;
const char* TestName   = "WeightsTest";
template <typename T, unsigned Dim>
using Vector_t = ippl::Vector<T, Dim>;

int main(int argc, char *argv[]) {
    ippl::initialize(argc, argv);
    {
        Inform msg(TestName);
        Inform msg2all(TestName, INFORM_ALL_NODES);
        
        std::stringstream fname;
        fname << "data/WeightTest.csv";
        Inform csvout(NULL, fname.str().c_str(), Inform::APPEND);
        csvout.precision(16);
        csvout.setf(std::ios::scientific, std::ios::floatfield);

        // Read input parameter
        const char * grid_filename = argv[1];

        // Create the grid
        vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
        reader->SetFileName(grid_filename);
        reader->Update();
        
        vtkSmartPointer<vtkUnstructuredGrid> grid = reader->GetOutput();
        
        // Create the locator
        vtkSmartPointer<vtkCellLocator> locator = vtkSmartPointer<vtkCellLocator>::New();
        locator->SetDataSet(grid);
        locator->BuildLocator();
        
        // Get gridpoints
        vtkSmartPointer<vtkPointData> points = grid->GetPointData();
        
        // Get the magnetic field
        vtkSmartPointer<vtkDataArray> B_field = grid->GetPointData()->GetArray("B_Field");

        // Create array to store where the interpolation weights are wrong
        vtkSmartPointer<vtkIntArray> wrongInterpolation = vtkSmartPointer<vtkIntArray>::New();
        wrongInterpolation->SetNumberOfComponents(1);
        wrongInterpolation->SetNumberOfTuples(grid->GetNumberOfPoints());
        wrongInterpolation->SetName("wrongInterpolation");

        // Loop over all points in the grid and interpolate the magnetic field
        for(int i = 0; i < grid->GetNumberOfPoints(); i++){
            // Get the point from the grid
            double point[3];
            grid->GetPoint(i, point);
            
            // Get the magnetic field at the point
            double B[3];
            B_field->GetTuple(i, B);
            
            // Interpolate the filed at the point
            
            // Get cell and Interpolationweights
            double tol2 = 1e-6;
            vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
            double pcoords[3];
            double weights[8] = {0, 0, 0, 0, 0, 0, 0, 0};
            vtkIdType cellId = locator->FindCell(point, tol2, cell, pcoords, weights);
            
            if(cellId == -1) {
                std::cout << "Point " << i << " not found." << std::endl;
                wrongInterpolation->SetComponent(i, 0, 2);
            }
            else{
            
                // Check if the interpolation weights are normalized
                double weightSum = 0;
                for(int j = 0; j < cell->GetNumberOfPoints(); j++){
                    weightSum += weights[j];
                    if(weights[j] < -1e-10) {
                        std::cout << "Weight a lot smaller than 0" << std::endl;
                    }
                }
                
                if(std::abs(weightSum - 1) > 1e-12){
                    std::cout << "Interpolation weights are not normalized for point " << i << std::endl;
                    wrongInterpolation->SetComponent(i, 0, 3);
                }
                else{
                    
                    // Interpolate the magnetic field
                    double B_interpolated[3] = {0, 0, 0};
                    double B_temp[3];
                    for(int j = 0; j < cell->GetNumberOfPoints(); j++){
                        B_field->GetTuple(cell->GetPointId(j), B_temp);
                        for(int k = 0; k < 3; k++){
                            B_interpolated[k] += weights[j] * B_temp[k];
                        }
                    }
                    
                    // Check if the interpolated field is the same as the field at the point
                    for(int j = 0; j < 3; j++){
                        double diff = std::abs(B[j] - B_interpolated[j]);
                        if(diff > 1e-12){
                            //std::cout << "Interpolated field is not the same as the field at the point for point " << i << ", error is " << diff << "." << std::endl;
                            wrongInterpolation->SetComponent(i, 0, 1);
                            break;
                        }
                        wrongInterpolation->SetComponent(i, 0, 0);
                    }
                }
                csvout << weights[0] << ", " << weights[1] << ", " << weights[2] << ", " << weights[3] << ", " << weights[4] << ", " << weights[5] << ", " << weights[6] << ", " << weights[7] << endl;
            }
        }

        // Add the array to the grid
        grid->GetPointData()->AddArray(wrongInterpolation);

        // Write the grid to a file
        vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
        writer->SetFileName("data/InterpolationTest.vtk");
        writer->SetInputData(grid);
        writer->Write();
    }
    ippl::finalize();

    return EXIT_SUCCESS;
}