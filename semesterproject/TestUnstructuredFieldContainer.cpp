// This file is used to test the UnstructuredFieldContainer class and its methods
// The file reads a grid from a .vtk file, calculates the curl and the magnitude of the B vector field at each point.
// The curl and the magnitude of the B vector field are then written to a .vtk file and a .csv file respectively.
// The file also tests the methond that finds the correct gridcell in which a particle is at
// and the interpolation of the magnetic field to a particle at a given point.

// Call the executable with the following arguments:
// ./FusionReactor <grid file name> <point file name> --Info 10

// The grid file is a .vtk file that contains the grid and the magnetic field
// The point file is a .csv file that contains the points at which the magnetic field should be interpolated as well as the the element number of the cell in which the point is located

#include "Ippl.h"
#include "Utility/IpplTimings.h"
#include "Utility/Inform.h"

#include <iostream>
#include <string>
#include <fstream>

#include "UnstructuredFieldContainer.hpp"
#include "datatypes.h"

constexpr unsigned Dim = 3;
using T                = double;
const char* TestName   = "FusionReactor";

using UnstructuredFieldContainer_t = UnstructuredFieldContainer<T, Dim>;

int main(int argc, char *argv[]) {
    ippl::initialize(argc, argv);
    {
        Inform msg(TestName);
        Inform msg2all(TestName, INFORM_ALL_NODES);

        // Read input parameter
        const char * grid_file_name = argv[1];
        const char * point_file_name = argv[2];

        // Create an instance of an unstructured field container
        std::shared_ptr<UnstructuredFieldContainer_t> ufcontainer = std::make_shared<UnstructuredFieldContainer_t>(grid_file_name);
        
        std::cout << "Grid read" << std::endl;

        // Calculate the curl and the magnitude of the B vector field at each point
        ufcontainer->calculateMagnitude("B_Field", "Mag_B_Field");
        ufcontainer->calculateCurl("B_Field", "Curl_B");
        
        ufcontainer->calculateMagnitude("Curl_B", "Mag_Curl_B_Field");
        
        std::cout << "Curl and Magnitude calculated" << std::endl;

        // Test writeGrid
        ufcontainer->writeGrid("outputGrid.vtk");
        
        std::cout << "Grid written" << std::endl;

        // Test writeField
        ufcontainer->writeField("B_Field.csv", "B_Field", true, "Mag_B_Field");
        ufcontainer->writeField("Curl_B.csv", "Curl_B", true, "Mag_Curl_B_Field");

        std::cout << "Field written" << std::endl;

        // Test getGridBounds
        ippl::Vector<T, Dim> min;
        ippl::Vector<T, Dim> max;
        ufcontainer->getGridBounds(min, max);

        std::cout << "Grid bounds: Min = " << min[0] << ", " << min[1] << ", " << min[2] << ", Max = " << max[0] << ", " << max[1] << ", " << max[2] << std::endl;

        // Test GetCellId, interpolateField and FindCellAndInterpolateField functions with the points and information given in the file point_file
        int pointId;
        ippl::Vector<T, Dim> R;
        vtkIdType LocatedCell;
        ippl::Vector<T, Dim> GivenB_Field;
        
        vtkIdType cellId;
        ippl::Vector<T, Dim> InterpolatedB_Field;
        
        
        // Get points and cellId from file point_file
        std::ifstream infile(point_file_name);
        if (!infile) {
            std::cerr << "Unable to open file points.txt";
            exit(1);
        }

        std::cout << "Reading points from file" << std::endl;

        while (infile >> pointId >> R[0] >> R[1] >> R[2] >> LocatedCell >> GivenB_Field[0] >> GivenB_Field[1] >> GivenB_Field[2]) {
            // Test GetCellId
            cellId = ufcontainer->GetCellId(R);
            if (cellId != LocatedCell) {
                std::cerr << "Error GetCellId: The cellId is not correct! Should be " << LocatedCell << " but is " << cellId << std::endl;
                exit(1);
            }

            std::cout << "CellId: " << cellId << std::endl;

            // Test FindCellAndInterpolateField

            Vector_t<double, 8> w;
            cellId = ufcontainer->FindCellAndInterpolateField(R, InterpolatedB_Field, w);
            if (cellId != LocatedCell) {
                std::cerr << "Error FindCellAndInterpolateField: The cellId is not correct! Should be " << LocatedCell << " but is " << cellId << std::endl;
                exit(1);
            }
            std::cout << "InterpolatedB_Field: " << InterpolatedB_Field[0] << ", " << InterpolatedB_Field[1] << ", " << InterpolatedB_Field[2] << std::endl;

        }

        infile.close();
    }
    ippl::finalize();

    return EXIT_SUCCESS;
}