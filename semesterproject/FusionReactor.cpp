// Fusion Reactor simulation
//   Usage:
//     ./FusionReactor <Np> <Nt> <step_method> <grid_file_name> <input_particles_file_name> <dt> <output_folder> <writeData> --info 10

//     Np       = Total no. of particles in the simulation
//     Nt       = Number of time steps
//     step_method = Time-stepping method used e.g. Leapfrog
//     grid_file_name = path to the grid file in VTK format
//     input_particles_file_name = path to the file containing the initial particle positions
//     dt      = time step size
//     output_folder = path to the folder where the output files will be written
//     writeData = true/false if the data should be written for each time step not only the last position of each particle

//     Example:
//     ./FusionReactor 1 20 Boris reactormesh/mesh.vtk reactormesh/InitialParticles.csv 1e-8 ./data true --info 10

constexpr unsigned Dim = 3;
using T                = double;
const char* TestName   = "FusionReactor";

#include "Ippl.h"

#include <Kokkos_MathematicalConstants.hpp>
#include <Kokkos_MathematicalFunctions.hpp>
#include <Kokkos_Random.hpp>
#include <chrono>
#include <iostream>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "datatypes.h"

#include "Utility/IpplTimings.h"

#include "FusionReactorManager.h"

int main(int argc, char* argv[]) {
    ippl::initialize(argc, argv);
    {
        Inform msg(TestName);
        Inform msg2all(TestName, INFORM_ALL_NODES);

        static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("total");
        IpplTimings::startTimer(mainTimer);

        // Read input parameters, assign them to the corresponding memebers of manager
        int arg = 1;
        size_type totalP   = std::atoll(argv[arg++]);
        int nt             = std::atoi(argv[arg++]);
        std::string step_method = argv[arg++];
        const char* grid_filename = argv[arg++];
        const char* particles_filename = argv[arg++];
        double dt = std::atof(argv[arg++]);
        
        const char* output_folder;
        if (argc > arg) {
            output_folder = argv[arg++];
            std::cout << "output_folder: " << output_folder << std::endl;
        } else {
            output_folder = "./data";
        }

        bool writeData;
        if (argc > arg) {
            std::string boolStr = argv[arg++];
            std::transform(boolStr.begin(), boolStr.end(), boolStr.begin(), ::tolower);
            if (boolStr == "true" || boolStr == "1") {
                writeData = true;
            } else if (boolStr == "false" || boolStr == "0") {
                writeData = false;
            } else {
                std::cerr << "Invalid boolean value: " << boolStr << ". Use 'true'/'false' or '1'/'0'.\n";
                exit(EXIT_FAILURE);
            }
        } else {
            writeData = false;
        }

        // Create an instance of a manger for the considered application
        FusionReactorManager<T, Dim> manager(totalP, nt, step_method, output_folder, dt, writeData);

        // Perform pre-run operations, including creating mesh, particles,...
        manager.pre_run(grid_filename, particles_filename);

        manager.setTime(0.0);

        msg << "Starting iterations ..." << endl;

        manager.run(manager.getNt());

        msg << "End." << endl;

        IpplTimings::stopTimer(mainTimer);
        IpplTimings::print();
        IpplTimings::print(std::string("timing.dat"));
    }
    ippl::finalize();

    return 0;
}