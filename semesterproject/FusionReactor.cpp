// Fusion Reactor simulation
//   Usage:
//     ./FusionReactor <Np> <Nt> <t_method> <grid_file_name> --info 10

//     Np       = Total no. of particles in the simulation
//     Nt       = Number of time steps
//     t_method = Time-stepping method used e.g. Leapfrog
//     grid_file_name = path to the grid file in VTK format

//     Example:
//     ./FusionReactor 1 20 LeapFrog reactormesh/mesh.vtk reactormesh/InitialParticles.csv --info 10

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

        // Create an instance of a manger for the considered application
        FusionReactorManager<T, Dim> manager(totalP, nt, step_method, dt);

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