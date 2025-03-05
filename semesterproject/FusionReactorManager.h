#ifndef IPPL_FUSION_REACTOR_MANAGER_H
#define IPPL_FUSION_REACTOR_MANAGER_H

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkGradientFilter.h>
#include <vtkArrayCalculator.h>
#include <vtkCellLocator.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>

#include <vector>

#include "UnstructuredFieldContainer.hpp"
#include "UnstructuredGridManager.h"
#include "UParticleContainer.h"

#include "datatypes.h"

template <typename T, unsigned Dim>
class FusionReactorManager
    : public UnstructuredGridManager<T, Dim, UParticleContainer<T, Dim>, UnstructuredFieldContainer<T, Dim>> {
public:
    using UParticleContainer_t = UParticleContainer<T, Dim>;
    using UnstructuredFieldContainer_t = UnstructuredFieldContainer<T, Dim>;

    FusionReactorManager(size_type totalP_, int nt_, std::string& stepMethod_, double dt_)
         : UnstructuredGridManager<T, Dim, UParticleContainer<T, Dim>, UnstructuredFieldContainer<T, Dim>>(totalP_, nt_, stepMethod_, dt_) {}

    ~FusionReactorManager(){}

    void pre_run(const char* grid_filename, const char* particles_filename) {
        Inform m("Pre Run");

        this->setUFieldContainer(std::make_shared<UnstructuredFieldContainer_t>(grid_filename));

        // Initialize variables for dummy mesh
        this->nr_m = 100;
        for (unsigned i = 0; i < Dim; i++) {
            this->domain_m[i] = ippl::Index(this->nr_m[i]);
        }
        this->ufcontainer_m->getGridBounds(this->rmin_m, this->rmax_m);         // set mesh bounds to the bounds of the unstructured grid TODO: check this
        this->hr_m = this->rmax_m / this->nr_m;                                 // calculate mesh spacing
        this->origin_m = this->rmin_m;

        // Create mesh
        this->mesh_m = Mesh_t<Dim>(this->domain_m, this->hr_m, this->origin_m);

        // Initialize variables for dummy field layout
        this->decomp_m = {false, false, false};                                 // No parallel decomposition needed
        this->isAllPeriodic_m = false;

        // Create dummy field layout
        this->fl_m = FieldLayout_t<Dim>(MPI_COMM_WORLD, this->domain_m, this->decomp_m, this->isAllPeriodic_m);

        this->setParticleContainer( std::make_shared<UParticleContainer_t>( this->mesh_m, this->fl_m));

        initializeParticles(particles_filename);
        
        //this->dt_m     = 1. / 50. * 2. * pi * this->pcontainer_m->mass(0) / (this->pcontainer_m->Q(0) * Kokkos::sqrt(ippl::dot(this->pcontainer_m->B(0), this->pcontainer_m->B(0)).apply()));  // TODO: how to choose timestep?
        //this->dt_m    = 1e-10;
        m << "Time step: " << this->dt_m << endl;
        this->it_m     = 0;
        this->time_m   = 0.0;


        // Calculate initial v_1/2
        std::shared_ptr<UParticleContainer_t> pc = this->pcontainer_m;
        std::shared_ptr<UnstructuredFieldContainer_t> ufc    = this->ufcontainer_m;

        pc->dV = pc->Q / pc->mass * ippl::cross(pc->V, pc->B);
        pc->V = pc->V + 0.5 * this->dt_m * pc->dV;

        dump();

        ufc->calculateCurl("B_Field", "Vorticity");
        ufc->calculateMagnitude("Vorticity", "VorticityMagnitude");
        ufc->writeField("data/curl.csv", "Vorticity", true, "VorticityMagnitude");
    }

    void initializeParticles(const char* particles_filename) {
        Inform m("Initialize Particles");

        std::shared_ptr<UParticleContainer_t> pc = this->pcontainer_m;
        std::shared_ptr<UnstructuredFieldContainer_t> ufc    = this->ufcontainer_m;

        static IpplTimings::TimerRef particleCreation = IpplTimings::getTimer("particlesCreation");
        IpplTimings::startTimer(particleCreation);

        pc->create(this->totalP_m);

        // Read particle data from a file
        std::ifstream file(particles_filename);  // Open the csv file containing particle data
        // Check if the file is opened successfully
        if (!file.is_open()) {
            std::cerr << "Error opening particle data file!" << std::endl;
        }

        // Vector to store data read from the CSV file
        Vector_t<double, Dim> RCyl;
        Vector_t<double, Dim> VelSplit;

        // Read the file line by line
        unsigned numP = 0;
        std::string header;
        std::string line;
        std::getline(file, line);   // Skip the header line
        std::stringstream ss(line);
        std::getline(ss, header, ',');

        bool isCylindrical = (header == "R" ? true : false);

        while (std::getline(file, line)) {
            std::stringstream ss(line);
            std::string value;

            // Read initial real position in cylindrical coordinates
            for(unsigned int i = 0; i < Dim; i++){
                std::getline(ss, value, ',');
                RCyl[i] = std::stod(value);
            }
            if(isCylindrical){
                // Calculate initial real position in cartesian coordinates
                pc->R(numP) = Vector_t<T, Dim>{
                    RCyl[0] * Kokkos::cos(RCyl[1]),
                    RCyl[0] * Kokkos::sin(RCyl[1]),
                    RCyl[2]};
            }
            else{
                pc->R(numP) = RCyl;
            }

            // Read initial Energy in eV
            std::getline(ss, value, ',');
            pc->Ek(numP) = std::stod(value) * 1.602e-19;    // convert eV to J

            // Read velocity split components
            for(unsigned int i = 0; i < Dim; i++){
                std::getline(ss, value, ',');
                VelSplit[i] = std::stod(value);
            }
            pc->V(numP) = VelSplit;
            
            ++numP;
        }

        // Ceck if the number of particles read from the file is equal to the total number of particles
        std::cout << "Number of particles read from the file: " << numP << ", Total of particles in the similation should be: " << this->totalP_m << std::endl;
        assert(numP == this->totalP_m);

        // Find cellId of the particle and interpolate the magnetic field at the initial position
        for(unsigned i = 0; i < pc->getTotalNum(); ++i){
            pc->cellId(i) = ufc->FindCellAndInterpolateField(pc->R(i), pc->B(i));
            assert(pc->cellId(i) != -1);    // Check if the initial particle position is inside the domain
            m << "position particle " << i << " = " << pc->R(i) << endl;
        }

        std::cout << "Particle data read from the file" << std::endl;

        // Set initial particle charge in C and mass in kg
        pc->Q = 1.602e-19;
        pc->mass = 1.673e-27;

        // Set initial particle velocity in m/s(
        double absv;
        double SplitSum;
        for(unsigned i = 0; i < pc->getTotalNum(); ++i){
            absv = Kokkos::sqrt(2.0 * pc->Ek(i) / pc->mass(i));
            SplitSum = 0.0;
            for(unsigned int j = 0; j < Dim; j++){
                SplitSum += pc->V(i)[j];
            }

            std::cout << absv << ", " << SplitSum << std::endl;
            pc->V(i) = Vector_t<T, Dim>{absv*std::sqrt(pc->V(i)[0]/SplitSum), absv*std::sqrt(pc->V(i)[1]/SplitSum), absv*std::sqrt(pc->V(i)[2]/SplitSum)};
            //pc->V(i) = Vector_t<T, Dim>{absv/Kokkos::sqrt(3), absv/Kokkos::sqrt(3), absv/Kokkos::sqrt(3)};
            //pc->V(i) = Vector_t<T, Dim>{absv/Kokkos::sqrt(2), absv/Kokkos::sqrt(2), 0.0};
            //pc->V(i) = Vector_t<T, Dim>{0.0, 0.0, absv};
            std::cout << "Particle " << i << " velocity: " << pc->V(i) << std::endl;
        }

        // Check if the velocity is set correctly
        assert((Kokkos::abs(ippl::dot(pc->V, pc->V)) - absv * absv) < 1e-10);

        // Calculate kinetic Energy
        pc->Ek = 0.5 * pc->mass * ippl::dot(pc->V, pc->V);

        IpplTimings::stopTimer(particleCreation);
        m << "particles created and initial conditions assigned " << endl;
    }

    void advance() override {
        if (this->stepMethod_m == "LeapFrog") {
            LeapFrogStep();
        }
        else if (this->stepMethod_m == "Boris") {
            BorisStep();
        }
        else{
                throw IpplException(TestName, "Step method is not set/recognized!");
            }
    }

    void LeapFrogStep(){
        // LeapFrog time stepping https://en.wikipedia.org/wiki/Leapfrog_integration
        static IpplTimings::TimerRef PTimer           = IpplTimings::getTimer("pushVelocity");
        static IpplTimings::TimerRef RTimer           = IpplTimings::getTimer("pushPosition");
        static IpplTimings::TimerRef UpdateTimer      = IpplTimings::getTimer("updatingParticles");

        double dt                                           = this->dt_m;
        std::shared_ptr<UParticleContainer_t> pc            = this->pcontainer_m;
        std::shared_ptr<UnstructuredFieldContainer_t> ufc   = this->ufcontainer_m;

        // drift
        IpplTimings::startTimer(RTimer);
        pc->R = pc->R + dt * pc->V;
        IpplTimings::stopTimer(RTimer);

        // Check if the particle position is inside the domain and interpolate the magnetic field at the new position
        IpplTimings::startTimer(UpdateTimer);
        for(unsigned i = 0; i < pc->getTotalNum(); ++i){
            pc->cellId(i) = ufc->FindCellAndInterpolateField(pc->R(i), pc->B(i));
        }
        IpplTimings::stopTimer(UpdateTimer);


        IpplTimings::startTimer(PTimer);
        // Calculate new change of velocity
        for(unsigned i = 0; i < pc->getTotalNum(); ++i){
            // If the particle is outside the grid set velocity and change of velocity to zero
            if(pc->cellId(i) == -1){
                pc->V(i) = 0.0;
                pc->dV(i) = 0.0;
            }
            else{
                pc->dV(i) = pc->Q(i) / pc->mass(i) * ippl::cross(pc->V(i), pc->B(i));
            }
        }

        // kick
        pc->V = pc->V + dt * pc->dV;
        IpplTimings::stopTimer(PTimer);

        // Calculate kinetic Energy
        pc->Ek = 0.5 * pc->mass * ippl::dot(pc->V, pc->V);
    }

    void BorisStep(){
        // LeapFrog time stepping https://en.wikipedia.org/wiki/Leapfrog_integration
        static IpplTimings::TimerRef PTimer           = IpplTimings::getTimer("pushVelocity");
        static IpplTimings::TimerRef RTimer           = IpplTimings::getTimer("pushPosition");
        static IpplTimings::TimerRef UpdateTimer      = IpplTimings::getTimer("updatingParticles");

        double dt                                           = this->dt_m;
        std::shared_ptr<UParticleContainer_t> pc            = this->pcontainer_m;
        std::shared_ptr<UnstructuredFieldContainer_t> ufc   = this->ufcontainer_m;

        // drift
        IpplTimings::startTimer(RTimer);
        pc->R = pc->R + dt * pc->V;
        IpplTimings::stopTimer(RTimer);

        // TODO: change this back after fixing the randomness issues

        int id = 0;
        while (pc->cellId(id) == -1) {
            id++;
        }
        IpplTimings::startTimer(UpdateTimer);
        pc->cellIdOld(id) = pc->cellId(id);
        pc->cellId(id) = ufc->FindCellAndInterpolateField(pc->R(id), pc->B(id), "B_field");
        IpplTimings::stopTimer(UpdateTimer);
        
        Kokkos::parallel_for("ParticleUpdate", pc->getTotalNum(), KOKKOS_LAMBDA (const int i) {
            // Check if the particle position is inside the domain and interpolate the magnetic field at the new position
            if(i != id){
                pc->B(i) = pc->B(id);
                pc->cellIdOld(i) = pc->cellId(i);
                pc->cellId(i) = ufc->GetGridCell(pc->R(i));
            }
            if(pc->cellId(i) != -1){
                pc->cellId(i) = ufc->FindCellAndInterpolateField(pc->R(i), pc->B(i));
            }
            IpplTimings::stopTimer(UpdateTimer);

                IpplTimings::startTimer(PTimer);
                // If the particle is outside the grid set velocity and change of velocity to zero
                if(pc->cellId(i) == -1){
                    dumpLostParticles(i);
                    pc->V(i) = 0.0;
                    pc->dV(i) = 0.0;
                    // Adjust number of particles in simulation
                    this->setTotalP(this->totalP_m - 1);
                }
    
                else{
                    Vector_t<T, Dim> t = 0.5 * pc->Q(i) * pc->B(i) / pc->mass(i) * dt;
                    Vector_t<T, Dim> s = 2.0 * t / (1.0 + ippl::dot(t, t).apply());
    
                    Vector_t<T, Dim> v_prime = pc->V(i) + ippl::cross(pc->V(i), t);
                    pc->dV(i) = ippl::cross(v_prime, s);
                }
                IpplTimings::stopTimer(PTimer);
            }
        });

        IpplTimings::startTimer(PTimer);
        // kick (multiplication by dt not needed, already done in dV)
        pc->V = pc->V + pc->dV;
        IpplTimings::stopTimer(PTimer);

        // Calculate kinetic Energy
        pc->Ek = 0.5 * pc->mass * ippl::dot(pc->V, pc->V);
    }

    void dump() override {
        static IpplTimings::TimerRef dumpDataTimer = IpplTimings::getTimer("dumpData");
        IpplTimings::startTimer(dumpDataTimer);
        dumpPositions();
        IpplTimings::stopTimer(dumpDataTimer);
    }

    void dumpPositions() {
        Inform m("dumpPositions");

        if (ippl::Comm->rank() == 0) {
            std::stringstream fname;
            fname << "data/Particles_";
            fname << ippl::Comm->size();
            fname << "_manager";
            fname << ".csv";
            Inform csvout(NULL, fname.str().c_str(), Inform::APPEND);
            csvout.precision(16);
            csvout.setf(std::ios::scientific, std::ios::floatfield);
            if ( std::fabs(this->time_m) < 1e-14 ) {
                csvout << "Time,Particle_id,Position_x,Position_y,Position_z,Cell_id,Velocity_x,Velocity_y,Velocity_z,E_kin,B_x,B_y,B_z,Mag_B,RotB_x,RotB_y,RotB_z,Mag_RotB" << endl;
            }
            for(unsigned i = 0; i < this->pcontainer_m->getTotalNum(); ++i){
                Vector_t<T, Dim> CrossB = Vector_t<T, Dim>{(this->pcontainer_m->B(i)[2] - this->pcontainer_m->B(i)[1])
                                            , (this->pcontainer_m->B(i)[0] - this->pcontainer_m->B(i)[2])
                                            , (this->pcontainer_m->B(i)[1] - this->pcontainer_m->B(i)[0])};

                csvout << this->time_m << "," << i << ","
                << this->pcontainer_m->R(i)[0] << "," << this->pcontainer_m->R(i)[1] << "," << this->pcontainer_m->R(i)[2] << ","
                << this->pcontainer_m->cellId(i) << ","
                << this->pcontainer_m->V(i)[0] << "," << this->pcontainer_m->V(i)[1] << "," << this->pcontainer_m->V(i)[2] << ","
                << this->pcontainer_m->Ek(i) << ","
                << this->pcontainer_m->B(i)[0] << "," << this->pcontainer_m->B(i)[1] << "," << this->pcontainer_m->B(i)[2] << ","
                << std::sqrt(ippl::dot(this->pcontainer_m->B(i), this->pcontainer_m->B(i)).apply()) << ","
                << CrossB[0] << "," << CrossB[1] << "," << CrossB[2] << ","
                << std::sqrt(ippl::dot(CrossB, CrossB).apply()) << endl;
            }
        }
    }

    void dumpLostParticles(int i) {
        Inform m("dumpLostParticles");
        if (ippl::Comm->rank() == 0) {
            std::string boundaryType;
            this->ufcontainer_m->getBoundaryInformation(this->pcontainer_m->cellIdOld(i), boundaryType);

            std::stringstream fname;
            fname << "data/LostParticles_";
            fname << ippl::Comm->size();
            fname << "_manager";
            fname << ".csv";
            Inform lostout(NULL, fname.str().c_str(), Inform::APPEND);
            lostout.precision(16);
            lostout.setf(std::ios::scientific, std::ios::floatfield);

            lostout << "Particle " << i << " left the grid at time " << this->time_m << " at " << boundaryType << ", last cell number is " << this->pcontainer_m->cellIdOld(i) << "." << endl;

            std::cout << "Particle " << i << " left the grid at " << boundaryType << ", last cell number is " << this->pcontainer_m->cellIdOld(i) << "." << std::endl;
        }
    }
};
#endif