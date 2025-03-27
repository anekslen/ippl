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

    FusionReactorManager(size_type totalP_, int nt_, std::string& stepMethod_, const char* output_folder_, double dt_)
         : UnstructuredGridManager<T, Dim, UParticleContainer<T, Dim>, UnstructuredFieldContainer<T, Dim>>(totalP_, nt_, stepMethod_, output_folder_, dt_) {}

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

        ufc->calculateCurl("B_Field", "Vorticity");
        ufc->calculateMagnitude("Vorticity", "VorticityMagnitude");

        std::string field_filename = std::string(this->output_folder_m) + "/curl.csv";
        const char* field_filename_cstr = field_filename.c_str();
        ufc->writeField(field_filename_cstr, "Vorticity", true, "VorticityMagnitude");

        // Create header for the output files

        // Print header for lost particles
        std::stringstream fname;
        fname << this->output_folder_m;
        fname << "/LostParticles_";
        fname << ippl::Comm->size();
        fname << "_manager";
        fname << ".csv";
        Inform lostout(NULL, fname.str().c_str(), Inform::APPEND);
        lostout.precision(16);
        lostout.setf(std::ios::scientific, std::ios::floatfield);
        lostout << "Time,Particle_id,Cell_id,Boundary_type" << endl;

        // Print header for exited particles
        std::stringstream efname;
        efname << this->output_folder_m;
        efname << "/ExitedParticles_";
        efname << ippl::Comm->size();
        efname << "_manager";
        efname << ".csv";
        Inform exitOut(NULL, efname.str().c_str(), Inform::APPEND);
        exitOut.precision(16);
        exitOut.setf(std::ios::scientific, std::ios::floatfield);

        exitOut << "Time,Particle_id,Cell_id" << endl;

        // Print header for missed cells
        std::stringstream mfname;
        mfname << this->output_folder_m;
        mfname << "/MissedCells_";
        mfname << ippl::Comm->size();
        mfname << "_manager";
        mfname << ".csv";
        Inform missOut(NULL, mfname.str().c_str(), Inform::APPEND);
        missOut.precision(16);
        missOut.setf(std::ios::scientific, std::ios::floatfield);

        missOut << "Time,Particle_id,Cell_id" << endl;

        // Print header for missed weights
        std::stringstream wfname;
        wfname << this->output_folder_m;
        wfname << "/MissedWeights_";
        wfname << ippl::Comm->size();
        wfname << "_manager";
        wfname << ".csv";
        Inform weightsOut(NULL, wfname.str().c_str(), Inform::APPEND);
        weightsOut.precision(16);
        weightsOut.setf(std::ios::scientific, std::ios::floatfield);

        weightsOut << "Time,Particle_id,Cell_id,weightSum,W0,W1,W2,W3,W4,W5,W6,W7" << endl;
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
            pc->Id(numP) = numP;
            
            ++numP;
        }

        // Ceck if the number of particles read from the file is equal to the total number of particles
        std::cout << "Number of particles read from the file: " << numP << ", Total of particles in the similation should be: " << this->totalP_m << std::endl;
        assert(numP == this->totalP_m);

        // Find cellId of the particle and interpolate the magnetic field at the initial position
        for(unsigned i = 0; i < pc->getTotalNum(); ++i){
            pc->cellId(i) = ufc->FindCellAndInterpolateField(pc->R(i), pc->B(i), pc->weights(i));
            assert(pc->cellId(i) != -1);    // Check if the initial particle position is inside the domain
            // m << "position particle " << i << " = " << pc->R(i) << endl;
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

            // std::cout << absv << ", " << SplitSum << std::endl;
            pc->V(i) = Vector_t<T, Dim>{absv*std::sqrt(pc->V(i)[0]/SplitSum), absv*std::sqrt(pc->V(i)[1]/SplitSum), absv*std::sqrt(pc->V(i)[2]/SplitSum)};
            //pc->V(i) = Vector_t<T, Dim>{absv/Kokkos::sqrt(3), absv/Kokkos::sqrt(3), absv/Kokkos::sqrt(3)};
            //pc->V(i) = Vector_t<T, Dim>{absv/Kokkos::sqrt(2), absv/Kokkos::sqrt(2), 0.0};
            //pc->V(i) = Vector_t<T, Dim>{0.0, 0.0, absv};
            // std::cout << "Particle " << i << " velocity: " << pc->V(i) << std::endl;
        }

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
            pc->cellId(i) = ufc->FindCellAndInterpolateField(pc->R(i), pc->B(i), pc->weights(i));
        }
        IpplTimings::stopTimer(UpdateTimer);


        IpplTimings::startTimer(PTimer);
        // Calculate new change of velocity
        for(unsigned i = 0; i < pc->getTotalNum(); ++i){
            // If the particle is outside the grid set velocity and change of velocity to zero
            if(pc->cellId(i) == -1 || pc->cellId(i) == -2 || pc->cellId(i) == -3){
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

        int lostNum = 0;

        Kokkos::View<bool*> invalid("invalid_particles", pc->getTotalNum());
        Kokkos::deep_copy(invalid, false);

        Kokkos::parallel_reduce("ParticleUpdate", Kokkos::RangePolicy<int>(0, pc->getTotalNum()), KOKKOS_CLASS_LAMBDA (const int& i, int& lLostNum) {
                IpplTimings::startTimer(UpdateTimer);
                // Check if the particle position is inside the domain and interpolate the magnetic field at the new position
                pc->cellIdOld(i) = pc->cellId(i);
                pc->cellId(i) = ufc->FindCellAndInterpolateField(pc->R(i), pc->B(i), pc->weights(i));
                IpplTimings::stopTimer(UpdateTimer);
            
                IpplTimings::startTimer(PTimer);
                // If the particle is outside the grid set velocity and change of velocity to zero
                if(pc->cellId(i) == -1 || pc->cellId(i) == -2 || pc->cellId(i) == -3){
                    invalid(i) = true;  // Mark particle as invalid
                    pc->V(i) = 0.0;
                    pc->dV(i) = 0.0;
                    lLostNum++;
                }
                // Calculate new calculated velocity using Boris method
                else{
                    Vector_t<T, Dim> t = 0.5 * pc->Q(i) * pc->B(i) / pc->mass(i) * dt;
                    Vector_t<T, Dim> s = 2.0 * t / (1.0 + ippl::dot(t, t).apply());

                    Vector_t<T, Dim> v_prime = pc->V(i) + ippl::cross(pc->V(i), t);
                    pc->dV(i) = ippl::cross(v_prime, s);
                }
                IpplTimings::stopTimer(PTimer);
        }, lostNum);

        // Print exit message for particles that left the grid
        for(unsigned i = 0; i < pc->getTotalNum(); ++i){
            if(invalid(i)) {
                dumpLostParticles(i);
            }
        }

        // In case all particle get destroyed, dump data to have final state
        if(pc->getTotalNum() - lostNum == 0) {
            dump();
        }

        // Destroy particles that left the grid
        pc->destroy(invalid, lostNum);

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
            fname << this->output_folder_m;
            fname << "/Particles_";
            fname << ippl::Comm->size();
            fname << "_manager";
            fname << ".csv";
            Inform csvout(NULL, fname.str().c_str(), Inform::APPEND);
            csvout.precision(16);
            csvout.setf(std::ios::scientific, std::ios::floatfield);
            
            csvout << "Time,Particle_id,Position_x,Position_y,Position_z,Cell_id,Velocity_x,Velocity_y,Velocity_z,E_kin,B_x,B_y,B_z,Mag_B,RotB_x,RotB_y,RotB_z,Mag_RotB,W0,W1,W2,W3,W4,W5,W6,W7" << endl;

            for(unsigned i = 0; i < this->pcontainer_m->getTotalNum(); ++i){
                Vector_t<T, Dim> CrossB = Vector_t<T, Dim>{(this->pcontainer_m->B(i)[2] - this->pcontainer_m->B(i)[1])
                                            , (this->pcontainer_m->B(i)[0] - this->pcontainer_m->B(i)[2])
                                            , (this->pcontainer_m->B(i)[1] - this->pcontainer_m->B(i)[0])};

                csvout << this->time_m << "," << this->pcontainer_m->Id(i) << ","
                << this->pcontainer_m->R(i)[0] << "," << this->pcontainer_m->R(i)[1] << "," << this->pcontainer_m->R(i)[2] << ","
                << this->pcontainer_m->cellId(i) << ","
                << this->pcontainer_m->V(i)[0] << "," << this->pcontainer_m->V(i)[1] << "," << this->pcontainer_m->V(i)[2] << ","
                << this->pcontainer_m->Ek(i) << ","
                << this->pcontainer_m->B(i)[0] << "," << this->pcontainer_m->B(i)[1] << "," << this->pcontainer_m->B(i)[2] << ","
                << std::sqrt(ippl::dot(this->pcontainer_m->B(i), this->pcontainer_m->B(i)).apply()) << ","
                << CrossB[0] << "," << CrossB[1] << "," << CrossB[2] << ","
                << std::sqrt(ippl::dot(CrossB, CrossB).apply()) << ","
                << this->pcontainer_m->weights(i)[0] << "," << this->pcontainer_m->weights(i)[1] << "," << this->pcontainer_m->weights(i)[2] << "," << this->pcontainer_m->weights(i)[3] << ","
                << this->pcontainer_m->weights(i)[4] << "," << this->pcontainer_m->weights(i)[5] << "," << this->pcontainer_m->weights(i)[6] << "," << this->pcontainer_m->weights(i)[7] << endl;
            }

        }
    }

    void dumpLostParticles(int i) {
        Inform m("dumpLostParticles");
        if (ippl::Comm->rank() == 0) {
            std::stringstream fname;
            fname << this->output_folder_m;
            fname << "/LostParticles_";
            fname << ippl::Comm->size();
            fname << "_manager";
            fname << ".csv";
            Inform lostout(NULL, fname.str().c_str(), Inform::APPEND);
            lostout.precision(16);
            lostout.setf(std::ios::scientific, std::ios::floatfield);

            int id = this->pcontainer_m->Id(i);
            
            std::string boundaryType;
            this->ufcontainer_m->getBoundaryInformation(this->pcontainer_m->cellIdOld(i), boundaryType);

            lostout << this->time_m << "," << this->pcontainer_m->Id(i) << "," << this->pcontainer_m->cellIdOld(i) << "," << boundaryType << endl;

            if(this->pcontainer_m->cellId(i) == -1) {
                std::stringstream efname;
                efname << this->output_folder_m;
                efname << "/ExitedParticles_";
                efname << ippl::Comm->size();
                efname << "_manager";
                efname << ".csv";
                Inform exitOut(NULL, efname.str().c_str(), Inform::APPEND);
                exitOut.precision(16);
                exitOut.setf(std::ios::scientific, std::ios::floatfield);

                exitOut << this->time_m << "," << id << "," << this->pcontainer_m->cellIdOld(i) << endl;
            }

            if(this->pcontainer_m->cellId(i) == -2) {
                std::stringstream mfname;
                mfname << this->output_folder_m;
                mfname << "/MissedCells_";
                mfname << ippl::Comm->size();
                mfname << "_manager";
                mfname << ".csv";
                Inform missOut(NULL, mfname.str().c_str(), Inform::APPEND);
                missOut.precision(16);
                missOut.setf(std::ios::scientific, std::ios::floatfield);

                missOut << this->time_m << "," << id << "," << this->pcontainer_m->cellIdOld(i) << endl;
            }

            if(this->pcontainer_m->cellId(i) == -3) {
                std::stringstream wfname;
                wfname << this->output_folder_m;
                wfname << "/MissedWeights_";
                wfname << ippl::Comm->size();
                wfname << "_manager";
                wfname << ".csv";
                Inform weightsOut(NULL, wfname.str().c_str(), Inform::APPEND);
                weightsOut.precision(16);
                weightsOut.setf(std::ios::scientific, std::ios::floatfield);

                double wSum = 0.0;
                for(unsigned j = 0; j < 8; j++){
                    wSum += this->pcontainer_m->weights(i)[j];
                }

                weightsOut << this->time_m << "," << id << "," << this->pcontainer_m->cellIdOld(i) << "," 
                << wSum << "," << this->pcontainer_m->weights(i)[0] << "," << this->pcontainer_m->weights(i)[1] << ","
                << this->pcontainer_m->weights(i)[2] << "," << this->pcontainer_m->weights(i)[3] << "," << this->pcontainer_m->weights(i)[4] << ","
                << this->pcontainer_m->weights(i)[5] << "," << this->pcontainer_m->weights(i)[6] << "," << this->pcontainer_m->weights(i)[7] << endl;
            }
        }
    }
};
#endif