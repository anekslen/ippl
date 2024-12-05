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

#include <vector>

#include "UnstructuredFieldContainer.hpp"
#include "UnstructuredGridManager.h"
#include "UParticleContainer.h"

#include "datatypes.h"

template <typename T>
KOKKOS_INLINE_FUNCTION Vector<T, 3> cross(const ippl::Vector<T, 3>& a, const ippl::Vector<T, 3>& b) {
    ippl::Vector<T, 3> ret;
    ret[0] = a[1] * b[2] - a[2] * b[1];
    ret[1] = a[2] * b[0] - a[0] * b[2];
    ret[2] = a[0] * b[1] - a[1] * b[0];
    return ret;
}

template <typename T, unsigned Dim>
class FusionReactorManager
    : public UnstructuredGridManager<T, Dim, UParticleContainer<T, Dim>, UnstructuredFieldContainer<T, Dim>> {
public:
    using UParticleContainer_t = UParticleContainer<T, Dim>;
    using UnstructuredFieldContainer_t = UnstructuredFieldContainer<T, Dim>;

    FusionReactorManager(size_type totalP_, int nt_, std::string& stepMethod_)
         : UnstructuredGridManager<T, Dim, UParticleContainer<T, Dim>, UnstructuredFieldContainer<T, Dim>>(totalP_, nt_, stepMethod_) {}

    ~FusionReactorManager(){}

    void pre_run(const char* grid_filename) {
        Inform m("Pre Run");

        this->setUFieldContainer(std::make_shared<UnstructuredFieldContainer_t>(grid_filename));

        // Initialize variables for dummy mesh
        this->nr_m = 100;
        for (unsigned i = 0; i < Dim; i++) {
            this->domain_m[i] = ippl::Index(this->nr_m[i]);
        }
        this->ufcontainer_m->getGridBounds(this->rmin_m, this->rmax_m);       // set mesh bounds to the bounds of the unstructured grid TODO: check this
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

        initializeParticles();

        this->dt_m     = 1e-8;  // TODO: how to choose timestep?
        this->it_m     = 0;
        this->time_m   = 0.0;

        dump();
    }

    void initializeParticles(){
        Inform m("Initialize Particles");

        std::shared_ptr<UParticleContainer_t> pc = this->pcontainer_m;
        std::shared_ptr<UnstructuredFieldContainer_t> ufc    = this->ufcontainer_m;

        static IpplTimings::TimerRef particleCreation = IpplTimings::getTimer("particlesCreation");
        IpplTimings::startTimer(particleCreation);

        pc->create(this->totalP_m);

        assert(this->totalP_m == 1);    // TODO: change this to read initial particle positions from a file

        // Set initial particle positions
        pc->RCylReal(0) = Vector_t<T, Dim>{0.318, (pi / 6.0), 0.43};    // Initial real position in cylindrical coordinates (30 °)
        pc->RReal(0) = Vector_t<T, Dim>{                                            // Initial real position in cartesian coordinates
            pc->RCylReal(0)[0] * Kokkos::cos(pc->RCylReal(0)[1]),
            pc->RCylReal(0)[0] * Kokkos::sin(pc->RCylReal(0)[1]),
            pc->RCylReal(0)[2]
        };

        std::cout << "Initial real position: " << pc->RReal(0) << std::endl;

        if(static_cast<int>(pc->RCylReal(0)[1] / (pi / 4.0)) % 2 == 0){   // Initial reference position in cylindrical coordinates
            pc->RCylRef(0) = Vector_t<T, Dim>{pc->RCylReal(0)[0], Kokkos::fmod(pc->RCylReal(0)[1], pi / 2.0), Kokkos::abs(pc->RCylReal(0)[2])};
        }
        else{
            pc->RCylRef(0) = Vector_t<T, Dim>{pc->RCylReal(0)[0], pi / 4.0 - Kokkos::fmod(pc->RCylReal(0)[1], pi / 4.0), Kokkos::abs(pc->RCylReal(0)[2])};
        }

        pc->R(0) = pc->RReal(0) = Vector_t<T, Dim>{                 // Initial reference position in cartesian coordinates
            pc->RCylRef(0)[0] * Kokkos::cos(pc->RCylRef(0)[1]),
            pc->RCylRef(0)[0] * Kokkos::sin(pc->RCylRef(0)[1]),
            pc->RCylRef(0)[2]
        };

        // Find cellId of the particle and interpolate the magnetic field at the initial position
        pc->cellId(0) = ufc->FindCellAndInterpolateField(pc->R(0), pc->B(0), "B_field");
        assert(pc->cellId(0) != -1);

        // Set initial particle charge in C and mass in kg
        pc->Q(0) = 1.602e-19;  
        pc->mass(0) = 1.673e-27;

        // Set initial particle energy in J (100 eV)
        pc->Ek(0) = 100 * 1.602e-19;

        // Set initial particle velocity in m/s
        double absv = Kokkos::sqrt(2.0 * pc->Ek(0) / pc->mass(0));
        pc->V(0) = Vector_t<T, Dim>{absv/Kokkos::sqrt(2), absv/Kokkos::sqrt(2), 0.0};
        m << "initial velocitly: " << pc->V(0) << endl;

        // Check if the velocity is set correctly
        assert(Kokkos::abs(pc->V(0).dot(pc->V(0)) - absv * absv) < 1e-10);

        m << "particle initialized" << endl;
        // Calculate initial change of velocity
        ufc->interpolateField(pc->R(0), pc->B(0), "B_field");
        m << "field interpolated" << pc->B(0) << endl;
        pc->dV(0) = pc->Q(0) * cross(pc->V(0), pc->B(0));
        m << "change of velocity calculated" << pc->dV(0) << endl;


        IpplTimings::stopTimer(particleCreation);
        m << "particles created and initial conditions assigned " << endl;
    }

    void advance() override {
        if (this->stepMethod_m == "LeapFrog") {
            LeapFrogStep();
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

        IpplTimings::startTimer(PTimer);
        pc->V = pc->V + 0.5 * dt * pc->dV;
        IpplTimings::stopTimer(PTimer);

        // drift
        IpplTimings::startTimer(RTimer);
        pc->RReal = pc->RReal + dt * pc->V;
        IpplTimings::stopTimer(RTimer);

        // Since the particles have moved calculate new reference position, check if particle is still in the grid
        IpplTimings::startTimer(UpdateTimer);
        // Position in cylindrical coordinates
        for(unsigned i = 0; i < pc->getTotalNum(); ++i){
            pc->RCylReal(i) = Vector_t<T, Dim>{
                Kokkos::sqrt(Kokkos::pow(pc->RReal(i)[0], 2) + Kokkos::pow(pc->RReal(i)[1], 2)),
                Kokkos::atan2(pc->RReal(i)[1], pc->RReal(i)[0]),
                pc->RReal(i)[2]
            };
        }

        // Reference position in cylindrical coordinates
        if(static_cast<int>(pc->RCylReal(0)[1] / (pi / 4.0)) % 2 == 0){
            for (unsigned i = 0; i < pc->getTotalNum(); ++i) {
                pc->RCylRef(i) = Vector_t<T, Dim>{pc->RCylReal(i)[0], Kokkos::fmod(pc->RCylReal(i)[1], pi / 2.0), Kokkos::abs(pc->RCylReal(i)[2])};
            }
        }
        else{
            for(unsigned i = 0; i < pc->getTotalNum(); ++i){
                pc->RCylRef(i) = {pc->RCylReal(i)[0], pi / 4.0 - Kokkos::fmod(pc->RCylReal(i)[1], pi / 4.0), Kokkos::abs(pc->RCylReal(i)[2])};
            }
        }

        // Reference position in cartesian coordinates
        for(unsigned i = 0; i < pc->getTotalNum(); ++i){
            pc->R(i) = Vector_t<T, Dim>{
                pc->RCylRef(i)[0] * Kokkos::cos(pc->RCylRef(i)[1]),
                pc->RCylRef(i)[0] * Kokkos::sin(pc->RCylRef(i)[1]),
                pc->RCylRef(i)[2]
            };
        }

        // Check if the particle position is inside the domain and interpolate the magnetic field at the new position
        for(unsigned i = 0; i < pc->getTotalNum(); ++i){
            pc->cellId(i) = ufc->FindCellAndInterpolateField(pc->R(i), pc->B(i), "B_field");
        }

        IpplTimings::stopTimer(UpdateTimer);

        // Calculate new change of velocity
        for(unsigned i = 0; i < pc->getTotalNum(); ++i){
            // If the particle is outside the grid set velocity and change of velocity to zero
            if(pc->cellId(i) == -1){
                pc->V(i) = 0.0;
                pc->dV(i) = 0.0;
            }
            else{
                pc->dV(i) = pc->Q(i) * cross(pc->V(i), pc->B(i));
                std::cout << "change of velocity calculated" << pc->dV(i) << ", Velocity = " << pc->V(i) << ", B Field = " << pc->B(i) << std::endl;
            }
        }

        // kick
        IpplTimings::startTimer(PTimer);
        pc->V = pc->V + 0.5 * dt * pc->dV;
        IpplTimings::stopTimer(PTimer);
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
                csvout << "Time,Particle_id,Position_x,Position_y,Position_z,Cell_id,Velocity_x,Velocity_y,Velocity_z" << endl;
            }
            for(unsigned i = 0; i < this->pcontainer_m->getTotalNum(); ++i){
                csvout << this->time_m << "," << i << "," << this->pcontainer_m->R(i)[0] << "," << this->pcontainer_m->R(i)[1] << "," << this->pcontainer_m->R(i)[2] << "," << this->pcontainer_m->cellId(i) << "," << this->pcontainer_m->V(i)[0] << "," << this->pcontainer_m->V(i)[1] << "," << this->pcontainer_m->V(i)[2] << endl;
            }
        }
    }

};
#endif