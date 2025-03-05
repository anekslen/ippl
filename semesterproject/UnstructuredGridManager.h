#ifndef IPPL_UNSTRUCTURED_GRID_MANAGER_H
#define IPPL_UNSTRUCTURED_GRID_MANAGER_H

#include <memory>

#include "Manager/BaseManager.h"
#include "UnstructuredFieldContainer.hpp"
#include "UParticleContainer.h"

template <typename T, unsigned Dim, class pc = UParticleContainer<T, Dim>, class ufc = UnstructuredFieldContainer<T, Dim>>
class UnstructuredGridManager : public ippl::BaseManager {
public:
    UnstructuredGridManager(size_type totalP_, int nt_, std::string& stepMethod_, double dt_)
        : ippl::BaseManager()
        , totalP_m(totalP_)
        , nt_m(nt_)
        , stepMethod_m(stepMethod_)
        , dt_m(dt_)
        , pcontainer_m(nullptr)
        , ufcontainer_m(nullptr) {}
    ~UnstructuredGridManager(){}

protected:
    size_type totalP_m;                         // total number of particles
    int nt_m;                                   // number of time steps
    std::string stepMethod_m;                   // time integration method (only LeapFrog is implemented)
    
    double time_m;                              // simulation time
    double dt_m;                                // time step
    int it_m;                                   // iteration number

    std::shared_ptr<pc> pcontainer_m;           // particle container
    std::shared_ptr<ufc> ufcontainer_m;         // unstructured field container


    // variables for dummy mesh used in particle container
    Vector_t<int, Dim> nr_m;
    ippl::NDIndex<Dim> domain_m;
    Vector_t<double, Dim> rmin_m;
    Vector_t<double, Dim> rmax_m;
    Vector_t<double, Dim> hr_m;
    Vector_t<double, Dim> origin_m;
    Mesh_t<Dim> mesh_m;

    // variables for dummy field layout used in particle container
    std::array<bool, Dim> decomp_m;
    bool isAllPeriodic_m;
    FieldLayout_t<Dim> fl_m;


public:
    size_type getTotalP() const { return totalP_m; }

    void setTotalP(size_type totalP_) { totalP_m = totalP_; }

    int getNt() const { return nt_m; }

    void setNt(int nt_) { nt_m = nt_; }

    const std::string& getStepMethod() const { return stepMethod_m; }

    void setStepMethod(const std::string& stepMethod_) { stepMethod_m = stepMethod_; }

    double getTime() const { return time_m; }

    void setTime(double time_) { time_m = time_; }

    virtual void dump() { /* default does nothing */ };

    void pre_step() override {
        /*
        Inform m("Pre-step");
        m << "Done" << endl;
        */
    }

    void post_step() override {
        // Update time
        this->time_m += this->dt_m;
        this->it_m++;
        // write solution to output file
        this->dump();

        if(this->it_m % 100 == 0){
            Inform m("Post-step:");
            m << "Finished time step: " << this->it_m << " time: " << this->time_m << endl;
        }
    }

    /*
    Not needed functions, TODO: possibly delete later
    void grid2par() override { gatherCIC(); }

    void gatherCIC() {
        //TODO: implement gatherCIC, even needed?
    }

    void par2grid() { scatterCIC(); }

    void scatterCIC() {
        //TODO: implement scatterCIC, even needed?
	}
    */

    std::shared_ptr<pc> getParticleContainer() {
        return pcontainer_m;
    }

    void setParticleContainer(std::shared_ptr<pc> pcontainer){
        pcontainer_m = pcontainer;
    }

    std::shared_ptr<ufc> getFieldContainer() {
        return ufcontainer_m;
    }

    void setUFieldContainer(std::shared_ptr<ufc> ufcontainer){
        ufcontainer_m = ufcontainer;
    }
};
#endif
