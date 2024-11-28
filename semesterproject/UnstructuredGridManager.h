#ifndef IPPL_UNSTRUCTURED_GRID_MANAGER_H
#define IPPL_UNSTRUCTURED_GRID_MANAGER_H

#include <memory>

#include "UnstructuredFieldContainer.hpp"
#include "Manager/BaseManager.h"
#include "UParticleContainer.hpp"

template <typename T, unsigned Dim, class pc = UParticleContainer<T, Dim>, class ufc = UnstructuredFieldContainer<T, Dim>>
class UnstructuredGridManager : public BaseManager {
public:
    using UParticleContainer_t = UParticleContainer<T, Dim>;
    using UnstructuredFieldContainer_t = UnstructuredFieldContainer<T, Dim>;
protected:
    size_type totalP_m;
    int nt_m;
    std::string stepMethod_m;
public:
    UnstructuredGridManager(size_type totalP_, int nt_, Vector_t<int, std::string& stepMethod_)
        : BaseManager()
        , ufcontainer_m(nullptr)
        , pcontainer_m(nullptr)
        , totalP_m(totalP_)
        , nt_m(nt_)
        , stepMethod_m(stepMethod_){}
    ~UnstructuredGridManager(){}

protected:
    double time_m;                              // simulation time
    double dt_m;                                // time step
    int it_m;                                   // iteration number
    std::shared_ptr<ufc> ufcontainer_m;         // unstructured field container
    std::shared_ptr<pc> pcontainer_m;           // particle container

    // variables for dummy mesh used in particle container
    Vector_t<int, Dim> nr_m;
    ippl::NDIndex<Dim> domain_m;
    Vector_t<double, Dim> rmin_m;
    Vector_t<double, Dim> rmax_m;
    Vector_t<double, Dim> hr_m;
    Vector_t<double, Dim> origin_m;
    std::array<bool, Dim> decomp_m;
    bool isAllPeriodic_m;   // TODO: two times
    Mesh_t mesh_m

    // variables for dummy field layout used in particle container
    std::array<bool, Dim> decomp_m; // TODO: rename
    bool isAllPeriodic_m;
    FieldLayout_t fl_m;


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
        Inform m("Pre-step");
        m << "Done" << endl;
    }

    void post_step() override {
        // Update time
        this->time_m += this->dt_m;
        this->it_m++;
        // write solution to output file
        this->dump();

        Inform m("Post-step:");
        m << "Finished time step: " << this->it_m << " time: " << this->time_m << endl;
    }

    void grid2par() override { gatherCIC(); }

    void gatherCIC() {
        /*TODO: implement gatherCIC, even needed?*/
    }

    void par2grid() { scatterCIC(); }

    void scatterCIC() {
        /*TODO: implement scatterCIC, even needed?*/
	}

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
