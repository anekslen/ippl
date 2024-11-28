#ifndef IPPL_PARTICLE_CONTAINER_H
#define IPPL_PARTICLE_CONTAINER_H

#include <memory>
#include "Manager/BaseManager.h"

template <typename T, unsigned Dim = 3>
class UParticleContainer : public ippl::ParticleBase<ippl::ParticleSpatialLayout<T, Dim>> {
    using Base = ippl::ParticleBase<ippl::ParticleSpatialLayout<T, Dim>>;

public:
    // Declare instances of the attribute containers
    typename Base::particle_position_type Ref_R;            // Reference position container
    typename Base::particle_position_type RReal;            // Real position container
    typename Base::particle_position_type RCylRef;          // Cylindrical coordinate reference position container
    typename Base::particle_position_type RCylReal;         // Cylindrical coordinate real position container
    ippl::ParticleAttrib<double> Q;                         // Charge container
    ippl::ParticleAttrib<double> mass;                      // Mass container
    ippl::ParticleAttrib<double> Ek;                        // Energy container
    typename Base::particle_position_type V;                // Velocity container
    typename Base::particle_position_type B;                // Magnetic field container
private:
    PLayout_t<T, Dim> pl_m;
public:
    // Constructor for the Bunch class, taking a PLayout reference
    UParticleContainer(Mesh_t<Dim>& mesh, FieldLayout_t<Dim>& FL)
    : pl_m(FL, mesh) {
        this->initialize(pl_m);     // Initialize ParticleSpatialLayout
        registerAttributes();
        setupBCs();
    }

    ~UParticleContainer(){}

    std::shared_ptr<PLayout_t<T, Dim>> getPL() { return pl_m; }
    void setPL(std::shared_ptr<PLayout_t<T, Dim>>& pl) { pl_m = pl; }

	void registerAttributes() {
		// Add attributes to the particle bunch
        // R is the reference position in the domain of the unstructured grid
        this->addAttribute(RReal);      // Real position attribute
        this->addAttribute(RCylRef);    // Cylindrical coordinate reference position attribute in the domain of the unstructured grid
        this->addAttribute(RCylReal);   // Cylindrical coordinate real position attribute
        this->addAttribute(Q);          // Charge attribute in C
        this->addAttribute(mass);       // Mass attribute in kg
        this->addAttribute(Ek);         // Energy attribute in J
        this->addAttribute(V);          // Velocity attribute in m/s
        this->addAttribute(B);          // Magnetic field attribute in T
	}
	void setupBCs() { setBCAllNo(); }

    private:
       void setBCAllNo() { this->setParticleBC(ippl::BC::NO); }

};
#endif
