#ifndef IPPL_UNSTRUCTURED_PARTICLE_CONTAINER_H
#define IPPL_UNSTRUCTURED_PARTICLE_CONTAINER_H

#include <memory>
#include "Manager/BaseManager.h"
#include "UnstructuredParticleBase.h"

template <typename T, unsigned Dim = 3>
class UnstructuredParticleContainer : public ippl::UnstructuredParticleBase<ippl::ParticleSpatialLayout<T, Dim>> {
    using scalar = _scalar;

    // Constructor for the Bunch class, taking a PLayout reference
    UnstructuredParticleContainer(ippl::ParticleSpatialLayout<T, Dim>& playout)
        : ippl::ParticleBase<ippl::ParticleSpatialLayout<T, Dim>>(playout) {
        
        // Add attributes to the particle bunch
        this->addAttribute(Ref_R);  // Reference position attribute in mesh
        this->addAttribute(Q);      // Charge attribute in C
        this->addAttribute(mass);   // Mass attribute in kg
        this->addAttribute(Ek);     // Energy attribute in eV
        this->addAttribute(V);      // Velocity attribute in km/s
        this->addAttribute(B);      // Magnetic field attribute in T
    }

    // Destructor for the Bunch class
    ~UnstructuredParticleContainer() {}

    // Declare instances of the attribute containers
    ippl::ParticleAttrib<ippl::Vector<scalar, Dim>> Ref_R;  // Reference position container
    ippl::ParticleAttrib<scalar> Q;                         // Charge container
    ippl::ParticleAttrib<scalar> mass;                      // Mass container
    ippl::ParticleAttrib<scalar> Ek;                        // Energy container
    ippl::ParticleAttrib<ippl::Vector<scalar, Dim>> V;      // Velocity container
    ippl::ParticleAttrib<ippl::Vector<scalar, Dim>> B;      // Magnetic field container
};

#endif