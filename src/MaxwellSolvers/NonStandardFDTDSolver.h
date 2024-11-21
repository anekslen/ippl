#ifndef IPPL_NON_STANDARD_FDTD_SOLVER_H
#define IPPL_NON_STANDARD_FDTD_SOLVER_H

#include <cstddef>
using std::size_t;
#include "Types/Vector.h"

#include "FieldLayout/FieldLayout.h"
#include "MaxwellSolvers/AbsorbingBC.h"
#include "MaxwellSolvers/Maxwell.h"
#include "MaxwellSolvers/FDTDSolverBase.h"
#include "Meshes/UniformCartesian.h"
#include "Particle/ParticleBase.h"

namespace ippl {

    template <typename EMField, typename SourceField, fdtd_bc boundary_conditions = periodic>
    class NonStandardFDTDSolver : public FDTDSolverBase<EMField, SourceField> {
    public:
        NonStandardFDTDSolver(SourceField& source, EMField& E, EMField& B);

        // TODO: Can this be inherited?
        constexpr static unsigned Dim = EMField::dim;
        using scalar                  = typename EMField::value_type::value_type;
        using Vector_t                = Vector<typename EMField::value_type::value_type, Dim>;
        using SourceVector_t            = typename SourceField::value_type;

        template <typename scalar>
        struct nondispersive {
            scalar a1;
            scalar a2;
            scalar a4;
            scalar a6;
            scalar a8;
        };

        void step() override;
        void initialize() override;
    };
} // namespace ippl

#include "NonStandardFDTDSolver.hpp"
#endif