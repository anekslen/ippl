#ifndef IPPL_FDTD_H
#define IPPL_FDTD_H
#include <cstddef>
using std::size_t;
#include "Types/Vector.h"

#include "FieldLayout/FieldLayout.h"
#include "MaxwellSolvers/AbsorbingBC.h"
#include "MaxwellSolvers/Maxwell.h"
#include "Meshes/UniformCartesian.h"
#include "Particle/ParticleBase.h"

namespace ippl {
    enum fdtd_bc {
        periodic,
        absorbing
    };

    template <typename EMField, typename SourceField, fdtd_bc boundary_conditions = periodic>
    class FDTDSolverBase : public Maxwell<EMField, SourceField> {
    public:
        constexpr static unsigned Dim = EMField::dim;
        using scalar                  = typename EMField::value_type::value_type;
        using Vector_t                = Vector<typename EMField::value_type::value_type, Dim>;
        using SourceVector_t            = typename SourceField::value_type;

        FDTDSolverBase(SourceField& source, EMField& E, EMField& B);
        void solve() override;

        void setPeriodicBoundaryConditions();

        SourceField A_n;
        SourceField A_np1;
        SourceField A_nm1;
        scalar dt;

    protected:
        void timeShift();
        void applyBCs();

    public:
        virtual void step() = 0;
        void evaluate_EB();

        bool periodic_bc;

        typename SourceField::Mesh_t* mesh_mp;
        FieldLayout<Dim>* layout_mp;
        NDIndex<Dim> domain_m;
        Vector_t hr_m;

        Vector<int, Dim> nr_m;

        virtual void initialize() = 0;
    };
}   // namespace ippl

#include "FDTDSolverBase.hpp"
#endif