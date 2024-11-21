//
// Class FDTDSolverBase
//  Base class for solvers for Maxwell's equations using the FDTD method
//

namespace ippl {
    
    template <typename EMField, typename SourceField, fdtd_bc boundary_conditions>
    NonStandardFDTDSolver<EMField, SourceField, boundary_conditions>::NonStandardFDTDSolver(SourceField& source, EMField& E, EMField& B) : FDTDSolverBase<EMField, SourceField>(source, E, B) {
        auto hx = source.get_mesh().getMeshSpacing();
        if ((hx[2] / hx[0]) * (hx[2] / hx[0]) + (hx[2] / hx[1]) * (hx[2] / hx[1]) >= 1) {
            std::cerr << "Dispersion-free CFL condition not satisfiable\n";
            std::abort();
        }
        initialize();
    }

    template <typename EMField, typename SourceField, fdtd_bc boundary_conditions>
    void NonStandardFDTDSolver<EMField, SourceField, boundary_conditions>::step() {
        const auto& ldom       = this->layout_mp->getLocalNDIndex();
        const int nghost       = this->A_n.getNghost();
        const auto aview       = this->A_n.getView();
        const auto anp1view    = this->A_np1.getView();
        const auto anm1view    = this->A_nm1.getView();
        const auto source_view = Maxwell<EMField, SourceField>::JN_mp->getView();

        const scalar calA = 0.25 * (1 + 0.02 / (Kokkos::pow(this->hr_m[2] / this->hr_m[0], 2) + Kokkos::pow(this->hr_m[2] / this->hr_m[1], 2)));
        nondispersive<scalar> ndisp{
            .a1 = 2
                    * (1 - (1 - 2 * calA) * Kokkos::pow(this->dt / this->hr_m[0], 2) - (1 - 2 * calA) * Kokkos::pow(this->dt / this->hr_m[1], 2)
                        - Kokkos::pow(this->dt / this->hr_m[2], 2)),
            .a2 = Kokkos::pow(this->dt / this->hr_m[0], 2),
            .a4 = Kokkos::pow(this->dt / this->hr_m[1], 2),
            .a6 = Kokkos::pow(this->dt / this->hr_m[2], 2) - 2 * calA * Kokkos::pow(this->dt / this->hr_m[0], 2) - 2 * calA * Kokkos::pow(this->dt / this->hr_m[1], 2),
            .a8 = Kokkos::pow(this->dt, 2)};
        Vector<uint32_t, Dim> true_nr = this->nr_m;
        true_nr += (nghost * 2);
        constexpr uint32_t one_if_absorbing_otherwise_0 =
            boundary_conditions == absorbing ? 1 : 0;
        Kokkos::parallel_for(
            ippl::getRangePolicy(aview, nghost), KOKKOS_LAMBDA(size_t i, size_t j, size_t k) {
                uint32_t ig = i + ldom.first()[0];
                uint32_t jg = j + ldom.first()[1];
                uint32_t kg = k + ldom.first()[2];
                uint32_t val =
                    uint32_t(ig == one_if_absorbing_otherwise_0)
                    + (uint32_t(jg == one_if_absorbing_otherwise_0) << 1)
                    + (uint32_t(kg == one_if_absorbing_otherwise_0) << 2)
                    + (uint32_t(ig == true_nr[0] - one_if_absorbing_otherwise_0 - 1) << 3)
                    + (uint32_t(jg == true_nr[1] - one_if_absorbing_otherwise_0 - 1) << 4)
                    + (uint32_t(kg == true_nr[2] - one_if_absorbing_otherwise_0 - 1) << 5);

                if (!val) {
                    anp1view(i, j, k) = -anm1view(i, j, k) + ndisp.a1 * aview(i, j, k)
                                        + ndisp.a2
                                                * (calA * aview(i + 1, j, k - 1)
                                                    + (1 - 2 * calA) * aview(i + 1, j, k)
                                                    + calA * aview(i + 1, j, k + 1))
                                        + ndisp.a2
                                                * (calA * aview(i - 1, j, k - 1)
                                                    + (1 - 2 * calA) * aview(i - 1, j, k)
                                                    + calA * aview(i - 1, j, k + 1))
                                        + ndisp.a4
                                                * (calA * aview(i, j + 1, k - 1)
                                                    + (1 - 2 * calA) * aview(i, j + 1, k)
                                                    + calA * aview(i, j + 1, k + 1))
                                        + ndisp.a4
                                                * (calA * aview(i, j - 1, k - 1)
                                                    + (1 - 2 * calA) * aview(i, j - 1, k)
                                                    + calA * aview(i, j - 1, k + 1))
                                        + ndisp.a6 * aview(i, j, k + 1)
                                        + ndisp.a6 * aview(i, j, k - 1)
                                        + ndisp.a8 * source_view(i, j, k);
                }
            });
        Kokkos::fence();
        this->A_np1.fillHalo();
        this->applyBCs();
    }

    template <typename EMField, typename SourceField, fdtd_bc boundary_conditions>
    void NonStandardFDTDSolver<EMField, SourceField, boundary_conditions>::initialize() {
        // get layout and mesh
        this->layout_mp = &(this->JN_mp->getLayout());
        this->mesh_mp   = &(this->JN_mp->get_mesh());

        // get mesh spacing, domain, and mesh size
        this->hr_m     = this->mesh_mp->getMeshSpacing();
        this->dt       = this->hr_m[2];
        this->domain_m = this->layout_mp->getDomain();
        for (unsigned int i = 0; i < Dim; ++i)
            this->nr_m[i] = this->domain_m[i].length();

        // initialize fields
        this->A_nm1.initialize(*this->mesh_mp, *this->layout_mp);
        this->A_n.initialize(*this->mesh_mp, *this->layout_mp);
        this->A_np1.initialize(*this->mesh_mp, *this->layout_mp);

        this->A_nm1 = 0.0;
        this->A_n   = 0.0;
        this->A_np1 = 0.0;
    };
}