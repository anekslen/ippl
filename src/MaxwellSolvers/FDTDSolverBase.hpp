//
// Class FDTDSolverBase
//  Base class for solvers for Maxwell's equations using the FDTD method
//

namespace ippl {

    template <typename EMField, typename SourceField, fdtd_bc boundary_conditions>
    FDTDSolverBase<EMField, SourceField, boundary_conditions>::FDTDSolverBase(SourceField& source, EMField& E, EMField& B) {
        Maxwell<EMField, SourceField>::setSources(source);
        Maxwell<EMField, SourceField>::setEMFields(E, B);
    }

    template <typename EMField, typename SourceField, fdtd_bc boundary_conditions>
    void FDTDSolverBase<EMField, SourceField, boundary_conditions>::solve() {
        step();
        timeShift();
        evaluate_EB();
    }

    template <typename EMField, typename SourceField, fdtd_bc boundary_conditions>
    void FDTDSolverBase<EMField, SourceField, boundary_conditions>::setPeriodicBoundaryConditions() {
        periodic_bc = true;
        typename SourceField::BConds_t vector_bcs;
        auto bcsetter_single = [&vector_bcs]<size_t Idx>(const std::index_sequence<Idx>&) {
            vector_bcs[Idx] = std::make_shared<ippl::PeriodicFace<SourceField>>(Idx);
            return 0;
        };
        auto bcsetter = [bcsetter_single]<size_t... Idx>(const std::index_sequence<Idx...>&) {
            int x = (bcsetter_single(std::index_sequence<Idx>{}) ^ ...);
            (void)x;
        };
        bcsetter(std::make_index_sequence<Dim * 2>{});
        A_n.setFieldBC(vector_bcs);
        A_np1.setFieldBC(vector_bcs);
        A_nm1.setFieldBC(vector_bcs);
    }

    template <typename EMField, typename SourceField, fdtd_bc boundary_conditions>
    void FDTDSolverBase<EMField, SourceField, boundary_conditions>::timeShift() {
        // Look into this, maybe cyclic swap is better
        Kokkos::deep_copy(this->A_nm1.getView(), this->A_n.getView());
        Kokkos::deep_copy(this->A_n.getView(), this->A_np1.getView());
    }

    template <typename EMField, typename SourceField, fdtd_bc boundary_conditions>
    void FDTDSolverBase<EMField, SourceField, boundary_conditions>::applyBCs() {
        if constexpr (boundary_conditions == periodic) {
            A_n.getFieldBC().apply(A_n);
            A_nm1.getFieldBC().apply(A_nm1);
            A_np1.getFieldBC().apply(A_np1);
        } else {
            Vector<uint32_t, Dim> true_nr = nr_m;
            true_nr += (A_n.getNghost() * 2);
            second_order_mur_boundary_conditions bcs{};
            bcs.apply(A_n, A_nm1, A_np1, this->dt, true_nr, layout_mp->getLocalNDIndex());
        }
    }

    template <typename EMField, typename SourceField, fdtd_bc boundary_conditions>
    void FDTDSolverBase<EMField, SourceField, boundary_conditions>::evaluate_EB() {
        *(Maxwell<EMField, SourceField>::En_mp)     = typename EMField::value_type(0);
        *(Maxwell<EMField, SourceField>::Bn_mp)     = typename EMField::value_type(0);
        ippl::Vector<scalar, 3> inverse_2_spacing = ippl::Vector<scalar, 3>(0.5) / hr_m;
        const scalar idt                          = scalar(1.0) / dt;
        auto A_np1 = this->A_np1.getView(), A_n = this->A_n.getView(),
                A_nm1  = this->A_nm1.getView();
        auto source = Maxwell<EMField, SourceField>::JN_mp->getView();
        auto Eview  = Maxwell<EMField, SourceField>::En_mp->getView();
        auto Bview  = Maxwell<EMField, SourceField>::Bn_mp->getView();

        Kokkos::parallel_for(
            this->A_n.getFieldRangePolicy(), KOKKOS_LAMBDA(size_t i, size_t j, size_t k) {
                ippl::Vector<scalar, 3> dAdt =
                    (A_n(i, j, k).template tail<3>() - A_nm1(i, j, k).template tail<3>()) * idt;
                ippl::Vector<scalar, 4> dAdx =
                    (A_n(i + 1, j, k) - A_n(i - 1, j, k)) * inverse_2_spacing[0];
                ippl::Vector<scalar, 4> dAdy =
                    (A_n(i, j + 1, k) - A_n(i, j - 1, k)) * inverse_2_spacing[1];
                ippl::Vector<scalar, 4> dAdz =
                    (A_n(i, j, k + 1) - A_n(i, j, k - 1)) * inverse_2_spacing[2];

                ippl::Vector<scalar, 3> grad_phi{dAdx[0], dAdy[0], dAdz[0]};
                ippl::Vector<scalar, 3> curlA{
                    dAdy[3] - dAdz[2],
                    dAdz[1] - dAdx[3],
                    dAdx[2] - dAdy[1],
                };
                Eview(i, j, k) = -dAdt - grad_phi;
                Bview(i, j, k) = curlA;
            });
        Kokkos::fence();
    }
    
} // namespace ippl