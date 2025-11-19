#ifndef IPPL_FINITE_ELEMENT_SPACE_TRAITS_H
#define IPPL_FINITE_ELEMENT_SPACE_TRAITS_H

#include "FEM/Entity.h"
#include <tuple>
#include <type_traits>

namespace ippl {

    ///////////////////////////////////////////////////////////////////////
    // Space Type Tags ////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    
    struct LagrangeSpaceTag {};
    struct NedelecSpaceTag {};
    struct RaviartThomasSpaceTag {};
    
    ///////////////////////////////////////////////////////////////////////
    // Order-Aware Space Traits ///////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    
    template <typename SpaceTag, unsigned Dim, unsigned Order>
    struct FiniteElementSpaceTraits {
        static constexpr unsigned dim = Dim;
        static constexpr unsigned order = Order;

        // EntityTypes and DOFNums need to be specialized for each space type
        using EntityTypes = std::tuple<>;
        using DOFNums = std::tuple<>;

        // Tag helper to check for entity presence and get index
        ippl::TagIndex<EntityTypes> tagIndex;

        template <typename EntityType>
        static constexpr bool hasEntity() {
            return tagIndex.template contains<EntityType>();
        }

        template <typename EntityType>
        static constexpr unsigned entityDOFCount() {
            static_assert(hasEntity<EntityType>(), "EntityType not found in this FiniteElementSpaceTraits");
            constexpr unsigned index = tagIndex.template index<EntityType>();
            return std::tuple_element_t<index, DOFNums>::value;
        }
    };
    
    ///////////////////////////////////////////////////////////////////////
    // Lagrange Space Specializations ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    
    // 1D Lagrange - Order 1 (Linear elements)
    template <>
    struct FiniteElementSpaceTraits<LagrangeSpaceTag, 1, 1> {
        // Only vertices have DOFs in first-order elements
        using EntityTypes = std::tuple<Vertex<1>>;
        using DOFNums = std::tuple<std::integral_constant<unsigned, 1>>;

        static constexpr unsigned dofsPerElement = 2; // 2 vertices
    };

    // 1D Lagrange - Order > 1 (Higher-order elements)
    template <unsigned Order>
    struct FiniteElementSpaceTraits<LagrangeSpaceTag, 1, Order> {
        // Vertices + internal edge DOFs for higher order
        using EntityTypes = std::tuple<Vertex<1>, EdgeX<1>>;
        using DOFNums = std::tuple<
            std::integral_constant<unsigned, 1>,           // 1 DOF per vertex
            std::integral_constant<unsigned, Order - 1>    // (Order-1) internal DOFs per edge
        >;

        static constexpr unsigned dofsPerElement = 2 * 1 + 1 * (Order - 1); // 2 vertices + 1 edge
    };

    // 2D Lagrange - Order 1 (Bilinear elements)
    template <>
    struct FiniteElementSpaceTraits<LagrangeSpaceTag, 2, 1> {

        // Only vertices have DOFs in first-order elements
        using EntityTypes = std::tuple<Vertex<2>>;
        using DOFNums = std::tuple<std::integral_constant<unsigned, 1>>;

        static constexpr unsigned dofsPerElement = 4; // 4 vertices
    };

    // 2D Lagrange - Order > 1 (Higher-order elements)
    template <unsigned Order>
    struct FiniteElementSpaceTraits<LagrangeSpaceTag, 2, Order> {

        // Vertices + edges + faces for higher order
        using EntityTypes = std::tuple<Vertex<2>, EdgeX<2>, EdgeY<2>, FaceXY<2>>;
        using DOFNums = std::tuple<
            std::integral_constant<unsigned, 1>,                           // 1 DOF per vertex
            std::integral_constant<unsigned, Order - 1>,                   // (Order-1) DOFs per edge
            std::integral_constant<unsigned, Order - 1>,                   // (Order-1) DOFs per edge
            std::integral_constant<unsigned, (Order - 1) * (Order - 1)>   // (Order-1)^2 DOFs per face
        >;

        static constexpr unsigned dofsPerElement = 4 * 1 + 2 * (Order - 1) + 2 * (Order - 1) + 1 * (Order - 1) * (Order - 1);
        // 4 vertices + 2 EdgeX + 2 EdgeY + 1 face
    };

    // 3D Lagrange - Order 1 (Trilinear elements)
    template <>
    struct FiniteElementSpaceTraits<LagrangeSpaceTag, 3, 1> {

        // Only vertices have DOFs in first-order elements
        using EntityTypes = std::tuple<Vertex<3>>;
        using DOFNums = std::tuple<std::integral_constant<unsigned, 1>>;

        static constexpr unsigned dofsPerElement = 8; // 8 vertices
    };

    // 3D Lagrange - Order > 1 (Higher-order elements)
    template <unsigned Order>
    struct FiniteElementSpaceTraits<LagrangeSpaceTag, 3, Order> {

        // All entity types for higher order
        using EntityTypes = std::tuple<
            Vertex<3>, EdgeX<3>, EdgeY<3>, EdgeZ<3>,
            FaceXY<3>, FaceXZ<3>, FaceYZ<3>, Hexahedron<3>
        >;
        using DOFNums = std::tuple<
            std::integral_constant<unsigned, 1>,                                       // 1 DOF per vertex
            std::integral_constant<unsigned, Order - 1>,                               // (Order-1) DOFs per edge
            std::integral_constant<unsigned, Order - 1>,                               // (Order-1) DOFs per edge
            std::integral_constant<unsigned, Order - 1>,                               // (Order-1) DOFs per edge
            std::integral_constant<unsigned, (Order - 1) * (Order - 1)>,              // (Order-1)^2 DOFs per face
            std::integral_constant<unsigned, (Order - 1) * (Order - 1)>,              // (Order-1)^2 DOFs per face
            std::integral_constant<unsigned, (Order - 1) * (Order - 1)>,              // (Order-1)^2 DOFs per face
            std::integral_constant<unsigned, (Order - 1) * (Order - 1) * (Order - 1)> // (Order-1)^3 DOFs per volume
        >;

        static constexpr unsigned dofsPerElement = 8 * 1 + 4 * (Order - 1) + 4 * (Order - 1) + 4 * (Order - 1)
                                                  + 2 * (Order - 1) * (Order - 1) + 2 * (Order - 1) * (Order - 1) + 2 * (Order - 1) * (Order - 1)
                                                  + 1 * (Order - 1) * (Order - 1) * (Order - 1);
        // 8 vertices + 4 EdgeX + 4 EdgeY + 4 EdgeZ + 2 FaceXY + 2 FaceXZ + 2 FaceYZ + 1 hexahedron
    };
    
    ///////////////////////////////////////////////////////////////////////
    // Nédélec Space Specializations /////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    
    // Nédélec elements always start from Order 1 and never have vertex DOFs
    
    // 2D Nédélec - Order 1 (Lowest-order H(curl))
    template <>
    struct FiniteElementSpaceTraits<NedelecSpaceTag, 2, 1> {

        // Only edges have DOFs in first-order Nédélec
        using EntityTypes = std::tuple<EdgeX<2>, EdgeY<2>>;
        using DOFNums = std::tuple<
            std::integral_constant<unsigned, 1>,  // 1 DOF per edge
            std::integral_constant<unsigned, 1>   // 1 DOF per edge
        >;

        static constexpr unsigned dofsPerElement = 2 * 1 + 2 * 1; // 2 EdgeX + 2 EdgeY
    };

    ///////////////////////////////////////////////////////////////////////
    // Convenience Functions and Type Aliases ////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Convenience aliases for creating FEMContainers with specific spaces
     */
    template <typename T, typename SpaceTag, unsigned Dim, unsigned Order>
    using SpaceTraits = FiniteElementSpaceTraits<SpaceTag, Dim, Order>;
    
    template <typename T, typename SpaceTag, unsigned Dim, unsigned Order>
    using SpaceFEMContainer = FEMContainer<T, Dim, 
                                         typename SpaceTraits<T, SpaceTag, Dim, Order>::EntityTypes,
                                         typename SpaceTraits<T, SpaceTag, Dim, Order>::DOFNums>;
    
    // Specific space aliases
    template <typename T, unsigned Dim, unsigned Order>
    using LagrangeFEMContainer = SpaceFEMContainer<T, LagrangeSpaceTag, Dim, Order>;
    
    template <typename T, unsigned Dim, unsigned Order>
    using NedelecFEMContainer = SpaceFEMContainer<T, NedelecSpaceTag, Dim, Order>;

}  // namespace ippl

#endif  // IPPL_FINITE_ELEMENT_SPACE_TRAITS_H