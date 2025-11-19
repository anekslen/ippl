# -----------------------------------------------------------------------------
# src/Field/CMakeLists.txt
#
# Adds core Field data structure headers to the IPPL target.
# -----------------------------------------------------------------------------

# Add this directory to public include path
target_include_directories(ippl
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
)

# Install Field-related headers
install(FILES
    # Elements
    Elements/Element.h
    Elements/HexahedralElement.h
    Elements/HexahedralElement.hpp
    Elements/EdgeElement.h
    Elements/EdgeElement.hpp
    Elements/QuadrilateralElement.h
    Elements/QuadrilateralElement.hpp

    # Quadrature
    Quadrature/GaussJacobiQuadrature.h
    Quadrature/GaussJacobiQuadrature.hpp
    Quadrature/MidpointQuadrature.h
    Quadrature/MidpointQuadrature.hpp
    Quadrature/Quadrature.h
    Quadrature/Quadrature.hpp

    # Spaces
    FiniteElementSpace.h
    FiniteElementSpace.hpp
    LagrangeSpaceFEMContainer.h
    LagrangeSpaceFEMContainer.hpp
    # LagrangeSpace.h
    # LagrangeSpace.hpp

    # Degrees of Freedom Container
    FEMContainer.h
    FEMContainer.hpp

    # Degrees of Freedom Handler
    DOFHandler.h
    DOFHandler.hpp

    # Entity and helper Structs
    Entity.h
    DOFArray.h
    FEMHelperStructs.h

    # SpaceTraits
    FiniteElementSpaceTraits.h

    DESTINATION include/FEM
)

