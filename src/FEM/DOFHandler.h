// Class DOFHandler
//   This class is responsible for handling the degrees of freedom (DOFs) in a finite
//   element mesh and space. It provides methods to access and manipulate DOFs,
//   including mapping between local DOF indices and their global position in the FEMContainer.

#ifndef IPPL_DOFHANDLER_H
#define IPPL_DOFHANDLER_H

#include "FEM/FiniteElementSpaceTraits.h"
#include "FEM/FEMContainer.h"
#include "Utility/IpplException.h"
#include <array>

namespace ippl {

    /**
     * @brief DOFHandler maps local element DOFs to entity types and local entity indices
     * 
     * The DOFHandler's main job is to answer: "Given an element and a local DOF number,
     * which entity type does this DOF belong to, and what is the local index within that entity?"
     * 
     * @tparam T The floating point type
     * @tparam SpaceTag The finite element space type (LagrangeSpaceTag, NedelecSpaceTag, etc.)
     * @tparam Dim The spatial dimension  
     * @tparam Order The polynomial order
     */
    template <typename T, typename SpaceTag, unsigned Dim, unsigned Order>
    class DOFHandler {
    public:
        // Space traits
        using SpaceTraits = FiniteElementSpaceTraits<SpaceTag, Dim, Order>;
        using EntityTypes = typename SpaceTraits::EntityTypes;
        using DOFNums = typename SpaceTraits::DOFNums;
        
        // Compatible FEMContainer type
        using FEMContainer_t = FEMContainer<T, Dim, EntityTypes, DOFNums>;
        
        static constexpr unsigned dim = Dim;
        static constexpr unsigned order = Order;
        static constexpr unsigned dofsPerElement = SpaceTraits::dofsPerElement;
        
        // Mesh types
        using Mesh_t = UniformCartesian<T, Dim>;
        using Layout_t = FieldLayout<Dim>;
        using indices_t = Vector<size_t, Dim>;

        /**
         * @brief Structure to hold DOF mapping information
         */
        struct DOFMapping {
            size_t entityTypeIndex;                         // Index in the EntityTypes tuple
            Kokkos::Array<size_t, Dim> entityLocalIndex;    // Offset from the NDIndex of the element to the NDIndex of the DOF (0 or 1 in each dimension)
            size_t entityLocalDOF;                          // Local DOF number within the entity
        };

        ///////////////////////////////////////////////////////////////////////
        // Constructors ///////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////
        
        DOFHandler();
        DOFHandler(Mesh_t& mesh, const Layout_t& layout);
        
        void initialize(Mesh_t& mesh, const Layout_t& layout);

        ///////////////////////////////////////////////////////////////////////
        // Core DOF Mapping Functions ////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////

        /**
         * @brief Map local element DOF to entity information
         * 
         * @param localElementDOF Local DOF number within the element (0 to dofsPerElement-1)
         * @return DOFMapping containing entity type index, entity local DOF, and entity local index
         */
        KOKKOS_FUNCTION DOFMapping getElementDOFMapping(const size_t& localElementDOF) const;

        /**
         * @brief Get element NDIndex from linear element index
         * 
         * @param elementIndex Linear element index  
         * @return NDIndex of the element in the mesh
         */
        KOKKOS_FUNCTION indices_t getElementNDIndex(const size_t& elementIndex) const;

        /**
         * @brief Get entity NDIndex for a specific entity within an element
         * 
         * @tparam EntityType The type of entity
         * @param elementNDIndex NDIndex of the element
         * @param entityLocalIndex Which instance of this entity in the element (e.g., which vertex: 0,1,2,3)
         * @return NDIndex of the entity in the mesh
         */
        template <typename EntityType>
        KOKKOS_FUNCTION indices_t getEntityNDIndex(const indices_t& elementNDIndex, 
                                                  const size_t& entityLocalIndex) const {
            if constexpr (std::is_same_v<EntityType, Vertex<Dim>>) {
                return getVertexNDIndex(elementNDIndex, entityLocalIndex);
            } else if constexpr (std::is_base_of_v<Edge<Dim>, EntityType>) {
                return getEdgeNDIndex<EntityType>(elementNDIndex, entityLocalIndex);
            } else if constexpr (std::is_base_of_v<Face<Dim>, EntityType>) {
                return getFaceNDIndex<EntityType>(elementNDIndex, entityLocalIndex);
            } else if constexpr (std::is_base_of_v<Volume<Dim>, EntityType>) {
                return getVolumeNDIndex<EntityType>(elementNDIndex, entityLocalIndex);
            }
        }

        ///////////////////////////////////////////////////////////////////////
        // FEMContainer Integration ///////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////

        /**
         * @brief Set a DOF value in a FEMContainer using local element DOF numbering
         * 
         * @param container The FEMContainer to update
         * @param elementIndex Linear element index
         * @param localElementDOF Local DOF number within the element (0 to dofsPerElement-1)
         * @param value The value to set
         */
        KOKKOS_FUNCTION void setElementDOF(FEMContainer_t& container,
                                         const size_t& elementIndex,
                                         const size_t& localElementDOF,
                                         const T& value) const {
            auto mapping = getElementDOFMapping(localElementDOF);
            indices_t elementNDIndex = getElementNDIndex(elementIndex);
            
            // Use template recursion to handle different entity types
            setDOFByEntityIndex(container, mapping, elementNDIndex, value, std::make_index_sequence<std::tuple_size_v<EntityTypes>>{});
        }

        /**
         * @brief Get a DOF value from a FEMContainer using local element DOF numbering
         * 
         * @param container The FEMContainer to read from
         * @param elementIndex Linear element index  
         * @param localElementDOF Local DOF number within the element (0 to dofsPerElement-1)
         * @return The DOF value
         */
        KOKKOS_FUNCTION auto getElementDOF(const FEMContainer_t& container,
                                         const size_t& elementIndex,
                                         const size_t& localElementDOF) const {
            auto mapping = getElementDOFMapping(localElementDOF);
            indices_t elementNDIndex = getElementNDIndex(elementIndex);
            
            // Use template recursion to handle different entity types
            return getDOFByEntityIndex(container, mapping, elementNDIndex, std::make_index_sequence<std::tuple_size_v<EntityTypes>>{});
        }

        /**
         * @brief Add to a DOF value in a FEMContainer (for assembly operations)
         * 
         * @param container The FEMContainer to update
         * @param elementIndex Linear element index
         * @param localElementDOF Local DOF number within the element
         * @param value The value to add
         */
        KOKKOS_FUNCTION void addToElementDOF(FEMContainer_t& container,
                                           const size_t& elementIndex,
                                           const size_t& localElementDOF,
                                           const T& value) const {
            auto currentValue = getElementDOF(container, elementIndex, localElementDOF);
            setElementDOF(container, elementIndex, localElementDOF, currentValue + value);
        }

        ///////////////////////////////////////////////////////////////////////
        // Direct Entity Access ///////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////

        /**
         * @brief Set DOF value using direct entity specification
         * 
         * @tparam EntityType The type of entity
         * @param container The FEMContainer to update
         * @param elementIndex Linear element index
         * @param entityLocalIndex Which instance of this entity in the element
         * @param entityLocalDOF Local DOF number within this entity instance
         * @param value The value to set
         */
        template <typename EntityType>
        KOKKOS_FUNCTION void setEntityDOF(FEMContainer_t& container,
                                        const size_t& elementIndex,
                                        const size_t& entityLocalIndex,
                                        const size_t& entityLocalDOF,
                                        const T& value) const {
            static_assert(SpaceTraits::template hasEntity<EntityType>(), 
                         "EntityType not supported by this finite element space");
            
            indices_t elementNDIndex = getElementNDIndex(elementIndex);
            indices_t entityNDIndex = getEntityNDIndex<EntityType>(elementNDIndex, entityLocalIndex);
            
            auto view = container.template getView<EntityType>();
            if constexpr (Dim == 1) {
                view(entityNDIndex[0])[entityLocalDOF] = value;
            } else if constexpr (Dim == 2) {
                view(entityNDIndex[0], entityNDIndex[1])[entityLocalDOF] = value;
            } else if constexpr (Dim == 3) {
                view(entityNDIndex[0], entityNDIndex[1], entityNDIndex[2])[entityLocalDOF] = value;
            }
        }

        /**
         * @brief Get DOF value using direct entity specification
         */
        template <typename EntityType>
        KOKKOS_FUNCTION auto getEntityDOF(const FEMContainer_t& container,
                                        const size_t& elementIndex,
                                        const size_t& entityLocalIndex,
                                        const size_t& entityLocalDOF) const {
            static_assert(SpaceTraits::template hasEntity<EntityType>(), 
                         "EntityType not supported by this finite element space");
            
            indices_t elementNDIndex = getElementNDIndex(elementIndex);
            indices_t entityNDIndex = getEntityNDIndex<EntityType>(elementNDIndex, entityLocalIndex);
            
            auto view = container.template getView<EntityType>();
            if constexpr (Dim == 1) {
                return view(entityNDIndex[0])[entityLocalDOF];
            } else if constexpr (Dim == 2) {
                return view(entityNDIndex[0], entityNDIndex[1])[entityLocalDOF];
            } else if constexpr (Dim == 3) {
                return view(entityNDIndex[0], entityNDIndex[1], entityNDIndex[2])[entityLocalDOF];
            }
        }

    private:
        ///////////////////////////////////////////////////////////////////////
        // Member Variables ///////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////
        
        Mesh_t* mesh_m;
        const Layout_t* layout_m;
        
        // Number of elements in each direction
        Vector<size_t, Dim> ne_m;

        // DOF mapping table (dual view for host and device access)
        Kokkos::DualView<DOFMapping*> dofMappingTable_m;

        ///////////////////////////////////////////////////////////////////////
        // Space-Specific DOF Mapping Implementation //////////////////////////
        ///////////////////////////////////////////////////////////////////////

        /**
         * @brief Compute DOF mapping for Lagrange elements on-the-fly
         */
        KOKKOS_FUNCTION DOFMapping computeLagrangeDOFMapping(const size_t& localElementDOF) const;

        /**
         * @brief Compute DOF mapping for Nédélec elements on-the-fly
         */
        KOKKOS_FUNCTION DOFMapping computeNedelecDOFMapping(const size_t& localElementDOF) const;

        ///////////////////////////////////////////////////////////////////////
        // Entity NDIndex Computation /////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////

        /**
         * @brief Get vertex NDIndex within an element
         */
        KOKKOS_FUNCTION indices_t getVertexNDIndex(const indices_t& elementNDIndex, 
                                                  const size_t& vertexLocalIndex) const {
            indices_t vertexNDIndex = elementNDIndex;
            
            // Convert vertex local index to offset
            if constexpr (Dim == 1) {
                vertexNDIndex[0] += vertexLocalIndex; // 0 or 1
            } else if constexpr (Dim == 2) {
                vertexNDIndex[0] += vertexLocalIndex & 1;        // 0,1,0,1 for vertices 0,1,2,3
                vertexNDIndex[1] += (vertexLocalIndex >> 1) & 1; // 0,0,1,1 for vertices 0,1,2,3
            } else if constexpr (Dim == 3) {
                vertexNDIndex[0] += vertexLocalIndex & 1;        // 0,1,0,1,0,1,0,1 for vertices 0-7
                vertexNDIndex[1] += (vertexLocalIndex >> 1) & 1; // 0,0,1,1,0,0,1,1 for vertices 0-7
                vertexNDIndex[2] += (vertexLocalIndex >> 2) & 1; // 0,0,0,0,1,1,1,1 for vertices 0-7
            }
            
            return vertexNDIndex;
        }

        /**
         * @brief Get edge NDIndex within an element
         */
        template <typename EdgeType>
        KOKKOS_FUNCTION indices_t getEdgeNDIndex(const indices_t& elementNDIndex,
                                               const size_t& edgeLocalIndex) const {
            indices_t edgeNDIndex = elementNDIndex;
            
            // Edge indexing depends on edge orientation and local index
            if constexpr (std::is_same_v<EdgeType, EdgeX<Dim>>) {
                // X-oriented edges
                if constexpr (Dim == 2) {
                    edgeNDIndex[1] += edgeLocalIndex; // Bottom (0) or top (1) edge
                } else if constexpr (Dim == 3) {
                    edgeNDIndex[1] += (edgeLocalIndex >> 1) & 1;
                    edgeNDIndex[2] += edgeLocalIndex & 1;
                }
            } else if constexpr (std::is_same_v<EdgeType, EdgeY<Dim>>) {
                // Y-oriented edges
                if constexpr (Dim == 2) {
                    edgeNDIndex[0] += edgeLocalIndex; // Left (0) or right (1) edge
                } else if constexpr (Dim == 3) {
                    edgeNDIndex[0] += (edgeLocalIndex >> 1) & 1;
                    edgeNDIndex[2] += edgeLocalIndex & 1;
                }
            } else if constexpr (std::is_same_v<EdgeType, EdgeZ<Dim>> && Dim == 3) {
                // Z-oriented edges (3D only)
                edgeNDIndex[0] += (edgeLocalIndex >> 1) & 1;
                edgeNDIndex[1] += edgeLocalIndex & 1;
            }
            
            return edgeNDIndex;
        }

        /**
         * @brief Get face NDIndex within an element
         */
        template <typename FaceType>
        KOKKOS_FUNCTION indices_t getFaceNDIndex(const indices_t& elementNDIndex,
                                               const size_t& faceLocalIndex) const {
            indices_t faceNDIndex = elementNDIndex;
            
            // Face indexing depends on face orientation
            if constexpr (std::is_same_v<FaceType, FaceXY<Dim>>) {
                if constexpr (Dim == 2) {
                    // In 2D, there's only one face per element (the element itself)
                    // faceLocalIndex should be 0
                } else if constexpr (Dim == 3) {
                    faceNDIndex[2] += faceLocalIndex; // Bottom (0) or top (1) face
                }
            } else if constexpr (std::is_same_v<FaceType, FaceXZ<Dim>> && Dim == 3) {
                // XZ-plane faces (3D only)
                faceNDIndex[1] += faceLocalIndex; // Front (0) or back (1) face
            } else if constexpr (std::is_same_v<FaceType, FaceYZ<Dim>> && Dim == 3) {
                // YZ-plane faces (3D only)
                faceNDIndex[0] += faceLocalIndex; // Left (0) or right (1) face
            }
            
            return faceNDIndex;
        }

        /**
         * @brief Get volume NDIndex within an element
         */
        template <typename VolumeType>
        KOKKOS_FUNCTION indices_t getVolumeNDIndex(const indices_t& elementNDIndex,
                                                 const size_t& volumeLocalIndex) const {
            // Volume entities are the element itself, so return element NDIndex
            return elementNDIndex;
        }

        ///////////////////////////////////////////////////////////////////////
        // Template Recursion Helpers /////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////

        template <typename ValueType, size_t... Is>
        KOKKOS_FUNCTION void setDOFByEntityIndex(FEMContainer_t& container,
                                               const DOFMapping& mapping,
                                               const indices_t& elementNDIndex,
                                               const ValueType& value,
                                               std::index_sequence<Is...>) const {
            ((Is == mapping.entityTypeIndex ? 
              setDOFForEntity<std::tuple_element_t<Is, EntityTypes>>(container, mapping, elementNDIndex, value) : 
              void()), ...);
        }

        template <typename ValueType, size_t... Is>
        KOKKOS_FUNCTION auto getDOFByEntityIndex(const FEMContainer_t& container,
                                               const DOFMapping& mapping,
                                               const indices_t& elementNDIndex,
                                               std::index_sequence<Is...>) const {
            // This is tricky with template recursion - you might need a different approach
            // For now, return a default-constructed value
            using ReturnType = typename std::tuple_element_t<0, EntityTypes>::value_type;
            ReturnType result{};
            
            ((Is == mapping.entityTypeIndex ? 
              (result = getDOFForEntity<std::tuple_element_t<Is, EntityTypes>>(container, mapping, elementNDIndex)) : 
              void()), ...);
              
            return result;
        }

        template <typename EntityType, typename ValueType>
        KOKKOS_FUNCTION void setDOFForEntity(FEMContainer_t& container,
                                           const DOFMapping& mapping,
                                           const indices_t& elementNDIndex,
                                           const ValueType& value) const {
            // Compute entity NDIndex from element NDIndex + offset
            indices_t entityNDIndex;
            for (unsigned d = 0; d < Dim; ++d) {
                entityNDIndex[d] = elementNDIndex[d] + mapping.entityLocalIndex[d];
            }

            auto view = container.template getView<EntityType>();
            if constexpr (Dim == 1) {
                view(entityNDIndex[0])[mapping.entityLocalDOF] = value;
            } else if constexpr (Dim == 2) {
                view(entityNDIndex[0], entityNDIndex[1])[mapping.entityLocalDOF] = value;
            } else if constexpr (Dim == 3) {
                view(entityNDIndex[0], entityNDIndex[1], entityNDIndex[2])[mapping.entityLocalDOF] = value;
            }
        }

        template <typename EntityType>
        KOKKOS_FUNCTION auto getDOFForEntity(const FEMContainer_t& container,
                                           const DOFMapping& mapping,
                                           const indices_t& elementNDIndex) const {
            // Compute entity NDIndex from element NDIndex + offset
            indices_t entityNDIndex;
            for (unsigned d = 0; d < Dim; ++d) {
                entityNDIndex[d] = elementNDIndex[d] + mapping.entityLocalIndex[d];
            }

            auto view = container.template getView<EntityType>();
            if constexpr (Dim == 1) {
                return view(entityNDIndex[0])[mapping.entityLocalDOF];
            } else if constexpr (Dim == 2) {
                return view(entityNDIndex[0], entityNDIndex[1])[mapping.entityLocalDOF];
            } else if constexpr (Dim == 3) {
                return view(entityNDIndex[0], entityNDIndex[1], entityNDIndex[2])[mapping.entityLocalDOF];
            }
        }
    };

    ///////////////////////////////////////////////////////////////////////
    // Type Aliases for Convenience ///////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////

    template <typename T, unsigned Dim, unsigned Order>
    using LagrangeDOFHandler = DOFHandler<T, LagrangeSpaceTag, Dim, Order>;

    template <typename T, unsigned Dim, unsigned Order>
    using NedelecDOFHandler = DOFHandler<T, NedelecSpaceTag, Dim, Order>;

}  // namespace ippl

#include "DOFHandler.hpp"

#endif  // IPPL_DOFHANDLER_H