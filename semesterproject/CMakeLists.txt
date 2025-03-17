file (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
message (STATUS "Adding index test found in ${_relPath}")

include_directories (
    ${CMAKE_SOURCE_DIR}/src
)

link_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${Kokkos_DIR}/..
)

set (IPPL_LIBS ippl ${MPI_CXX_LIBRARIES})
set (COMPILE_FLAGS ${OPAL_CXX_FLAGS})

add_executable(FusionReactor FusionReactor.cpp)
target_link_libraries(FusionReactor PRIVATE ${IPPL_LIBS} ${VTK_LIBRARIES})
vtk_module_autoinit( TARGETS FusionReactor MODULES ${VTK_LIBRARIES})

add_executable(TestUnstructuredFieldContainer TestUnstructuredFieldContainer.cpp)
target_link_libraries(TestUnstructuredFieldContainer PRIVATE ${IPPL_LIBS} ${VTK_LIBRARIES})
vtk_module_autoinit( TARGETS TestUnstructuredFieldContainer MODULES ${VTK_LIBRARIES})

add_executable(WeightsTest WeightsTest.cpp)
target_link_libraries(WeightsTest PRIVATE ${IPPL_LIBS} ${VTK_LIBRARIES})
vtk_module_autoinit( TARGETS WeightsTest MODULES ${VTK_LIBRARIES})

# vi: set et ts=4 sw=4 sts=4:

# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# require-final-newline: nil
# End: