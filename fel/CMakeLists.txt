# -----------------------------------------------------------------------------
# fel/CMakeLists.txt
#
# This module provides the implementation of a free electron laser in ippl.
# -----------------------------------------------------------------------------

set(COMMON_LIBS
    ippl
    ${MPI_CXX_LIBRARIES}
)

set(COMMON_FLAGS
    ${OPAL_CXX_FLAGS}
)

function(add_fel_example name)
    add_executable(${name} ${name}.cpp)

    target_compile_options(${name} PRIVATE ${COMMON_FLAGS})

    target_link_libraries(${name}
        PRIVATE ${COMMON_LIBS}
    )

    target_include_directories(${name}
        PRIVATE
            $<TARGET_PROPERTY:ippl,INTERFACE_INCLUDE_DIRECTORIES>
            ${CMAKE_CURRENT_SOURCE_DIR}/..
            ${DOWNLOADED_HEADERS_DIR}
    )
endfunction()

add_fel_example(FreeElectronLaser)
add_fel_example(TestFDTDSolverWithParticles)