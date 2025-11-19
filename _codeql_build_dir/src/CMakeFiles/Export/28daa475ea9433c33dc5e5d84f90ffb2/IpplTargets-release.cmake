#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Ippl::ippl" for configuration "Release"
set_property(TARGET Ippl::ippl APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Ippl::ippl PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libippl.a"
  )

list(APPEND _cmake_import_check_targets Ippl::ippl )
list(APPEND _cmake_import_check_files_for_Ippl::ippl "${_IMPORT_PREFIX}/lib/libippl.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
