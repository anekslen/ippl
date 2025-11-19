# Install script for directory: /home/runner/work/ippl/ippl/src/Communicate

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Communicate" TYPE FILE FILES
    "/home/runner/work/ippl/ippl/src/Communicate/LogEntry.h"
    "/home/runner/work/ippl/ippl/src/Communicate/BufferHandler.h"
    "/home/runner/work/ippl/ippl/src/Communicate/BufferHandler.hpp"
    "/home/runner/work/ippl/ippl/src/Communicate/LoggingBufferHandler.h"
    "/home/runner/work/ippl/ippl/src/Communicate/LoggingBufferHandler.hpp"
    "/home/runner/work/ippl/ippl/src/Communicate/Archive.h"
    "/home/runner/work/ippl/ippl/src/Communicate/Archive.hpp"
    "/home/runner/work/ippl/ippl/src/Communicate/Buffers.hpp"
    "/home/runner/work/ippl/ippl/src/Communicate/Communicator.h"
    "/home/runner/work/ippl/ippl/src/Communicate/Environment.h"
    "/home/runner/work/ippl/ippl/src/Communicate/DataTypes.h"
    "/home/runner/work/ippl/ippl/src/Communicate/Operations.h"
    "/home/runner/work/ippl/ippl/src/Communicate/Collectives.hpp"
    "/home/runner/work/ippl/ippl/src/Communicate/Serializable.h"
    "/home/runner/work/ippl/ippl/src/Communicate/Request.h"
    "/home/runner/work/ippl/ippl/src/Communicate/Status.h"
    "/home/runner/work/ippl/ippl/src/Communicate/TagMaker.h"
    "/home/runner/work/ippl/ippl/src/Communicate/Tags.h"
    "/home/runner/work/ippl/ippl/src/Communicate/Wait.h"
    "/home/runner/work/ippl/ippl/src/Communicate/Window.h"
    "/home/runner/work/ippl/ippl/src/Communicate/Window.hpp"
    "/home/runner/work/ippl/ippl/src/Communicate/PointToPoint.hpp"
    "/home/runner/work/ippl/ippl/src/Communicate/CommunicatorLogging.hpp"
    )
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/home/runner/work/ippl/ippl/_codeql_build_dir/src/Communicate/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
