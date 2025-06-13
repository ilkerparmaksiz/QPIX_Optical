## ---------------------------------------------------------
##  G4_QPIX | CMakeLists.txt
##
##  Main CMake build script.
##   * Author: Justo Martin-Albo
##   * Creation date: 14 Aug 2019
## ---------------------------------------------------------
#set( ENV{PATH} $ENV{NOCONDA_PATH} )
#MESSAGE("PATH VAR $ENV{PATH}")

cmake_minimum_required(VERSION 3.12)
project(QPixG4)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

## Find Geant4 package, activating all available UI and Vis
## drivers by default. You can set WITH_GEANT4_UIVIS to OFF
## via the command line or ccmake/cmake-gui to build a batch
## mode only executable.
option(WITH_GEANT4_UIVIS "Build with Geant4 UI/VIS drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

option (With_Opticks OFF)



## Setup Geant4 include directories and compile definitions.
include(${Geant4_USE_FILE})

## find the ROOT package
find_package(ROOT REQUIRED)

## include ROOT header files
include(${ROOT_USE_FILE})

## link ROOT libraries
link_libraries(${ROOT_LIBRARIES})

if(With_Opticks)
  message(STATUS "WITH_OPTICKS is set ON")
  option (With_G4OpticksTest OFF)  # This is used to compare simulated photon events between  G4 and Opticks by producing them in the same event.

  if(With_G4OpticksTest)
    message(STATUS "Warning ! With_G4OpticksTest Option is set ON !")
    message(STATUS "With_G4OpticksTest Option allows you to save both G4 and Opticks photon hits into same file so hits can be compared between two ")
  endif ()

  list(APPEND CMAKE_MODULE_PATH "$ENV{OPTICKS_HOME}/cmake/Modules")
  #list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/Modules")
  set(OPTICKS_PREFIX "$ENV{OPTICKS_PREFIX}")
  ### we need to Make sure this OPTICKS_PREFIX matches with the one on OpticksBuildOptions
  message(STATUS "OPTICKS_PREFIX : $ENV{OPTICKS_PREFIX} ")
  include(OpticksBuildOptions)
  message(STATUS "---- OPTICKS ----")
  message(STATUS "CMAKE_MODULE_PATH --> ${CMAKE_MODULE_PATH}")
  message(STATUS "OPTICKS HOME --> $ENV{OPTICKS_HOME}")
  message(STATUS "---- OPTICKS END ----")
  #find_dependency(PLog) # For Some Reason it was not finding PLog , doing this compiles
  ## These are needed for opticks
  find_package(G4CX REQUIRED CONFIG)
  find_package(U4 REQUIRED CONFIG)

endif ()

# include a config file
configure_file(cfg/config.h.in cfg/config.h)

## Recurse through sub-directories
add_subdirectory(src)
add_subdirectory(app)
add_subdirectory(cfg)

