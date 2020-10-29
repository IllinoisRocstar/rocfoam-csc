cmake_minimum_required(VERSION 3.1)

# Set project metadata #########################################################

project(RocFOAM CXX)

set(RocFOAM_MAJOR_VERSION 1)
set(RocFOAM_MINOR_VERSION 0)
set(RocFOAM_PATCH_VERSION 2)
set(RocFOAM_EXTRA_VERSION "")
set(RocFOAM_VERSION "${RocFOAM_MAJOR_VERSION}.${RocFOAM_MINOR_VERSION}")
set(RocFOAM_VERSION "${RocFOAM_VERSION}.${RocFOAM_PATCH_VERSION}${RocFOAM_EXTRA_VERSION}")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set project defaults #########################################################

# package module
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/Modules/)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
endif()
# Set the possible values of build type for cmake-gui
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")

# creating code compilation database
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Installation Directories
include(GNUInstallDirs)

# Setting configure options ####################################################

option(ENABLE_TESTING "Enable Testing of RocFOAM." ON)
option(ENABLE_VERBOSITY "Enable Verbosity of RocFOAM." OFF)

# Setting compile and link options #############################################

# Turn on default warnings for GNU compilers.
if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  # Update if necessary
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wpedantic")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Wno-unused-parameter")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-invalid-offsetof -Wno-attributes")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wold-style-cast -Wnon-virtual-dtor")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wstrict-null-sentinel -Wredundant-decls")
  # -Wshadow gives false positives for gcc versions <= 4.8
  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.9)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wshadow")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wfloat-conversion")
  endif(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.9)
endif()

# set up RPATH
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Finding libraries ############################################################

# Find OpenFOAM
# checking OpenFOAM version
string(FIND $ENV{WM_PROJECT_VERSION} "2006" OF20)
if(OF20 GREATER -1)
    set(OF_VERSION 20)
endif()

string(FIND $ENV{WM_PROJECT_VERSION} "7" OF7)
if(OF7 GREATER -1)
    set(OF_VERSION 7)
endif()

string(FIND $ENV{WM_PROJECT_VERSION} "8" OF8)
if(OF8 GREATER -1)
    set(OF_VERSION 8)
endif()

message(STATUS "OpenFOAM version ${OF_VERSION}")
if(OF_VERSION EQUAL 20)
  find_package(OpenFOAM REQUIRED COMPONENTS OpenFOAM finiteVolume fvOptions
      compressibleTransportModels
      fluidThermophysicalModels specie
      turbulenceModels compressibleTurbulenceModels
      dynamicFvMesh
      topoChangerFvMesh meshTools Pstream)
elseif(OF_VERSION EQUAL 7)
  find_package(OpenFOAM REQUIRED COMPONENTS OpenFOAM finiteVolume fvOptions
      compressibleTransportModels
      fluidThermophysicalModels specie
      turbulenceModels compressibleTurbulenceModels
      dynamicFvMesh
      topoChangerFvMesh meshTools Pstream)
elseif(OF_VERSION EQUAL 8)
  find_package(OpenFOAM REQUIRED COMPONENTS OpenFOAM finiteVolume fvOptions
      momentumTransportModels fluidThermoMomentumTransportModels
      thermophysicalTransportModels
      fluidThermophysicalModels specie
      dynamicFvMesh sampling
      topoChangerFvMesh meshTools Pstream)
else()
  message(STATUS "OpenFOAM version not identified/supported")
  #  # > 7.x series
  #  find_package(OpenFOAM REQUIRED COMPONENTS OpenFOAM finiteVolume fvOptions
  #      compressibleTransportModels
  #      compressibleTurbulenceModels
  #      fluidThermophysicalModels specie
  #      turbulenceModels  dynamicFvMesh
  #      topoChangerFvMesh meshTools Pstream)
endif()
message(STATUS "OpenFOAM include directory ${OPNF_INCLUDE_DIRS}")
message(STATUS "OpenFOAM link directory ${OPNF_LIBRARY_DIRS}")

# Find IMPACT
find_package(IMPACT REQUIRED)

# Find Boost
#set(Boost_USE_MULTITHREADED ON)
find_package(Boost 1.45.0 REQUIRED COMPONENTS filesystem)

# Setting library files ########################################################

# Setting include files
file(GLOB_RECURSE RocFOAM_INC_FILES *.H)

set(RocFOAM_INC_DIRS "")
foreach(file_path ${RocFOAM_INC_FILES})
  get_filename_component(dir_path ${file_path} PATH)
  list(APPEND RocFOAM_INC_DIRS ${dir_path})
endforeach()

list(REMOVE_DUPLICATES RocFOAM_INC_DIRS)

# Setting RocFOAM source
set(RocFOAM_SRCS
    src/foam/rocFoam.C
    src/foam/comFoam.C
    src/foam/fileContainer.C
    src/foam/surfaceDisplacementPointPatchVectorFieldFSI.C

    src/foam/statusMethods.C
    src/foam/volumeMethods.C
    src/foam/faceMethods.C
    src/foam/surfaceMethods.C
    src/foam/reconstMethods.C
    src/foam/filesMethods.C
    )
    
set(rocrhopimple_SRCS
    ${RocFOAM_SRCS}
    src/rhoPimple/rocRhoPimple.C)

# Setting rhoPimple driver source
set(rocPmpDrv_SRCS src/rhoPimple/rocRhoPimple_drv.C)

# Setting rocrhocentral source
set(rocrhocentral_SRCS
    ${RocFOAM_SRCS}
    src/rhoCentral/rocRhoCentral.C)
    
# Setting rhoCentral driver source
set(rocCntDrv_SRCS src/rhoCentral/rocRhoCentral_drv.C)

# Setting rhoFoam driver source
set(rocFoamDrv_SRCS src/foam/rocFoam_drv.C)

# Building library #############################################################

# Build rocrhopimple library ^^^^^^^^^^^^^^^^^^^^^^^^
add_library(rocrhopimple ${rocrhopimple_SRCS})
target_include_directories(rocrhopimple
    PUBLIC
        Boost::headers
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/rocfoam>
)

foreach(include_dir ${RocFOAM_INC_DIRS})
  target_include_directories(rocrhopimple
      PUBLIC
          $<BUILD_INTERFACE:${include_dir}>
  )
endforeach()

target_compile_options(rocrhopimple
    PUBLIC
        -Wall
        -Wextra
        -Wold-style-cast
        -Wnon-virtual-dtor
        -Wno-unused-parameter
        -Wno-invalid-offsetof
    PUBLIC
        -ftemplate-depth-100
)

# Link OpenFOAM
target_link_libraries(rocrhopimple
    PUBLIC
        ${OPNF_LIBRARIES}
        Boost::filesystem
        IMPACT::SITCOM
)

target_include_directories(rocrhopimple SYSTEM
    PUBLIC
        ${OPNF_INCLUDE_DIRS}
)

target_compile_definitions(rocrhopimple
    PUBLIC
        ${OPNF_COMPILE_DEFINITIONS}
)

if(ENABLE_VERBOSITY)
  target_compile_definitions(rocrhopimple PUBLIC VERBOSE)
endif()

if(OF_VERSION EQUAL 20)
  target_compile_definitions(rocrhopimple PUBLIC HAVE_OFE20)
elseif(OF_VERSION EQUAL 7)
  target_compile_definitions(rocrhopimple PUBLIC HAVE_OF7)
elseif(OF_VERSION EQUAL 8)
  target_compile_definitions(rocrhopimple PUBLIC HAVE_OF8)
endif()
message(STATUS "OPENFOAM VERSION ${OF_VERSION}")
#----------------------------------------------------

# Build rocRhoCentral library ^^^^^^^^^^^^^^^^^^^^^^^
add_library(rocrhocentral ${rocrhocentral_SRCS})
target_include_directories(rocrhocentral
    PUBLIC
        Boost::headers
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/rocfoam>
)

foreach(include_dir ${RocFOAM_INC_DIRS})
  target_include_directories(rocrhocentral
      PUBLIC
        $<BUILD_INTERFACE:${include_dir}>
  )
endforeach()

target_compile_options(rocrhocentral
    PUBLIC
        -Wall
        -Wextra
        -Wold-style-cast
        -Wnon-virtual-dtor
        -Wno-unused-parameter
        -Wno-invalid-offsetof
    PUBLIC
        -ftemplate-depth-100
)

# Link OpenFOAM
target_link_libraries(rocrhocentral
    PUBLIC
        ${OPNF_LIBRARIES}
        Boost::filesystem
        IMPACT::SITCOM
)

target_include_directories(rocrhocentral SYSTEM
    PUBLIC
        ${OPNF_INCLUDE_DIRS}
)

target_compile_definitions(rocrhocentral
    PUBLIC
        ${OPNF_COMPILE_DEFINITIONS}
)

if(ENABLE_VERBOSITY)
  target_compile_definitions(rocrhocentral PUBLIC VERBOSE)
endif()

if(OF_VERSION EQUAL 20)
  target_compile_definitions(rocrhocentral PUBLIC HAVE_OFE20)
elseif(OF_VERSION EQUAL 7)
  target_compile_definitions(rocrhocentral PUBLIC HAVE_OF7)
elseif(OF_VERSION EQUAL 8)
  target_compile_definitions(rocrhocentral PUBLIC HAVE_OF8)
endif()
message(STATUS "OPENFOAM VERSION ${OF_VERSION}")
#----------------------------------------------------

# Set soname version ###########################################################

# Set soname version
message(STATUS "RocFOAM_VERSION: ${RocFOAM_VERSION}")
message(STATUS "RocFOAM_MAJOR_VERSION: ${RocFOAM_MAJOR_VERSION}")
set_target_properties(rocrhopimple PROPERTIES VERSION ${RocFOAM_VERSION}
    SOVERSION ${RocFOAM_MAJOR_VERSION})
# Set soname version
set_target_properties(rocrhocentral PROPERTIES VERSION ${RocFOAM_VERSION}
    SOVERSION ${RocFOAM_MAJOR_VERSION})

# Install libraries ############################################################

# rocrhopimple ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Headers
install(FILES ${RocFOAM_INC_FILES}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Libraries
install(TARGETS rocrhopimple
    EXPORT rocrhopimple
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
#------------------------------------------------

# rocrhocentral ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Headers
install(FILES ${RocFOAM_INC_FILES}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Libraries
install(TARGETS rocrhocentral
    EXPORT rocrhocentral
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
#------------------------------------------------

# Building utilities ###########################################################
add_executable(rocRhoCentral ${rocCntDrv_SRCS})
target_link_libraries(rocRhoCentral rocrhocentral IMPACT::SITCOM)
target_include_directories(rocRhoCentral SYSTEM
    PUBLIC
        ${OPNF_INCLUDE_DIRS}
)

add_executable(rocRhoPimple ${rocPmpDrv_SRCS})
target_link_libraries(rocRhoPimple rocrhopimple IMPACT::SITCOM)
target_include_directories(rocRhoPimple SYSTEM
    PRIVATE
        ${OPNF_INCLUDE_DIRS}
)

##add_executable(rocFoam ${rocFoamDrv_SRCS})
##target_link_libraries(rocFoam rocfoam IMPACT::SITCOM)
##target_include_directories(rocFoam SYSTEM
##    PRIVATE
##        ${OPNF_INCLUDE_DIRS}
##)

install(TARGETS rocRhoCentral rocRhoPimple
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)

# Perform Testing ##############################################################

if(ENABLE_TESTING)
  enable_testing()
  add_subdirectory(testing)
endif()
