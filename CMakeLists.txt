cmake_minimum_required(VERSION 3.1)

project(rocfoam CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set project metadata #########################################################
set(rocfoam_MAJOR_VERSION 1)
set(rocfoam_MINOR_VERSION 0)
set(rocfoam_PATCH_VERSION 0)
set(rocfoam_EXTRA_VERSION "")
set(rocfoam_VERSION "${rocfoam_MAJOR_VERSION}.${rocfoam_MINOR_VERSION}")
set(rocfoam_VERSION
    "${rocfoam_VERSION}.${rocfoam_PATCH_VERSION}${rocfoam_EXTRA_VERSION}")

# package module
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/Modules/)

# Finding libraries ############################################################

option(ENABLE_TESTS "Enabling Testing of RocFOAM" OFF)

if(ENABLE_TESTS)
   enable_testing()
endif()

# Find OpenFOAM
# checking OpenFOAM version
string(FIND $ENV{WM_PROJECT_VERSION} "7" OF7)
if(OF7 GREATER -1)
    set(OF_VERSION 7)
##    target_compile_definitions(rocfoam PUBLIC HAVE_OF7)
endif()

string(FIND $ENV{WM_PROJECT_VERSION} "8" OF8)
if(OF8 GREATER -1)
    set(OF_VERSION 8)
##    target_compile_definitions(rocfoam PUBLIC HAVE_OF8)
endif()

message(STATUS "OpenFOAM version ${OF_VERSION}")
if(OF_VERSION EQUAL 7)
    find_package(OpenFOAM REQUIRED COMPONENTS 
                compressibleTransportModels
                fluidThermophysicalModels
                specie
                turbulenceModels
                compressibleTurbulenceModels
                finiteVolume
                dynamicFvMesh
                topoChangerFvMesh
                meshTools
                sampling
                fvOptions
                Pstream)
elseif(OF_VERSION EQUAL 8)
    find_package(OpenFOAM REQUIRED COMPONENTS 
                fluidThermophysicalModels
                specie
                momentumTransportModels
                fluidThermoMomentumTransportModels
                thermophysicalTransportModels
                finiteVolume
                dynamicFvMesh
                topoChangerFvMesh
                meshTools
                sampling
                fvOptions
                Pstream)
else()
  # > 7.x series
  find_package(OpenFOAM REQUIRED COMPONENTS OpenFOAM finiteVolume fvOptions
      compressibleTransportModels
      compressibleTurbulenceModels
      fluidThermophysicalModels specie
      turbulenceModels  dynamicFvMesh
      topoChangerFvMesh meshTools Pstream)
endif()
message(STATUS "OpenFOAM include directory ${OPNF_INCLUDE_DIRS}")
message(STATUS "OpenFOAM link directory ${OPNF_LIBRARY_DIRS}")



#find_package(IMPACT COMPONENTS)
find_package(IMPACT REQUIRED)

#set(Boost_USE_MULTITHREADED ON)
find_package(Boost 1.45.0 REQUIRED COMPONENTS filesystem)
if(Boost_FOUND)
   include_directories(${Boost_INCLUDE_DIR})
   message(STATUS "Boost found.")
endif()

# Setting library files ########################################################
# Setting include files
file(GLOB_RECURSE rocfoam_INC_FILES *.H)

set(rocfoam_INC_DIRS "")
foreach(file_path ${rocfoam_INC_FILES})
  get_filename_component(dir_path ${file_path} PATH)
  list(APPEND rocfoam_INC_DIRS ${dir_path})
endforeach()

list(REMOVE_DUPLICATES rocfoam_INC_DIRS)

# Setting cfMesh source
set(rocfoam_SRCS
    src/foam/rocFoam.C
    src/foam/comFoam.C
    src/foam/fileContainer.C
    src/foam/surfaceDisplacementPointPatchVectorFieldFSI.C
    src/rhoCentral/rocRhoCentral.C
    src/rhoPimple/rocRhoPimple.C)
    
# Setting rhoCentral driver source
set(rocCntDrv_SRCS src/rhoCentral/rocRhoCentral_drv.C)

# Setting rhoPimple driver source
set(rocPmpDrv_SRCS src/rhoPimple/rocRhoPimple_drv.C)

# Setting rhoFoam driver source
set(rocFoamDrv_SRCS src/foam/rocFoam_drv.C)

# Building library #############################################################

# Build rocfoam library
add_library(rocfoam ${rocfoam_SRCS})
target_include_directories(rocfoam
    PUBLIC
        $<INSTALL_INTERFACE:include/rocfoam>
        Boost::headers
)

foreach(include_dir ${rocfoam_INC_DIRS})
  target_include_directories(rocfoam
      PUBLIC
        $<BUILD_INTERFACE:${include_dir}>
  )
endforeach()

target_compile_options(rocfoam
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
target_link_libraries(rocfoam
    PUBLIC
        ${OPNF_LIBRARIES}
        Boost::filesystem
        IMPACT::SITCOM
)

target_include_directories(rocfoam SYSTEM
    PUBLIC
        ${OPNF_INCLUDE_DIRS}
)

target_compile_definitions(rocfoam
    PUBLIC
        ${OPNF_COMPILE_DEFINITIONS}
)

# Set soname version ###########################################################

set_target_properties(rocfoam PROPERTIES VERSION ${rocfoam_VERSION}
                                        SOVERSION ${rocfoam_MAJOR_VERSION})
                                        
# Install libraries ############################################################

# Headers
install(FILES ${rocfoam_INC_FILES}
    DESTINATION include)

# Libraries
install(TARGETS rocfoam
    EXPORT rocfoam
    RUNTIME DESTINATION bin
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)


if(OF_VERSION EQUAL 7)
    target_compile_definitions(rocfoam PUBLIC HAVE_OF7)
message(STATUS "OPENFOAM VERSION ${OF_VERSION}")
endif() 

if(OF_VERSION EQUAL 8)
    target_compile_definitions(rocfoam PUBLIC HAVE_OF8)

message(STATUS "OPENFOAM VERSION ${OF_VERSION}")

endif()

# Building utilities #############################################################
#add_executable(rocRhoCentral ${rocCntDrv_SRCS})
#target_link_libraries(rocRhoCentral rocfoam)
#target_include_directories(rocRhoCentral SYSTEM
#    PUBLIC
#        ${OPNF_INCLUDE_DIRS}
#)

add_executable(rocRhoPimple ${rocPmpDrv_SRCS})
target_link_libraries(rocRhoPimple rocfoam)
target_include_directories(rocRhoPimple SYSTEM
    PRIVATE
        ${OPNF_INCLUDE_DIRS}
)

add_executable(rocFoam ${rocFoamDrv_SRCS})
target_link_libraries(rocFoam rocfoam IMPACT::SITCOM)
target_include_directories(rocFoam SYSTEM
    PRIVATE
        ${OPNF_INCLUDE_DIRS}
)

#install(TARGETS rocRhoCentral rocRhoPimple rocFoam
#         RUNTIME DESTINATION bin
#         PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
#         GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE

install(TARGETS rocRhoPimple rocFoam
         RUNTIME DESTINATION bin
         PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
         GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)

if(ENABLE_TESTS)
   add_subdirectory(testing)
endif()

