cmake_minimum_required(VERSION 3.1)

project(rocfoam CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set project metadata #########################################################

set(rocfoam_MAJOR_VERSION 0)
set(rocfoam_MINOR_VERSION 1)
set(rocfoam_PATCH_VERSION 0)
set(rocfoam_EXTRA_VERSION "")
set(rocfoam_VERSION "${rocfoam_MAJOR_VERSION}.${rocfoam_MINOR_VERSION}")
set(rocfoam_VERSION
    "${rocfoam_VERSION}.${rocfoam_PATCH_VERSION}${rocfoam_EXTRA_VERSION}")

# package module
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/Modules/)

# Finding libraries ############################################################

# Find OpenFOAM
# checking OpenFOAM version
string(FIND $ENV{WM_PROJECT_VERSION} "4" OF4)
if(OF4 GREATER -1)
  # OpenFOAM 4.x
  # TODO: rmains to be tested on OF4
  find_package(OpenFOAM REQUIRED COMPONENTS OpenFOAM finiteVolume fvOptions
      compressibleTransportModels fluidThermophysicalModels specie
      turbulenceModels compressibleTurbulenceModels dynamicFvMesh
      topoChangerFvMesh meshTools)
else()
  # > 4.x series
  find_package(OpenFOAM REQUIRED COMPONENTS OpenFOAM finiteVolume fvOptions
      compressibleTransportModels fluidThermophysicalModels specie
      turbulenceModels compressibleTurbulenceModels dynamicFvMesh
      topoChangerFvMesh meshTools)
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
    src/rhoCentral/rocRhoCentral.C
    src/rhoPimple/rocRhoPimple.C)

# Setting rhocentral driver source
set(rocCntDrv_SRCS src/rhoCentral/rocRhoCentral_drv.C)

# Setting rhocentral driver source
set(rocPmpDrv_SRCS src/rhoPimple/rocRhoPimple_drv.C)

# Building library #############################################################

# Build rocfoam library
add_library(rocfoam ${rocfoam_SRCS})
target_include_directories(rocfoam
    PUBLIC
        $<INSTALL_INTERFACE:include/rocfoam>
)
foreach(include_dir ${rocfoam_INC_DIRS})
  target_include_directories(rocfoam
      PUBLIC
        $<BUILD_INTERFACE:${include_dir}>
)
endforeach()
target_compile_options(rocfoam
    PRIVATE
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
    PRIVATE
        ${OPNF_LIBRARIES}
)
target_include_directories(rocfoam SYSTEM
    PRIVATE
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
    DESTINATION include/rocfoam)

# Libraries
install(TARGETS rocfoam
    EXPORT rocfoam
    RUNTIME DESTINATION bin
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)


# Building utilities #############################################################

add_executable(rocRhoCentral ${rocCntDrv_SRCS})
target_link_libraries(rocRhoCentral rocfoam)
target_include_directories(rocRhoCentral SYSTEM
    PRIVATE
        ${OPNF_INCLUDE_DIRS}
)

add_executable(rocRhoPimple ${rocPmpDrv_SRCS})
target_link_libraries(rocRhoPimple rocfoam)
target_include_directories(rocRhoPimple SYSTEM
    PRIVATE
        ${OPNF_INCLUDE_DIRS}
)

install(TARGETS rocRhoCentral rocRhoPimple
         RUNTIME DESTINATION bin
         PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
         GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
