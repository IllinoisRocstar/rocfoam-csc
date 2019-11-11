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

# Find OpenFOAM
# checking OpenFOAM version
string(FIND $ENV{WM_PROJECT_VERSION} "7" OF7)
if(OF7 GREATER -1)
  # OpenFOAM 7.x
  # TODO: rmains to be tested on OF7
  find_package(OpenFOAM REQUIRED COMPONENTS OpenFOAM finiteVolume fvOptions
      compressibleTransportModels fluidThermophysicalModels specie
      turbulenceModels compressibleTurbulenceModels dynamicFvMesh
      topoChangerFvMesh meshTools)
else()
  # > 7.x series
  find_package(OpenFOAM REQUIRED COMPONENTS OpenFOAM finiteVolume fvOptions
      compressibleTransportModels fluidThermophysicalModels specie
      turbulenceModels compressibleTurbulenceModels dynamicFvMesh
      topoChangerFvMesh meshTools)
endif()

#find_package(IMPACT COMPONENTS)
find_package(IMPACT REQUIRED)
FIND_LIBRARY(IMPACT_LIB SITCOM)
FIND_FILE(IMPACT_HDR com.h)
GET_FILENAME_COMPONENT(IMPACT_LIBPATH ${IMPACT_LIB} PATH)
GET_FILENAME_COMPONENT(IMPACT_INCPATH ${IMPACT_HDR} PATH)


# Setting library files ########################################################

# Setting include files
file(GLOB_RECURSE rocfoam_INC_FILES *.H)

set(rocfoam_INC_DIRS "")
foreach(file_path ${rocfoam_INC_FILES})
  get_filename_component(dir_path ${file_path} PATH)
  list(APPEND rocfoam_INC_DIRS ${dir_path})
endforeach()

list(APPEND rocfoam_INC_DIRS ${IMPACT_HDR})

list(REMOVE_DUPLICATES rocfoam_INC_DIRS)

# Setting cfMesh source
set(rocfoam_SRCS
    src/foam/rocFoam.C
    src/rhoCentral/rocRhoCentral.C
    src/rhoPimple/rocRhoPimple.C)
    
# Setting rhoCentral driver source
set(rocCntDrv_SRCS src/rhoCentral/rocRhoCentral_drv.C)

# Setting rhoPimple driver source
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
        ${IMPACT_LIB}
)

target_include_directories(rocfoam SYSTEM
    PRIVATE
        ${OPNF_INCLUDE_DIRS}
        ${IMPACT_HDR}
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

# Building utilities #############################################################
add_executable(rocRhoCentral ${rocCntDrv_SRCS})
target_link_libraries(rocRhoCentral rocfoam IMPACT::SITCOM)
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




#2 comfoam library ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


# Set project metadata #########################################################
set(comfoam_MAJOR_VERSION 1)
set(comfoam_MINOR_VERSION 0)
set(comfoam_PATCH_VERSION 0)
set(comfoam_EXTRA_VERSION "")
set(comfoam_VERSION "${comfoam_MAJOR_VERSION}.${comfoam_MINOR_VERSION}")
set(comfoam_VERSION
    "${comfoam_VERSION}.${comfoam_PATCH_VERSION}${comfoam_EXTRA_VERSION}")


# Finding libraries ############################################################




# Setting library files ########################################################
# Setting include files
file(GLOB_RECURSE comfoam_INC_FILES *.H)

set(comfoam_INC_DIRS "")
foreach(file_path ${comfoam_INC_FILES})
  get_filename_component(dir_path ${file_path} PATH)
  list(APPEND comfoam_INC_DIRS ${dir_path})
endforeach()

list(APPEND comfoam_INC_DIRS ${IMPACT_HDR})

list(REMOVE_DUPLICATES comfoam_INC_DIRS)

# Setting cfMesh source
set(comfoam_SRCS
    src/com/comFoam.C
    src/foam/rocFoam.C
    src/rhoCentral/rocRhoCentral.C
    src/rhoPimple/rocRhoPimple.C)

# Setting comFoam driver source
set(comFoamDrv_SRCS
    src/com/comFoam_drv.C)

# Building library #############################################################
# Build comfoam library
add_library(comfoam ${comfoam_SRCS})
target_include_directories(comfoam
    PUBLIC
        $<INSTALL_INTERFACE:include/comFoam>
)

foreach(include_dir ${comfoam_INC_DIRS})
  target_include_directories(comfoam
      PUBLIC
        $<BUILD_INTERFACE:${include_dir}>
)
endforeach()

target_compile_options(comfoam
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
target_link_libraries(comfoam
    PRIVATE
        ${OPNF_LIBRARIES}
        ${IMPACT_LIB}
)

target_include_directories(comfoam SYSTEM
    PRIVATE
        ${OPNF_INCLUDE_DIRS}
        ${IMPACT_HDR}
)

target_compile_definitions(comfoam
    PUBLIC
        ${OPNF_COMPILE_DEFINITIONS}
)

# Set soname version ###########################################################

set_target_properties(comfoam PROPERTIES VERSION ${comfoam_VERSION}
                                        SOVERSION ${comfoam_MAJOR_VERSION})

# Install libraries ############################################################

# Headers
install(FILES ${comfoam_INC_FILES}
    DESTINATION include)

# Libraries
install(TARGETS comfoam
    EXPORT comfoam
    RUNTIME DESTINATION bin
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

# Building utilities #############################################################
add_executable(comFoam ${comFoamDrv_SRCS})
target_link_libraries(comFoam comfoam IMPACT::SITCOM IMPACT::SITCOMF) 
target_include_directories(comFoam SYSTEM
    PRIVATE
        ${OPNF_INCLUDE_DIRS}
)

install(TARGETS comFoam
         RUNTIME DESTINATION bin
         PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
         GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

