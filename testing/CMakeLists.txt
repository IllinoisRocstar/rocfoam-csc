cmake_minimum_required(VERSION 3.10)

message(STATUS "ACCESSING TEST DIRECTORY: ${CMAKE_CURRENT_BINARY_DIR}")

set (ENABLE_MPI TRUE CACHE BOOL "Build with MPI Support")
set (TEST_ENV_PATH_OPTIONS "PATH=${CMAKE_CURRENT_BINARY_DIR}/:${CMAKE_CURRENT_BINARY_DIR}/../:$ENV{PATH}")

set (TEST_ENV_LD_OPTIONS
     "LD_LIBRARY_PATH=${CMAKE_CURRENT_BINARY_DIR}/../:\
     $ENV{LD_LIBRARY_PATH}")

add_subdirectory(lib/googletest)

#Copies data from the source directory to the build directory
file (COPY ${CMAKE_CURRENT_SOURCE_DIR}/data
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

#add subdir in for the test results
#file (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data/testResults)
#SET (TEST_DATA ${CMAKE_CURRENT_BINARY_DIR}/data)
#SET (TEST_RESULTS ${TEST_DATA}/testResults)

#finding libraries
#MESSAGE(STATUS "Prefix path = ${CMAKE_PREFIX_PATH}")

#adding gtest subdir with its own cmake file
#find_package(GTest REQUIRED)
#FIND_LIBRARY(GTEST_LIB gtest)
#FIND_FILE(GTEST_INC gtest.h)

#find_library(GTEST_LIB gtest)

# Setting include files
file(GLOB_RECURSE INC_FILES *.H *.h)

set(TEST_INC_DIRS "")
foreach(file_path ${INC_FILES})
  get_filename_component(dir_path ${file_path} PATH)
  list(APPEND TEST_INC_DIRS ${dir_path})
endforeach()

#list(APPEND TEST_INC_DIRS ${GTEST_INC})
list(REMOVE_DUPLICATES TEST_INC_DIRS)

# Setting cfMesh source
set(rhocentraltests_SRCS
    rocRhoCentralTests.C)

#MESSAGE(STATUS "gtest_SOURCE_DIR = ${GTest_SOURCE_DIR}")
#include gtest library. gtest_SOURCE_DIR is available from subdir addition
include_directories( ${gtest_SOURCE_DIR}/include
                     ${gtest_SOURCE_DIR}/include/gtest
                     ${gtest_SOURCE_DIR})


#MESSAGE(STATUS "Test ENV_PATH_OPTIONS: ${TEST_ENV_PATH_OPTIONS}")
#MESSAGE(STATUS "Test ENV_LD_OPTIONS: ${TEST_ENV_LD_OPTIONS}")
#--------------- COM Test Executables ---------------
add_executable(runRhoCentralTest src/rocRhoCentralTests.C)
target_link_libraries(runRhoCentralTest rocfoam gtest gtest_main)
target_include_directories(runRhoCentralTest
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../>
)

add_executable(runRhoPimpleTest src/rocRhoPimpleTests.C)
target_link_libraries(runRhoPimpleTest rocfoam gtest gtest_main)
target_include_directories(runRhoPimpleTest
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../>
)

add_executable(runRhoFoamTest src/rocFoamTests.C)
target_link_libraries(runRhoFoamTest rocfoam gtest gtest_main)
target_include_directories(runRhoFoamTest
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../>
)

add_executable(normtest src/normTest.C)
target_link_libraries(normtest rocfoam gtest gtest_main)
target_include_directories(normtest
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../>
)


#--------------- rocFoam serial tests ---------------

#add subdir in for the test results
#file (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data/testResults)

set (TEST_SOURCE_RHOCENTRAL ${CMAKE_CURRENT_BINARY_DIR}/data/forwardStep/rhoCentral/run)
set (TEST_SOURCE_RHOPIMPLE  ${CMAKE_CURRENT_BINARY_DIR}/data/forwardStep/rhoPimple/run)

set (TEST_SOURCE_ROCFOM_RHOCENTRAL ${CMAKE_CURRENT_BINARY_DIR}/data/forwardStep/rocFoam/rhoCentral/run)
set (TEST_SOURCE_ROCFOM_RHOPIMPLE  ${CMAKE_CURRENT_BINARY_DIR}/data/forwardStep/rocFoam/rhoPimple/run)


#MESSAGE(STATUS "TEST_SOURCE_RHOCENTRAL: ${TEST_SOURCE_RHOCENTRAL}")
#MESSAGE(STATUS "TEST_SOURCE_RHOPimple : ${TEST_SOURCE_RHOPIMPLE}")

#SET (TEST_RESULTS ${TEST_DATA}/testResults)
#add_test(NAME Usage COMMAND runRhoCentralTest)


#^^^ Tesing rocRhoCentral ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
add_test(NAME rocRhoCentralSerial
         COMMAND
         runRhoCentralTest -case ${TEST_SOURCE_RHOCENTRAL}
         WORKING_DIRECTORY ${TEST_RESULTS})

add_test(NAME rocRhoCentralParallel
         COMMAND
         mpirun -np 2 runRhoCentralTest -case ${TEST_SOURCE_RHOCENTRAL} -parallel
         WORKING_DIRECTORY ${TEST_RESULTS})
#--------------------------------------------------------------------

#^^^ Tesing rocRhoPimple ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
add_test(NAME rocRhoPimpleSerial
         COMMAND
         ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}" 
         runRhoPimpleTest -case ${TEST_SOURCE_RHOPIMPLE}
         WORKING_DIRECTORY ${TEST_RESULTS})

add_test(NAME rocRhoPimpleParallel
         COMMAND
         mpirun -np 2 runRhoPimpleTest -case ${TEST_SOURCE_RHOPIMPLE} -parallel
         WORKING_DIRECTORY ${TEST_RESULTS})
#--------------------------------------------------------------------

#^^^ Tesing rocFoam ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
add_test(NAME rocFoamRhoCentralSerial
         COMMAND
         runRhoFoamTest -case ${TEST_SOURCE_ROCFOM_RHOCENTRAL} -rocRhoCentral
         WORKING_DIRECTORY ${TEST_RESULTS})

add_test(NAME rocFoamRhoCentralParallel
         COMMAND
         mpirun -np 2 runRhoFoamTest -case ${TEST_SOURCE_ROCFOM_RHOCENTRAL} -rocRhoCentral -parallel
         WORKING_DIRECTORY ${TEST_RESULTS})


add_test(NAME rocFoamRhoPimpleSerial
         COMMAND
         runRhoFoamTest -case ${TEST_SOURCE_ROCFOM_RHOPIMPLE} -rocRhoPimple
         WORKING_DIRECTORY ${TEST_RESULTS})

add_test(NAME rocFoamRhoPimpleParallel
         COMMAND
         mpirun -np 2 runRhoFoamTest -case ${TEST_SOURCE_ROCFOM_RHOPIMPLE} -rocRhoPimple -parallel
         WORKING_DIRECTORY ${TEST_RESULTS})
#--------------------------------------------------------------------

#^^^ Tesing infinity norms ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
add_test(NAME rocRhoCentralNormTest
         COMMAND
         normtest -rocRhoCentral -case ${TEST_SOURCE_RHOCENTRAL}
         WORKING_DIRECTORY ${TEST_RESULTS})

add_test(NAME rocRhoPimpleNormTest
         COMMAND
         normtest -rocRhoPimple -case ${TEST_SOURCE_RHOPIMPLE}
         WORKING_DIRECTORY ${TEST_RESULTS})
#--------------------------------------------------------------------
