cmake_minimum_required(VERSION 3.10)

message(STATUS "ACCESSING TEST DIRECTORY: ${CMAKE_CURRENT_BINARY_DIR}")

# Add googletest subdir with its own cmake file
add_subdirectory(lib/googletest)

# Set location of test scripts, test data, and run directory
file(COPY data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
set(TEST_DATA ${CMAKE_CURRENT_BINARY_DIR}/data)
set(TEST_RESULTS ${TEST_DATA}/testResults)
file(MAKE_DIRECTORY ${TEST_RESULTS})
set(TEST_SOURCE_RHOCENTRAL ${TEST_DATA}/forwardStep/rhoCentral/run)
set(TEST_SOURCE_RHOPIMPLE ${TEST_DATA}/forwardStep/rhoPimple/run)
#set(TEST_SOURCE_ROCFOAM_RHOCENTRAL ${TEST_DATA}/forwardStep/rocFoam/rhoCentral/run)
#set(TEST_SOURCE_ROCFOAM_RHOPIMPLE ${TEST_DATA}/forwardStep/rocFoam/rhoPimple/run)

# Set environment variables for testing.
set(TEST_ENV_PATH_OPTIONS "PATH=${CMAKE_CURRENT_BINARY_DIR}:\
${CMAKE_CURRENT_BINARY_DIR}/..:\
$ENV{PATH}")
set(TEST_ENV_LD_OPTIONS "LD_LIBRARY_PATH=${CMAKE_CURRENT_BINARY_DIR}/..:\
$ENV{LD_LIBRARY_PATH}")

message(STATUS "Test Env Options: ${TEST_ENV_PATH_OPTIONS} ${TEST_ENV_LD_OPTIONS}")

#include gtest library. gtest_SOURCE_DIR is available from subdir addition
include_directories(${gtest_SOURCE_DIR}/include
    ${gtest_SOURCE_DIR}/include/gtest
    ${gtest_SOURCE_DIR})

#--------------- COM Test Executables ---------------
add_executable(runRhoCentralTest src/rocRhoCentralTests.C)
target_link_libraries(runRhoCentralTest rocrhocentral gtest gtest_main)

add_executable(runRhoPimpleTest src/rocRhoPimpleTests.C)
target_link_libraries(runRhoPimpleTest rocrhopimple gtest gtest_main)

#add_executable(runRhoFoamTest src/rocFoamTests.C)
#target_link_libraries(runRhoFoamTest rocfoam gtest gtest_main)

add_executable(normtestCentral src/normTestCentral.C)
target_link_libraries(normtestCentral rocrhocentral gtest gtest_main)

add_executable(normtestPimple src/normTestPimple.C)
target_link_libraries(normtestPimple rocrhopimple gtest gtest_main)

#--------------- rocFoam serial tests ---------------

#^^^ Testing rocRhoCentral ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
add_test(NAME rocRhoCentralSerial
    COMMAND
    ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
    runRhoCentralTest -case ${TEST_SOURCE_RHOCENTRAL}
    WORKING_DIRECTORY ${TEST_RESULTS})

add_test(NAME rocRhoCentralParallel
    COMMAND
    ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
    mpirun -np 2 runRhoCentralTest -case ${TEST_SOURCE_RHOCENTRAL} -parallel
    WORKING_DIRECTORY ${TEST_RESULTS})
#--------------------------------------------------------------------

#^^^ Testing rocRhoPimple ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
add_test(NAME rocRhoPimpleSerial
    COMMAND
    ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
    runRhoPimpleTest -case ${TEST_SOURCE_RHOPIMPLE}
    WORKING_DIRECTORY ${TEST_RESULTS})

add_test(NAME rocRhoPimpleParallel
    COMMAND
    ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
    mpirun -np 2 runRhoPimpleTest -case ${TEST_SOURCE_RHOPIMPLE} -parallel
    WORKING_DIRECTORY ${TEST_RESULTS})
#--------------------------------------------------------------------

#^^^ Testing infinity norms ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
add_test(NAME rocRhoCentralNormTest
    COMMAND
    ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
    normtestCentral -rocRhoCentral -case ${TEST_SOURCE_RHOCENTRAL}
    WORKING_DIRECTORY ${TEST_RESULTS})

add_test(NAME rocRhoPimpleNormTest
    COMMAND
    ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
    normtestPimple -rocRhoPimple -case ${TEST_SOURCE_RHOPIMPLE}
    WORKING_DIRECTORY ${TEST_RESULTS})
#--------------------------------------------------------------------
