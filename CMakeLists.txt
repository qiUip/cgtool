cmake_minimum_required(VERSION 3.12)

project(CGTOOL)

# Find git to create version info
find_package(Git)
# Get latest Git tag
execute_process(
  COMMAND ${GIT_EXECUTABLE} describe --tags --abbrev=0
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_LATEST_TAG
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Parse the version components from the tag
string(REGEX MATCH "v([0-9]+)\\.([0-9]+)(\\.([0-9]+))?((pre)?([0-9]+))?" VERSION_MATCHES "${GIT_LATEST_TAG}")

# Set version
if (CMAKE_MATCH_1)
  set(VERSION_MAJOR ${CMAKE_MATCH_1})
else()
  set(VERSION_MAJOR 0)
endif()
if (CMAKE_MATCH_2)
  set(VERSION_MINOR ${CMAKE_MATCH_2})
else()
  set(VERSION_MINOR 0)
endif()
if (CMAKE_MATCH_4)
  set(VERSION_PATCH ${CMAKE_MATCH_4})
else()
  set(VERSION_PATCH 0)
endif()

MESSAGE("CMAKE ${CMAKE_VERSION}")
MESSAGE("CXX   ${CMAKE_SYSTEM_NAME} ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
MESSAGE("C     ${CMAKE_SYSTEM_NAME} ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -ffast-math -march=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -march=native")


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: Release Debug ." FORCE)
endif()
# Toggleable build type
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Release;Debug")

if(NOT CMAKE_CXX_FLAGS_Release)
    set(CMAKE_CXX_FLAGS_Release "-O3")
endif()
if(NOT CMAKE_CXX_FLAGS_Debug)
    set(CMAKE_CXX_FLAGS_Debug "-g -O0")
endif()

if(${CMAKE_BUILD_TYPE} MATCHES "Release")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_Release}")
elseif(CMAKE_BUILD_TYPE MATCHES "Debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_Debug}")
endif()

# Should we ignore warnings - default ignore a few - only on gcc
option(IGNORE_SOME_WARNINGS "Don't ignore warnings" ON)
set(IGNORED_WARNINGS " -Wno-sign-compare -Wno-unused-variable -Wno-comment -Wno-unused-but-set-variable")
set(IGNORED_WARNINGS_ALWAYS " -Wno-sign-compare -Wno-reorder")

# Platform specifics
if("${CMAKE_SYSTEM_NAME}" STREQUAL "Darwin")
    # Apple OSX - don't use OPENMP
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

elseif("${CMAKE_SYSTEM_NAME}" STREQUAL "Linux")
    # GCC and Clang use different OPENMP flags
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp -pthread")

    elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        if(IGNORE_SOME_WARNINGS)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IGNORED_WARNINGS}")
            MESSAGE("IGNORED")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IGNORED_WARNINGS_ALWAYS}")
            MESSAGE("ALWAYS")
        endif()

    elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp")

    else()
        MESSAGE("Unrecognised compiler - will attempt to compile")
    endif()
endif()

message("Using " ${CMAKE_CXX_FLAGS})

# Find libxdrfile - compile it if not present
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/core)
find_library(XDRFILE_LIBRARY NAMES xdrfile
             PATHS /usr/lib /usr/local/lib ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty)
if(XDRFILE_LIBRARY_FOUND)
    MESSAGE("Found libxdrfile: ${XDRFILE_LIBRARY}")
    set (EXTRA_LIBS ${EXTRA_LIBS} ${XDRFILE_LIBRARY})
else()
    MESSAGE("Cannot find libxdrfile, including in compile")
    set(THIRDPARTY_OUTPUT_DIRCTORY ThirdParty)
    add_library(xdrfile STATIC ThirdParty/xdrfile/xdrfile.c ThirdParty/xdrfile/xdrfile_xtc.c)
    set (EXTRA_LIBS ${EXTRA_LIBS} xdrfile)
    set_target_properties(xdrfile PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
endif()

# Find Boost library - required
set(Boost_USE_STATIC_LIBS ON)
find_package(Boost REQUIRED)
find_package(Boost COMPONENTS program_options)
include_directories(${Boost_INCLUDE_DIRS})
set(EXTRA_LIBS ${EXTRA_LIBS} ${Boost_LIBRARIES})

if(Boost_PROGRAM_OPTIONS_FOUND)
    set(CMD_SRC "src/core/cmd.cpp")
else()
    set(CMD_SRC "src/core/cmd_simple.cpp")
    add_definitions(-DCMD_SIMPLE)
endif()

# Add target for documentation if we have Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc
        COMMENT "Generating API documentation with Doxygen" VERBATIM)
endif(DOXYGEN_FOUND)

# Get Git revision
execute_process(
  COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get Git remote and branch
execute_process(
  COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref --symbolic-full-name @{u}
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_REMOTE_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get Git remote URL and branch
execute_process(
  COMMAND ${GIT_EXECUTABLE} remote get-url origin
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_REMOTE_URL
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get current date and time
string(TIMESTAMP BUILD_DATETIME "%Y-%m-%d %H:%M:%S")

# Get C++ compiler name and version
set(FULL_COMPILER_VERSION "${CMAKE_SYSTEM_NAME} ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
# # Combine compiler name and version
# set(FULL_COMPILER_VERSION "${COMPILER_NAME} ${COMPILER_VERSION}")

# Configure a header file to pass the version number to the source code
configure_file(
  "${CMAKE_SOURCE_DIR}/src/version.h.in"
  "${CMAKE_SOURCE_DIR}/src/version.h"
)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Find source files
set(CGTOOLCORE_FILES
    "src/core/common.cpp"
    "src/core/frame.cpp"
    "src/core/cg_map.cpp"
    "src/core/parser.cpp"
    "src/core/residue.cpp"
    "src/core/small_functions.cpp"
    "src/core/gro_input.cpp"
    "src/core/gro_output.cpp"
    "src/core/xtc_input.cpp"
    "src/core/xtc_output.cpp"
    ${CMD_SRC})

set(CGTOOL_FILES
    "src/cgtool/cgtool.cpp"
    "src/cgtool/boltzmann_inverter.cpp"
    "src/cgtool/bond_struct.cpp"
    "src/cgtool/bondset.cpp"
    "src/cgtool/histogram.cpp"
    "src/cgtool/itp_writer.cpp"
    "src/cgtool/rdf.cpp"
    "src/cgtool/lammps_data_output.cpp"
    "src/cgtool/lammps_trj_output.cpp")

set(RAMSI_FILES
    "src/ramsi/ramsi.cpp"
    "src/ramsi/membrane.cpp")

# Add CGTOOLCORE library
add_library(cgtoolcore ${CGTOOLCORE_FILES})
target_link_libraries(cgtoolcore ${EXTRA_LIBS})

# Conditionally set the output directory for the library based on BUILD_SHARED_LIBS
if(BUILD_SHARED_LIBS)
    set_target_properties(cgtoolcore PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
else()
    set_target_properties(cgtoolcore PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
endif()

# Add CGTOOL target
add_executable(cgtool ${CGTOOL_FILES})
target_link_libraries(cgtool cgtoolcore)

# Add RAMSi target
add_executable(ramsi ${RAMSI_FILES})
target_link_libraries(ramsi cgtoolcore)

# Add xtc-length target
add_executable(xtc-length "src/core/xtc-length.c")

# Add install option
install(TARGETS xtc-length RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
install(TARGETS ramsi RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
install(TARGETS cgtool RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

# Testing
# Set the output directory for tests
set(TEST_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Testing)

# Check if the tests directory exists, and create it if not
if(NOT EXISTS ${TEST_OUTPUT_DIRECTORY})
    file(MAKE_DIRECTORY ${TEST_OUTPUT_DIRECTORY})
    message(STATUS "Created tests directory: ${TEST_OUTPUT_DIRECTORY}")
endif()# Set the output directory for test executables
set(TEST_BINARY_DIRECTORY ${TEST_OUTPUT_DIRECTORY}/bin)
#
add_subdirectory(ThirdParty/gtest-1.7.0 EXCLUDE_FROM_ALL)
include_directories(ThirdParty/gtest-1.7.0/include ThirdParty/gtest-1.7.0/include)

# Test parser
add_executable(gtest_parser EXCLUDE_FROM_ALL src/tests/parser_test.cpp)
target_link_libraries(gtest_parser gtest gtest_main cgtoolcore)
add_test(GTestParserAll ${TEST_BINARY_DIRECTORY}/gtest_parser)
# Test bondset
add_executable(gtest_bondset EXCLUDE_FROM_ALL src/tests/bondset_test.cpp
    src/cgtool/bondset.cpp src/cgtool/boltzmann_inverter.cpp src/cgtool/histogram.cpp src/cgtool/bond_struct.cpp)
target_link_libraries(gtest_bondset gtest gtest_main cgtoolcore)
add_test(GTestBondSetAll ${TEST_BINARY_DIRECTORY}/gtest_bondset)
# Test lightarray
add_executable(gtest_light_array EXCLUDE_FROM_ALL src/tests/light_array_test.cpp)
target_link_libraries(gtest_light_array gtest gtest_main)
add_test(GTestLightArrayAll ${TEST_BINARY_DIRECTORY}/gtest_light_array)
# Test small_functions
add_executable(gtest_small_functions EXCLUDE_FROM_ALL src/tests/small_functions_test.cpp)
target_link_libraries(gtest_small_functions gtest gtest_main)
add_test(GTestSmallFunctionsAll ${TEST_BINARY_DIRECTORY}/gtest_small_functions)
set_target_properties(gtest_parser gtest_bondset gtest_light_array gtest_small_functions PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TEST_BINARY_DIRECTORY})

# Integration test - does it run
add_test(NAME IntegrationRUNCGTOOL
  COMMAND ${CMAKE_COMMAND} -E env
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/cgtool -c ${CMAKE_CURRENT_SOURCE_DIR}/test_data/ALLA/cg.cfg -x ${CMAKE_CURRENT_SOURCE_DIR}/test_data/ALLA/md.xtc -g ${CMAKE_CURRENT_SOURCE_DIR}/test_data/ALLA/md.gro -i ${CMAKE_CURRENT_SOURCE_DIR}/test_data/ALLA/topol.top
  WORKING_DIRECTORY ${TEST_OUTPUT_DIRECTORY}/Temporary)
add_test(NAME IntegrationRUNRAMSi
  COMMAND ${CMAKE_COMMAND} -E env
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ramsi -c ${CMAKE_CURRENT_SOURCE_DIR}/test_data/staph/mem.cfg -x ${CMAKE_CURRENT_SOURCE_DIR}/test_data/staph/md.xtc -g ${CMAKE_CURRENT_SOURCE_DIR}/test_data/staph/md.gro
  WORKING_DIRECTORY ${TEST_OUTPUT_DIRECTORY}/Temporary)
# Integration test - check output files
add_test(IntegrationITP ${CMAKE_COMMAND} -E compare_files
  ${TEST_OUTPUT_DIRECTORY}/Temporary/ALLA.itp ${CMAKE_SOURCE_DIR}/test_data/ALLA/ALLACG.itp)
add_test(IntegrationGRO ${CMAKE_COMMAND} -E compare_files
  ${TEST_OUTPUT_DIRECTORY}/Temporary/ALLA.gro ${CMAKE_SOURCE_DIR}/test_data/ALLA/ALLACG.gro)
add_test(IntegrationTOP ${CMAKE_COMMAND} -E compare_files
  ${TEST_OUTPUT_DIRECTORY}/Temporary/ALLA.top ${CMAKE_SOURCE_DIR}/test_data/ALLA/ALLACG.top)

# Custom target for cleaning up temporary files after all tests
add_custom_target(clean_temporary_files ALL
  COMMAND ${CMAKE_COMMAND} -E remove_directory ${TEST_OUTPUT_DIRECTORY}/Temporary
  COMMENT "Cleaning up temporary files after tests"
  VERBATIM
)
# Set POST_BUILD option for the cleanup target
set_target_properties(clean_temporary_files PROPERTIES POST_BUILD COMMAND ${CMAKE_COMMAND} -E remove_directory ${TEST_OUTPUT_DIRECTORY}/Temporary)

enable_testing()
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}
                  DEPENDS gtest_parser gtest_bondset gtest_light_array gtest_small_functions cgtool ramsi clean_temporary_files)
add_custom_target(check-v COMMAND ${CMAKE_CTEST_COMMAND} "-V"
                  DEPENDS gtest_parser gtest_bondset gtest_light_array gtest_small_functions cgtool ramsi clean_temporary_files)
