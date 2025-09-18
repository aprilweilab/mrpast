cmake_minimum_required(VERSION 3.12)
project(mrpast)
include(CheckCXXCompilerFlag)

option(ENABLE_TESTS "Enable automated test execution" ON)
option(ENABLE_NATIVE "Enable the compiler to use architecture-specific functionality; speeds up numerical computations but makes the resulting binaries less portable" OFF)
option(TRACE_MATRICES "Enable verbose output showing matrix values; debuging feature" OFF)

# By default we use static libraries, which makes Python packaging easier.
option(BUILD_SHARED_LIBS "Build using shared libraries" OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

# Eigen needs C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/third-party/nlopt/)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/third-party/eigen/)

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(FATAL_ERROR "In-source builds not allowed. Please make a build directory and run CMake from there.\n")
endif()

# Native build results in a solver that is faster, but less portable.
if(${ENABLE_NATIVE})
  add_compile_options(-march=native -mtune=native)
  check_cxx_compiler_flag(-mfma COMPILER_SUPPORTS_FMA)
  if(${COMPILER_SUPPORTS_FMA})
    add_compile_options(-mfma)
  endif()
endif()

if(${TRACE_MATRICES})
  add_compile_options(-DTRACE_MATRICES=1)
endif()

set(MRPAST_CORE_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/src/objective.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/solve.cpp
  )

set(MRPAST_TEST_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/test/test_derivs.cpp
  )

include_directories(${CMAKE_CURRENT_LIST_DIR}/third-party/nlohmann/
                    ${CMAKE_CURRENT_LIST_DIR}/third-party/nlopt/
                    ${CMAKE_CURRENT_LIST_DIR}/third-party/eigen/
                    ${CMAKE_CURRENT_LIST_DIR}/src/
                    )

# Library that contains the objective function and solver. Used by all the
# executables below.
add_library(mrpast STATIC ${MRPAST_CORE_SOURCES})
target_compile_options(mrpast PRIVATE -fPIC)
target_link_libraries(mrpast nlopt)

# The main MLE solver
add_executable(mrp_solver ${CMAKE_CURRENT_LIST_DIR}/src/mrp_solver.cpp)
set_target_properties(mrp_solver PROPERTIES OUTPUT_NAME "mrp-solver")
target_link_libraries(mrp_solver nlopt mrpast)

# Tool for evaluating MLE results
add_executable(mrp_eval ${CMAKE_CURRENT_LIST_DIR}/src/mrp_eval.cpp)
set_target_properties(mrp_eval PROPERTIES OUTPUT_NAME "mrp-eval")
target_link_libraries(mrp_eval nlopt mrpast)

if(ENABLE_TESTS)
  # https://github.com/google/googletest/tree/master/googletest#incorporating-into-an-existing-cmake-project

  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/5376968f6948923e2411081fd9372e71a59d8e77.zip
  )
  # For Windows: Prevent overriding the parent project's compiler/linker settings
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)

  include_directories(${CMAKE_CURRENT_LIST_DIR}/test/third-party/autodiff/)

  # Now simply link against gtest or gtest_main as needed. Eg
  add_executable(mrpast_test ${MRPAST_TEST_SOURCES})
  target_link_libraries(mrpast_test gtest_main mrpast)
  add_test(NAME mrpast_test COMMAND mrpast_test)
endif()
