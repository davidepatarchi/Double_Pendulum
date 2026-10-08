cmake_minimum_required(VERSION 3.28)
project(Double_Pendulum VERSION 0.1.0)

# enable_testing() must be called in the top-level directory for
# `ctest --test-dir build` to find the tests.
include(CTest)

# Put the executables directly in the build directory (build/Debug/progetto
# with Ninja Multi-Config), as documented in the README.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Project-wide settings (compiler flags, sanitizers), the physics library, the
# application and the tests live in simulation/CMakeLists.txt, which can also
# still be configured on its own: cmake -S simulation -B simulation/build
add_subdirectory(simulation)
