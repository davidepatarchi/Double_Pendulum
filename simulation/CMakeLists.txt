cmake_minimum_required(VERSION 3.28)
project(Double_Pendulum VERSION 0.1.0)

include(CTest)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

string(APPEND CMAKE_CXX_FLAGS
      " -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion"
      " -Wshadow -Wimplicit-fallthrough -Wextra-semi -Wold-style-cast"
      " -fno-omit-frame-pointer")

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  string(APPEND CMAKE_CXX_FLAGS " -D_GLIBCXX_ASSERTIONS")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  string(APPEND CMAKE_CXX_FLAGS " -D_LIBCPP_HARDENING_MODE=_LIBCPP_HARDENING_MODE_EXTENSIVE")
endif()

string(APPEND CMAKE_CXX_FLAGS_DEBUG " -fsanitize=address,undefined")
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  string(APPEND CMAKE_CXX_FLAGS " -D_GLIBCXX_SANITIZE_STD_ALLOCATOR")
endif()
string(APPEND CMAKE_EXE_LINKER_FLAGS_DEBUG " -fsanitize=address,undefined")

# Physics engine: no graphics dependency, so that it can be unit-tested without
# SFML and without opening any window.
add_library(pendulum STATIC pendulum.cpp)
target_include_directories(pendulum PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Graphical application. Turn OFF to build only the engine and the tests on a
# machine without SFML (cmake -S . -B build -DDP_BUILD_APP=OFF).
option(DP_BUILD_APP "Build the SFML graphical application" ON)
if (DP_BUILD_APP)
  find_package(SFML 3 REQUIRED COMPONENTS Graphics Window System)

  add_executable(double_pendulum main.cpp ../rendering/render.cpp)
  target_link_libraries(double_pendulum PRIVATE pendulum SFML::Graphics SFML::Window
                                         SFML::System)
endif()

# Unit tests (run with: ctest --test-dir <build-dir> --output-on-failure)
if (BUILD_TESTING)
  add_subdirectory(../tests tests)
endif()
