# Unit tests of the physics engine (no SFML, no window: they run headless).
#
# This directory is added from simulation/CMakeLists.txt when BUILD_TESTING is
# ON (the default of include(CTest)), after the project-wide compiler flags have
# been set, so the tests are compiled with the same strict warnings and
# sanitizers as the rest of the project.

include(FetchContent)

# doctest: a single-header testing framework. An installed copy (>= 2.4.11) is
# used if found, otherwise the pinned release is downloaded at configure time.
# SYSTEM: its warnings are not reported against the project's strict flags.
FetchContent_Declare(
  doctest
  GIT_REPOSITORY https://github.com/doctest/doctest.git
  GIT_TAG v2.4.12
  GIT_SHALLOW TRUE
  SYSTEM
  FIND_PACKAGE_ARGS 2.4.11 CONFIG)
FetchContent_MakeAvailable(doctest)

add_executable(double_pendulum_tests
  main.cpp
  test_dynamics.cpp
  test_energy.cpp
  test_physics.cpp)
target_link_libraries(double_pendulum_tests PRIVATE pendulum doctest::doctest)

# One ctest entry per suite (TEST_SUITE_BEGIN in each file). On failure
# `ctest --output-on-failure` prints the doctest report with the offending
# assertion, the values involved and the INFO/CAPTURE context.
foreach(suite IN ITEMS dynamics energy physics)
  add_test(NAME ${suite} COMMAND double_pendulum_tests --test-suite=${suite})
endforeach()
