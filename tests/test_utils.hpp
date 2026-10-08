#ifndef TEST_UTILS_HPP
#define TEST_UTILS_HPP

// Helpers shared by the test files.
//
// Everything here is written WITHOUT reusing the formulas of the code under
// test (pnd::DoublePendulum::derive / kinetic / potential): the reference
// implementations below follow different derivations, so that a mistake in
// the production code cannot be silently reproduced by the tests.

#include "pendulum.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>

namespace testutil {

inline constexpr double pi = std::numbers::pi;
inline constexpr double g = pnd::constants::g;

// The reference values in the tests (Python/sympy, see tests/reference) and the
// tolerances calibrated on the truncation error of RK4 assume these values.
// Fail at compile time, with an explicit message, if they ever change.
static_assert(pnd::constants::g == 9.81,
              "reference solutions assume g = 9.81 m/s^2: regenerate them with "
              "tests/reference/generate_reference.py");
static_assert(pnd::constants::dt == 1e-3,
              "tolerances of the default-dt tests were calibrated for "
              "dt = 1e-3 s: recalibrate them if the default time step changes");

struct Params {
  double m1;
  double m2;
  double L1;
  double L2;
};

// Several different mass/length ratios, so that no accidental symmetry
// (m1 = m2 or L1 = L2) can hide a mistake such as swapping the two rods.
inline constexpr std::array<Params, 4> parameter_sets{{
    {1.0, 1.0, 1.0, 1.0},
    {2.0, 0.5, 1.3, 0.7},
    {10.0, 10.0, 0.5, 0.5}, // the values used by main.cpp
    {0.1, 5.0, 0.4, 2.2},
}};

// DoublePendulum's constructor takes a non-const lvalue reference.
inline pnd::DoublePendulum make_pendulum(Params const &p, pnd::State state) {
  return pnd::DoublePendulum{p.L1, p.m1, p.L2, p.m2, state};
}

// Number of steps of size dt that cover a time t.
inline int steps_for(double t, double dt) {
  return static_cast<int>(std::lround(t / dt));
}

// Signed difference of two angles, in (-pi, pi]. Angles are defined modulo
// 2 pi (the simulation wraps them), so they must never be subtracted directly.
inline double angle_difference(double a, double b) {
  return std::remainder(a - b, 2. * pi);
}

// Max-norm distance between two states (angles compared modulo 2 pi).
inline double state_distance(pnd::State const &a, pnd::State const &b) {
  return std::max({std::abs(angle_difference(a.theta1, b.theta1)),
                   std::abs(angle_difference(a.theta2, b.theta2)),
                   std::abs(a.omega1 - b.omega1),
                   std::abs(a.omega2 - b.omega2)});
}

// Angular accelerations from the Lagrange equations written in matrix form,
//
//     M(theta) alpha = b(theta, omega),      delta = theta1 - theta2
//
//     M = | (m1+m2) L1^2          m2 L1 L2 cos(delta) |
//         | m2 L1 L2 cos(delta)   m2 L2^2             |
//
//     b1 = -m2 L1 L2 omega2^2 sin(delta) - (m1+m2) g L1 sin(theta1)
//     b2 = +m2 L1 L2 omega1^2 sin(delta) - m2 g L2 sin(theta2)
//
// solved with Cramer's rule. This is algebraically equivalent to the
// closed-form expressions in pendulum.cpp, but follows a different route
// (no division by m1 + m2 sin^2(delta), no manual elimination).
inline std::array<double, 2> reference_accelerations(Params const &p,
                                                     pnd::State const &s) {
  double const delta = s.theta1 - s.theta2;
  double const c = std::cos(delta);
  double const sn = std::sin(delta);

  double const m11 = (p.m1 + p.m2) * p.L1 * p.L1;
  double const m12 = p.m2 * p.L1 * p.L2 * c;
  double const m22 = p.m2 * p.L2 * p.L2;

  double const b1 = -p.m2 * p.L1 * p.L2 * s.omega2 * s.omega2 * sn -
                    (p.m1 + p.m2) * g * p.L1 * std::sin(s.theta1);
  double const b2 = p.m2 * p.L1 * p.L2 * s.omega1 * s.omega1 * sn -
                    p.m2 * g * p.L2 * std::sin(s.theta2);

  double const det = m11 * m22 - m12 * m12;
  return {(b1 * m22 - m12 * b2) / det, (m11 * b2 - m12 * b1) / det};
}

// A small deterministic grid of generic states (no RNG: reproducible on every
// platform). It includes zero and non-zero velocities and angles on both
// sides of the vertical.
inline constexpr std::array<double, 4> grid_angles{-2.8, -1.1, 0.4, 2.0};
inline constexpr std::array<double, 3> grid_omega1{-3.0, 0.0, 0.7};
inline constexpr std::array<double, 3> grid_omega2{-2.0, 0.0, 1.9};

template <typename Callable> void for_each_grid_state(Callable &&callable) {
  for (double th1 : grid_angles) {
    for (double th2 : grid_angles) {
      for (double w1 : grid_omega1) {
        for (double w2 : grid_omega2) {
          callable(pnd::State{th1, th2, w1, w2});
        }
      }
    }
  }
}

} // namespace testutil

#endif
