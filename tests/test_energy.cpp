// Tests of the mechanical energy E = T + V:
//   - potential() and kinetic() against independent geometric derivations,
//   - consistency of the energy with the equations of motion (dE/dt = 0),
//   - conservation during a simulation and scaling of the drift with dt.

#include "test_utils.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>

namespace tu = testutil;

namespace {

// Relative comparison for quantities obtained with a handful of floating-point
// operations: a few ulp are expected (measured on the grids below: < 1e-15),
// 1e-13 leaves a margin of ~100 and still catches any wrong term.
constexpr double few_flops_tolerance = 1e-13;

doctest::Approx exactly(double expected) {
  return doctest::Approx(expected).epsilon(few_flops_tolerance).scale(0.);
}

double total_energy(pnd::DoublePendulum &pend) {
  return pend.kinetic() + pend.potential();
}

// Largest |E(t) - E(0)| / E(0) over a simulation of the given duration.
double max_relative_energy_drift(tu::Params const &params,
                                 pnd::State const &initial, double dt,
                                 double duration) {
  pnd::State state = initial;
  auto pend = tu::make_pendulum(params, state);
  double const e0 = total_energy(pend);

  double max_drift = 0.;
  int const steps = tu::steps_for(duration, dt);
  for (int i = 0; i < steps; ++i) {
    pend.evolution(dt);
    max_drift = std::max(max_drift, std::abs(total_energy(pend) - e0) / e0);
  }
  return max_drift;
}

} // namespace

TEST_SUITE_BEGIN("energy");

// ---------------------------------------------------------------------------
// Potential energy
// ---------------------------------------------------------------------------
TEST_CASE("potential energy is m g h with zero at the rest position") {
  // V = g * sum_i m_i h_i, where h_i is the height of bob i above its height
  // in the rest configuration (-L1 for bob 1, -(L1 + L2) for bob 2).
  auto const p = tu::parameter_sets[1];
  pnd::State state{0., 0., 0., 0.};
  auto pend = tu::make_pendulum(p, state);

  SUBCASE("rest configuration") {
    // 1 - cos(0) is exactly 0 in floating point: exact comparison is valid.
    CHECK(pend.potential() == 0.);
  }

  SUBCASE("both rods straight up: heights 2 L1 and 2 (L1 + L2)") {
    pend.state() = {tu::pi, tu::pi, 0., 0.};
    CHECK(pend.potential() ==
          exactly(tu::g * (p.m1 * 2. * p.L1 + p.m2 * 2. * (p.L1 + p.L2))));
  }

  SUBCASE("both rods horizontal: heights L1 and L1 + L2") {
    pend.state() = {tu::pi / 2., tu::pi / 2., 0., 0.};
    CHECK(pend.potential() ==
          exactly(tu::g * (p.m1 * p.L1 + p.m2 * (p.L1 + p.L2))));
  }

  SUBCASE("first rod horizontal, second hanging: both bobs at height L1") {
    pend.state() = {tu::pi / 2., 0., 0., 0.};
    CHECK(pend.potential() == exactly(tu::g * (p.m1 + p.m2) * p.L1));
  }

  SUBCASE("first rod hanging, second horizontal: only bob 2 at height L2") {
    pend.state() = {0., tu::pi / 2., 0., 0.};
    CHECK(pend.potential() == exactly(tu::g * p.m2 * p.L2));
  }
}

TEST_CASE("potential energy is consistent with the bob positions") {
  // Same quantity computed from the Cartesian positions returned by
  // pos_y1/pos_y2 (the ones used for rendering): it ties together the angle
  // convention of the geometry and the one of the energy.
  for (auto const &p : tu::parameter_sets) {
    pnd::State state{0., 0., 0., 0.};
    auto pend = tu::make_pendulum(p, state);

    tu::for_each_grid_state([&](pnd::State const &s) {
      pend.state() = s;
      double const y1 = pend.pos_y1(0.);
      double const y2 = pend.pos_y2(y1);
      double const expected =
          tu::g * (p.m1 * (y1 + p.L1) + p.m2 * (y2 + p.L1 + p.L2));
      CHECK(std::abs(pend.potential() - expected) / (1. + expected) <
            few_flops_tolerance);
    });
  }
}

// ---------------------------------------------------------------------------
// Kinetic energy
// ---------------------------------------------------------------------------
TEST_CASE("kinetic energy vanishes at zero angular velocity") {
  // Every term of T contains a factor omega: with omega1 = omega2 = 0 the
  // result is exactly zero whatever the angles are.
  for (auto const &p : tu::parameter_sets) {
    for (double th1 : tu::grid_angles) {
      for (double th2 : tu::grid_angles) {
        pnd::State state{th1, th2, 0., 0.};
        auto pend = tu::make_pendulum(p, state);
        CHECK(pend.kinetic() == 0.);
      }
    }
  }
}

TEST_CASE("kinetic energy equals the sum of the Cartesian kinetic energies") {
  // T = 1/2 m1 |v1|^2 + 1/2 m2 |v2|^2 with
  //   v1 = L1 omega1 (cos th1, sin th1)
  //   v2 = v1 + L2 omega2 (cos th2, sin th2)
  // This independent derivation checks in particular the coupling term
  // m2 L1 L2 omega1 omega2 cos(th1 - th2), and that T is never negative.
  for (auto const &p : tu::parameter_sets) {
    pnd::State state{0., 0., 0., 0.};
    auto pend = tu::make_pendulum(p, state);

    tu::for_each_grid_state([&](pnd::State const &s) {
      double const v1x = p.L1 * s.omega1 * std::cos(s.theta1);
      double const v1y = p.L1 * s.omega1 * std::sin(s.theta1);
      double const v2x = v1x + p.L2 * s.omega2 * std::cos(s.theta2);
      double const v2y = v1y + p.L2 * s.omega2 * std::sin(s.theta2);
      double const expected = 0.5 * p.m1 * (v1x * v1x + v1y * v1y) +
                              0.5 * p.m2 * (v2x * v2x + v2y * v2y);

      pend.state() = s;
      double const kinetic = pend.kinetic();
      CHECK(kinetic >= 0.);
      CHECK(std::abs(kinetic - expected) / (1. + expected) <
            few_flops_tolerance);
    });
  }
}

TEST_CASE("kinetic energy in the two single-rotation limits") {
  auto const p = tu::parameter_sets[1];
  constexpr double w = 1.7;
  pnd::State state{0.4, -1.1, 0., 0.};
  auto pend = tu::make_pendulum(p, state);

  // omega2 = 0: rod 2 translates rigidly with bob 1, so both bobs move with
  // velocity L1 omega1 whatever theta2 is: T = 1/2 (m1 + m2) L1^2 omega1^2.
  // (The first term of T in the project README lacks the m2: that expression
  // would fail this check, the code is the correct one.)
  pend.state() = {0.4, -1.1, w, 0.};
  CHECK(pend.kinetic() == exactly(0.5 * (p.m1 + p.m2) * p.L1 * p.L1 * w * w));

  // omega1 = 0: bob 1 is at rest and bob 2 moves on a circle of radius L2
  // around it: T = 1/2 m2 L2^2 omega2^2.
  pend.state() = {0.4, -1.1, 0., w};
  CHECK(pend.kinetic() == exactly(0.5 * p.m2 * p.L2 * p.L2 * w * w));
}

// ---------------------------------------------------------------------------
// Total energy: consistency with the independent reference values
// ---------------------------------------------------------------------------
TEST_CASE("total energy of the reference states matches the symbolic value") {
  // E(0) of the states used in the regression tests, computed independently
  // with sympy from the Lagrangian (tests/reference/generate_reference.py),
  // with V = 0 at the rest position. The asymmetric case has non-zero initial
  // velocities, so it checks T as well as V.
  struct Reference {
    tu::Params params;
    pnd::State state;
    double energy;
  };
  std::array<Reference, 3> const references{{
      {tu::parameter_sets[2],
       {0.7853981633974483, 1.0471975511965976, 0., 0.},
       53.257824765599693},
      {tu::parameter_sets[1], {2.0, -1.0, 0.5, -0.3}, 47.335390536908989},
      {tu::parameter_sets[0], {0., 0., 10., 6.}, 178.},
  }};

  for (auto const &r : references) {
    pnd::State state = r.state;
    auto pend = tu::make_pendulum(r.params, state);
    CHECK(total_energy(pend) == exactly(r.energy));
  }
}

// ---------------------------------------------------------------------------
// Energy and equations of motion must be mutually consistent
// ---------------------------------------------------------------------------
TEST_CASE("total energy is a first integral of the equations of motion") {
  // For the exact dynamics dE/dt = 0. Along the vector field f = derive(s),
  //   dE/dt = (E(s + h f) - E(s - h f)) / (2 h) + O(h^2).
  // It holds only if energy() and derive() describe the SAME system, so it
  // catches a sign or coefficient error in either of them even when each one
  // looks plausible alone.
  //
  // h = 1e-5: truncation error O(h^2) ~ 1e-10 and rounding error
  // ~ eps E / h ~ 1e-9 (measured: ~2e-9 on a scale of ~25). Tolerance: 1e-6 of
  // the natural energy rate (m1 + m2) g (L1 + L2) (1 + |omega1| + |omega2|):
  // a relative error of 1e-6 in any coefficient of the equations would
  // already be detected, a sign error gives O(1).
  constexpr double h = 1e-5;
  constexpr double relative_tolerance = 1e-6;

  for (auto const &p : tu::parameter_sets) {
    pnd::State dummy{0., 0., 0., 0.};
    auto pend = tu::make_pendulum(p, dummy);

    tu::for_each_grid_state([&](pnd::State const &s) {
      // State's operators are non-const member functions: use mutable copies.
      pnd::State s0 = s;
      pnd::State f = pend.derive(s0);
      pnd::State f_plus = f * h;
      pnd::State f_minus = f * (-h);

      pend.state() = s0 + f_plus;
      double const e_plus = total_energy(pend);
      pend.state() = s0 + f_minus;
      double const e_minus = total_energy(pend);

      double const dedt = (e_plus - e_minus) / (2. * h);
      double const scale = (p.m1 + p.m2) * tu::g * (p.L1 + p.L2) *
                           (1. + std::abs(s.omega1) + std::abs(s.omega2));
      INFO("state=(" << s.theta1 << ", " << s.theta2 << ", " << s.omega1 << ", "
                     << s.omega2 << ")");
      CHECK(std::abs(dedt) < relative_tolerance * scale);
    });
  }
}

// ---------------------------------------------------------------------------
// Energy conservation during the simulation
// ---------------------------------------------------------------------------
TEST_CASE("energy is conserved over 5 s at the default time step") {
  // The system is conservative, so E(t) = E(0) exactly. The integrator is RK4,
  // which is NOT symplectic: the energy error is not zero, it is the
  // truncation error ~ C dt^4 (its scaling is verified in the next test).
  // Tolerances are therefore tied to the measured RK4 error at dt = 1e-3 s
  // (max over 5 s) with a margin of ~40-70x:
  //
  //   case                         measured    tolerance
  //   main.cpp parameters, 45/60   2.8e-12     1e-10
  //   asymmetric, large amplitude  1.2e-10     5e-9
  //   chaotic, near inverted       3.5e-9      2e-7
  //   fast rotation (wraps angles) 7.0e-10     5e-8
  //
  // The error grows with the energy because the local truncation error grows
  // with the size of the derivatives (~ (omega dt)^4). A real mistake (wrong
  // equation, wrong integrator, sign error) gives relative errors of order
  // 1e-3 to 1 on these time scales.
  struct Case {
    char const *name;
    tu::Params params;
    pnd::State initial;
    double tolerance;
  };
  std::array<Case, 4> const cases{{
      {"main.cpp parameters, 45/60 deg",
       tu::parameter_sets[2],
       {0.7853981633974483, 1.0471975511965976, 0., 0.},
       1e-10},
      {"asymmetric, large amplitude",
       tu::parameter_sets[1],
       {2.0, -1.0, 0.5, -0.3},
       5e-9},
      {"chaotic, near the inverted position",
       tu::parameter_sets[0],
       {3.0, 2.5, 0., 0.},
       2e-7},
      {"fast rotation", tu::parameter_sets[0], {0., 0., 6., -4.}, 5e-8},
  }};

  for (auto const &c : cases) {
    double const drift =
        max_relative_energy_drift(c.params, c.initial, pnd::constants::dt, 5.);
    INFO(c.name);
    CHECK_LT(drift, c.tolerance);
  }
}

TEST_CASE("energy drift scales as dt^4 (RK4 truncation error)") {
  // Halving dt must reduce the error of a 4th-order method by 2^4 = 16. This
  // needs no external reference and separates the unavoidable numerical error
  // from an implementation error: measured ratios are 15.8 and 15.9 here
  // (15.3-16.3 on four different trajectories). Euler would give 2, RK2 4,
  // RK3 8. The window [12, 20] corresponds to an order between 3.6 and 4.3.
  auto const p = tu::parameter_sets[1];
  pnd::State const initial{2.0, -1.0, 0.5, -0.3};

  double const drift_coarse = max_relative_energy_drift(p, initial, 2e-3, 5.);
  double const drift_middle = max_relative_energy_drift(p, initial, 1e-3, 5.);
  double const drift_fine = max_relative_energy_drift(p, initial, 5e-4, 5.);

  // Must stay well above the rounding noise floor of E (~1e-15), otherwise the
  // ratio would measure rounding instead of truncation error.
  REQUIRE(drift_fine > 1e-12);

  double const ratio_coarse = drift_coarse / drift_middle;
  double const ratio_fine = drift_middle / drift_fine;
  CAPTURE(drift_coarse);
  CAPTURE(drift_middle);
  CAPTURE(drift_fine);
  CHECK(ratio_coarse > 12.);
  CHECK(ratio_coarse < 20.);
  CHECK(ratio_fine > 12.);
  CHECK(ratio_fine < 20.);
}

TEST_SUITE_END();
