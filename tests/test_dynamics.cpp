// Tests of the equations of motion and of the numerical integrator
// (pnd::DoublePendulum::derive and pnd::DoublePendulum::evolution, which is a
// classical fourth-order Runge-Kutta scheme).

#include "test_utils.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace tu = testutil;

TEST_SUITE_BEGIN("dynamics");

// ---------------------------------------------------------------------------
// State: the building block of the Runge-Kutta stages
// ---------------------------------------------------------------------------
TEST_CASE("State arithmetic is componentwise") {
  // RK4 is a linear combination of States, so a wrong operator would corrupt
  // the integrator. The values are exactly representable in binary floating
  // point and each operation is a single IEEE operation, so the comparison is
  // exact.
  pnd::State a{1., 2., 3., 4.};
  pnd::State b{0.5, -2., 8., 0.25};

  pnd::State const sum = a + b;
  CHECK(sum.theta1 == 1.5);
  CHECK(sum.theta2 == 0.);
  CHECK(sum.omega1 == 11.);
  CHECK(sum.omega2 == 4.25);

  pnd::State const scaled = a * 0.5;
  CHECK(scaled.theta1 == 0.5);
  CHECK(scaled.theta2 == 1.);
  CHECK(scaled.omega1 == 1.5);
  CHECK(scaled.omega2 == 2.);
}

// ---------------------------------------------------------------------------
// derive(): the vector field  d/dt (theta1, theta2, omega1, omega2)
// ---------------------------------------------------------------------------
TEST_CASE("derive: kinematic equations and absence of side effects") {
  // d(theta_i)/dt = omega_i is an identity of the model: the first two
  // components of the vector field are copies of the velocities (exact).
  pnd::State const s{0.3, -0.7, 1.2, -0.4};
  pnd::State initial{0.1, 0.2, 0.3, 0.4};
  auto pend = tu::make_pendulum(tu::parameter_sets[1], initial);

  pnd::State const f = pend.derive(s);
  CHECK(f.theta1 == s.omega1);
  CHECK(f.theta2 == s.omega2);

  // derive() evaluates the vector field at an arbitrary state: it must not
  // modify the state of the pendulum (RK4 relies on this for its stages).
  CHECK(pend.state().theta1 == 0.1);
  CHECK(pend.state().theta2 == 0.2);
  CHECK(pend.state().omega1 == 0.3);
  CHECK(pend.state().omega2 == 0.4);
}

TEST_CASE("derive reproduces the Lagrange equations of motion") {
  // Reference: the same equations obtained from the Lagrangian and solved as a
  // 2x2 linear system (tu::reference_accelerations). Both are exact algebra,
  // so they may only differ by floating-point rounding, amplified by the
  // conditioning of the mass matrix (measured on this grid: ~1e-13 relative).
  // A wrong sign or coefficient would give O(1) differences.
  constexpr double relative_tolerance = 1e-11;

  for (auto const &params : tu::parameter_sets) {
    pnd::State dummy{0., 0., 0., 0.};
    auto pend = tu::make_pendulum(params, dummy);

    tu::for_each_grid_state([&](pnd::State const &s) {
      auto const [alpha1, alpha2] = tu::reference_accelerations(params, s);
      pnd::State const f = pend.derive(s);

      INFO("m1=" << params.m1 << " m2=" << params.m2 << " L1=" << params.L1
                 << " L2=" << params.L2 << " state=(" << s.theta1 << ", "
                 << s.theta2 << ", " << s.omega1 << ", " << s.omega2 << ")");
      // Relative error with a floor of 1 rad/s^2, to stay meaningful where the
      // acceleration happens to be close to zero.
      CHECK(std::abs(f.omega1 - alpha1) / (1. + std::abs(alpha1)) <
            relative_tolerance);
      CHECK(std::abs(f.omega2 - alpha2) / (1. + std::abs(alpha2)) <
            relative_tolerance);
    });
  }
}

// ---------------------------------------------------------------------------
// evolution(): order of the integrator
// ---------------------------------------------------------------------------
TEST_CASE("RK4 converges to the exact solution with fourth order") {
  // The global error of a p-th order method is C dt^p, so halving dt must
  // reduce the error by 2^p. This distinguishes the expected truncation error
  // of RK4 (p = 4) from an implementation error: Euler (p = 1), midpoint or
  // Heun (p = 2) or a wrongly weighted RK4 would give a different order.
  //
  // Reference: independent high-accuracy solution (README example, 45/60 deg,
  // main.cpp parameters) at t = 1 s, accurate to ~1e-12 (see
  // tests/reference/generate_reference.py).
  pnd::State const expected{-0.81030120523121407, -0.98635741329736215,
                            0.72103875145363916, -0.40155544466652598};
  constexpr double t_final = 1.0;
  constexpr std::array<double, 3> dts{1e-2, 5e-3, 2.5e-3};

  std::array<double, 3> errors{};
  for (std::size_t i = 0; i < dts.size(); ++i) {
    pnd::State initial{0.7853981633974483, 1.0471975511965976, 0., 0.};
    auto pend = tu::make_pendulum(tu::parameter_sets[2], initial);
    int const n = tu::steps_for(t_final, dts[i]);
    for (int k = 0; k < n; ++k) {
      pend.evolution(dts[i]);
    }
    errors[i] = tu::state_distance(pend.state(), expected);
  }

  // The finest error must be far above the accuracy of the reference
  // (~1e-12), otherwise the measured order would be polluted by it.
  REQUIRE(errors[2] > 1e-10);

  // Observed order p = log2(e(dt) / e(dt/2)). For RK4 it is 4 up to O(dt)
  // corrections: measured 3.99 and 4.00 here; the window [3.8, 4.2] still
  // excludes any method of a different order.
  double const order_coarse = std::log2(errors[0] / errors[1]);
  double const order_fine = std::log2(errors[1] / errors[2]);
  CAPTURE(errors[0]);
  CAPTURE(errors[1]);
  CAPTURE(errors[2]);
  CHECK(order_coarse > 3.8);
  CHECK(order_coarse < 4.2);
  CHECK(order_fine > 3.8);
  CHECK(order_fine < 4.2);
}

// ---------------------------------------------------------------------------
// Regression tests against independent high-accuracy solutions
// ---------------------------------------------------------------------------
namespace {

struct Checkpoint {
  double time;
  pnd::State expected;
};

struct ReferenceCase {
  char const *name;
  tu::Params params;
  pnd::State initial;
  std::vector<Checkpoint> checkpoints;
  // Max-norm tolerance on (theta1, theta2, omega1, omega2) [rad, rad/s].
  double tolerance;
};

// Expected states come from tests/reference/generate_reference.py: equations
// of motion derived symbolically from the Lagrangian (sympy), integrated with
// DOP853 at rtol = atol = 1e-13 and cross-checked with the implicit Radau
// method (the two agree to better than 1.4e-12). They are NOT produced by the
// code under test. Angles are unwrapped there; the comparison is modulo 2 pi.
//
// Tolerances: the discrepancy at the default dt = 1e-3 s is the truncation
// error of RK4, which was measured (max over the checkpoints) as
//     main_params_45_60            4.7e-10   -> tolerance 5e-9   (~10x)
//     asymmetric_large_amplitude   6.5e-10   -> tolerance 1e-8   (~15x)
//     rotating_wraps_angles        4.1e-8    -> tolerance 5e-7   (~12x)
// The error is larger for the last case because the truncation error grows
// like (omega dt)^4 and omega reaches ~10 rad/s there (E = 178 J).
// Any error in the equations (sign, coefficient, swapped L1/L2 or m1/m2)
// changes the trajectory by orders of magnitude more than these tolerances.
std::vector<ReferenceCase> reference_cases() {
  return {
      {"main_params_45_60",
       {10.0, 10.0, 0.5, 0.5},
       {0.7853981633974483, 1.0471975511965976, 0., 0.},
       {{0.5,
         {0.0017409703121903236, -0.099135963489334566, -2.1361322966078569,
          -4.0232341291465676}},
        {1.0,
         {-0.81030120523121407, -0.98635741329736215, 0.72103875145363916,
          -0.40155544466652598}},
        {2.0,
         {0.84884393319858942, 0.87897714794292081, -0.939587589320238,
          -0.0089799087111530982}},
        {3.0,
         {-0.80286049061813614, -0.90413162560055849, 0.87709775885337549,
          0.82743896090646585}}},
       5e-9},

      {"asymmetric_large_amplitude",
       {2.0, 0.5, 1.3, 0.7},
       {2.0, -1.0, 0.5, -0.3},
       {{0.5,
         {1.3739029843709019, -1.1632794236103667, -3.0051794411676109,
          0.45635074960826871}},
        {1.0,
         {-0.74952270742334171, 1.3427987856060182, -4.0743756597128646,
          1.0163779948925225}},
        {2.0,
         {-1.2137195232465132, -1.8670096282212885, 3.4007848573823023,
          -2.4427106955775217}}},
       1e-8},

      // Angles go beyond +-pi several times (up to ~15 rad ~ 2.4 turns):
      // exercises the wrap-around in evolution().
      {"rotating_wraps_angles",
       {1.0, 1.0, 1.0, 1.0},
       {0., 0., 10., 6.},
       {{0.5,
         {3.875206316491671, 3.606138237190434, 6.9396255617883469,
          7.4129659142231272}},
        {1.0,
         {7.6579590106421174, 7.8980522025732558, 7.3129470684490041,
          8.6262552627156008}},
        {2.0,
         {15.40584168465389, 15.263083209953972, 6.4547083849104157,
          7.7339999774505657}}},
       5e-7},
  };
}

} // namespace

TEST_CASE("trajectories agree with independent high-accuracy references") {
  for (auto const &c : reference_cases()) {
    pnd::State initial = c.initial;
    auto pend = tu::make_pendulum(c.params, initial);

    int done = 0; // steps already performed
    for (auto const &checkpoint : c.checkpoints) {
      int const target = tu::steps_for(checkpoint.time, pnd::constants::dt);
      for (; done < target; ++done) {
        pend.evolution();
      }

      double const error =
          tu::state_distance(pend.state(), checkpoint.expected);
      INFO(c.name << " at t = " << checkpoint.time << " s");
      CHECK_LT(error, c.tolerance);
    }
  }
}

// ---------------------------------------------------------------------------
// Angle wrapping
// ---------------------------------------------------------------------------
TEST_CASE("angles are wrapped into [-pi, pi] without changing the physics") {
  SUBCASE("range of the angles during a fast rotation") {
    // Both rods spin through several full turns: the angles must stay in
    // [-pi, pi] at every step.
    //
    // Staying in range is not enough: shifting an angle by pi instead of 2 pi
    // would also land in range, but would change the physical state. A correct
    // wrap leaves sin(theta) and cos(theta) untouched, so between consecutive
    // steps they can change by at most |Delta theta| ~ |omega| dt (|sin'| and
    // |cos'| are <= 1). The bound uses the largest |omega| of the run, with a
    // 10% margin for the variation of omega inside a step.
    pnd::State initial{0., 0., 10., 6.};
    auto pend = tu::make_pendulum(tu::parameter_sets[0], initial);

    bool in_range = true;
    bool wrapped_at_least_once = false;
    double max_trig_change = 0.;
    double max_omega = 0.;
    for (int i = 0; i < 2000; ++i) {
      pnd::State const before = pend.state();
      pend.evolution();
      pnd::State const after = pend.state();

      in_range = in_range && std::abs(after.theta1) <= tu::pi &&
                 std::abs(after.theta2) <= tu::pi;
      // A jump of about 2 pi between consecutive steps can only be the wrap.
      wrapped_at_least_once = wrapped_at_least_once ||
                              std::abs(after.theta1 - before.theta1) > tu::pi ||
                              std::abs(after.theta2 - before.theta2) > tu::pi;

      max_trig_change = std::max(
          {max_trig_change,
           std::abs(std::sin(after.theta1) - std::sin(before.theta1)),
           std::abs(std::cos(after.theta1) - std::cos(before.theta1)),
           std::abs(std::sin(after.theta2) - std::sin(before.theta2)),
           std::abs(std::cos(after.theta2) - std::cos(before.theta2))});
      max_omega =
          std::max({max_omega, std::abs(before.omega1), std::abs(before.omega2),
                    std::abs(after.omega1), std::abs(after.omega2)});
    }
    CHECK(in_range);
    CHECK(wrapped_at_least_once); // the scenario really exercises the wrap
    CAPTURE(max_omega);
    CHECK_LT(max_trig_change, 1.1 * max_omega * pnd::constants::dt);
  }

  SUBCASE("theta and theta + 2 pi describe the same physical state") {
    // The equations depend on the angles only through sin and cos, hence two
    // initial conditions differing by 2 pi must evolve into the same physical
    // state. The first one starts close to pi and crosses it, so the wrap is
    // triggered. Rounding of 3.0 - 2 pi (~4e-16) is amplified only mildly in
    // this short, regular run: 1e-12 leaves a wide margin.
    auto const params = tu::parameter_sets[1];
    pnd::State a{3.0, 0.5, 2.0, 0.};
    pnd::State b{3.0 - 2. * tu::pi, 0.5, 2.0, 0.};
    auto pend_a = tu::make_pendulum(params, a);
    auto pend_b = tu::make_pendulum(params, b);

    for (int i = 0; i < 500; ++i) {
      pend_a.evolution();
      pend_b.evolution();
    }
    CHECK_LT(tu::state_distance(pend_a.state(), pend_b.state()), 1e-12);
  }
}

// ---------------------------------------------------------------------------
// Input validation (the simulation must fail loudly, not silently diverge)
// ---------------------------------------------------------------------------
TEST_CASE("non-physical input is rejected") {
  double const nan = std::numeric_limits<double>::quiet_NaN();
  double const inf = std::numeric_limits<double>::infinity();
  pnd::State const ok{0.1, 0.2, 0., 0.};

  auto construct = [](double L1, double m1, double L2, double m2,
                      pnd::State state) {
    pnd::DoublePendulum pend{L1, m1, L2, m2, state};
    (void)pend;
  };

  CHECK_NOTHROW(construct(1., 1., 1., 1., ok));

  // masses: strictly positive and finite
  CHECK_THROWS_AS(construct(1., 0., 1., 1., ok), std::invalid_argument);
  CHECK_THROWS_AS(construct(1., 1., 1., -1., ok), std::invalid_argument);
  CHECK_THROWS_AS(construct(1., nan, 1., 1., ok), std::invalid_argument);
  CHECK_THROWS_AS(construct(1., 1., 1., inf, ok), std::invalid_argument);

  // lengths: strictly positive and finite
  CHECK_THROWS_AS(construct(0., 1., 1., 1., ok), std::invalid_argument);
  CHECK_THROWS_AS(construct(1., 1., -2., 1., ok), std::invalid_argument);
  CHECK_THROWS_AS(construct(nan, 1., 1., 1., ok), std::invalid_argument);
  CHECK_THROWS_AS(construct(1., 1., inf, 1., ok), std::invalid_argument);

  // state: all four components finite
  CHECK_THROWS_AS(construct(1., 1., 1., 1., {nan, 0., 0., 0.}),
                  std::invalid_argument);
  CHECK_THROWS_AS(construct(1., 1., 1., 1., {0., inf, 0., 0.}),
                  std::invalid_argument);
  CHECK_THROWS_AS(construct(1., 1., 1., 1., {0., 0., nan, 0.}),
                  std::invalid_argument);
  CHECK_THROWS_AS(construct(1., 1., 1., 1., {0., 0., 0., -inf}),
                  std::invalid_argument);

  // A velocity so large that omega^2 overflows makes the integration diverge:
  // evolution() has to detect it instead of propagating inf/NaN.
  pnd::State huge{0., 0., 1e200, 0.};
  auto pend = tu::make_pendulum(tu::parameter_sets[0], huge);
  CHECK_THROWS_AS(pend.evolution(), std::invalid_argument);
}

TEST_SUITE_END();
