// Tests of physical properties of the double pendulum that can be compared with
// exact analytic results: equilibria, small oscillations, the single-pendulum
// limit, symmetries and scaling laws, and the geometry of the model.

#include "test_utils.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>

namespace tu = testutil;

namespace {

// ---------------------------------------------------------------------------
// Normal modes of the linearised system
// ---------------------------------------------------------------------------
struct NormalMode {
  double omega; // angular frequency
  double u1;    // unit-norm eigenvector (amplitudes of theta1 and theta2)
  double u2;
};

// For small angles (sin x ~ x, cos x ~ 1, omega^2 terms dropped) the equations
// of motion of the double pendulum with absolute angles reduce to
//
//     (m1 + m2) L1 theta1'' + m2 L2 theta2'' + (m1 + m2) g theta1 = 0
//     L1 theta1'' + L2 theta2'' + g theta2 = 0
//
// With theta = u cos(omega t) and lambda = omega^2 this is a 2x2 homogeneous
// system, whose determinant vanishes for
//
//     m1 L1 L2 lambda^2 - (m1 + m2) g (L1 + L2) lambda + (m1 + m2) g^2 = 0.
//
// The second equation gives the eigenvector: u1 / u2 = (g - L2 lambda) /
// (L1 lambda).
std::array<NormalMode, 2> normal_modes(tu::Params const &p) {
  double const a = p.m1 * p.L1 * p.L2;
  double const b = -(p.m1 + p.m2) * tu::g * (p.L1 + p.L2);
  double const c = (p.m1 + p.m2) * tu::g * tu::g;
  double const disc = std::sqrt(b * b - 4. * a * c);

  std::array<NormalMode, 2> modes{};
  std::array<double, 2> const lambdas{(-b - disc) / (2. * a),
                                      (-b + disc) / (2. * a)};
  for (std::size_t i = 0; i < 2; ++i) {
    double const lambda = lambdas[i];
    double const u1 = (tu::g - p.L2 * lambda) / (p.L1 * lambda);
    double const norm = std::hypot(u1, 1.);
    modes[i] = {std::sqrt(lambda), u1 / norm, 1. / norm};
  }
  return modes;
}

// ---------------------------------------------------------------------------
// Exact quarter period of a simple pendulum released at rest from theta0
// ---------------------------------------------------------------------------
// T/4 = K(k) / omega0 with omega0 = sqrt(g / L), k = sin(theta0 / 2) and K the
// complete elliptic integral of the first kind, K(k) = pi / (2 AGM(1, k')),
// k' = sqrt(1 - k^2) = cos(theta0 / 2). The arithmetic-geometric mean converges
// quadratically: 30 iterations are far more than enough.
double simple_pendulum_quarter_period(double theta0, double length) {
  double a = 1.;
  double b = std::cos(theta0 / 2.);
  for (int i = 0; i < 30; ++i) {
    double const next_a = 0.5 * (a + b);
    b = std::sqrt(a * b);
    a = next_a;
  }
  double const elliptic_k = tu::pi / (2. * a);
  return elliptic_k / std::sqrt(tu::g / length);
}

} // namespace

TEST_SUITE_BEGIN("physics");

// ---------------------------------------------------------------------------
// Geometry of the model
// ---------------------------------------------------------------------------
TEST_CASE("geometry: positions follow the documented convention") {
  // Convention (README): angles from the DOWNWARD vertical, x to the right,
  // y upwards; bob 2 is positioned relative to bob 1 (pos_x2/pos_y2 take the
  // position of bob 1 as origin, as main.cpp does).
  auto const p = tu::parameter_sets[1];
  pnd::State state{0., 0., 0., 0.};
  auto pend = tu::make_pendulum(p, state);

  SUBCASE("rest configuration: both bobs straight below the pivot") {
    double const x1 = pend.pos_x1(0.);
    double const y1 = pend.pos_y1(0.);
    CHECK(x1 == doctest::Approx(0.).epsilon(1e-15));
    CHECK(y1 == doctest::Approx(-p.L1).epsilon(1e-15));
    CHECK(pend.pos_x2(x1) == doctest::Approx(0.).epsilon(1e-15));
    CHECK(pend.pos_y2(y1) == doctest::Approx(-(p.L1 + p.L2)).epsilon(1e-15));
  }

  SUBCASE("positive theta moves the bob to the right (+x)") {
    // theta1 = pi/2: first rod horizontal to the right. theta2 = -pi/2: second
    // rod horizontal to the LEFT of bob 1.
    pend.state() = {tu::pi / 2., -tu::pi / 2., 0., 0.};
    double const x1 = pend.pos_x1(0.);
    double const y1 = pend.pos_y1(0.);
    CHECK(x1 == doctest::Approx(p.L1).epsilon(1e-15));
    CHECK(std::abs(y1) < 1e-15);
    CHECK(pend.pos_x2(x1) == doctest::Approx(p.L1 - p.L2).epsilon(1e-15));
    CHECK(std::abs(pend.pos_y2(y1)) < 1e-15);
  }

  SUBCASE("the origin only translates the picture") {
    pend.state() = {0.7, -1.9, 0., 0.};
    CHECK(pend.pos_x1(2.5) == doctest::Approx(2.5 + pend.pos_x1(0.)));
    CHECK(pend.pos_y1(-4.) == doctest::Approx(-4. + pend.pos_y1(0.)));
  }
}

TEST_CASE("geometry: the rods are rigid for every configuration") {
  // |bob 1 - pivot| = L1 and |bob 2 - bob 1| = L2 follow from
  // sin^2 + cos^2 = 1 (a handful of operations: ~1e-16 expected).
  for (auto const &p : tu::parameter_sets) {
    pnd::State state{0., 0., 0., 0.};
    auto pend = tu::make_pendulum(p, state);

    for (double th1 : tu::grid_angles) {
      for (double th2 : tu::grid_angles) {
        pend.state() = {th1, th2, 0., 0.};
        double const x1 = pend.pos_x1(0.);
        double const y1 = pend.pos_y1(0.);
        double const x2 = pend.pos_x2(x1);
        double const y2 = pend.pos_y2(y1);

        CHECK(std::hypot(x1, y1) == doctest::Approx(p.L1).epsilon(1e-13));
        CHECK(std::hypot(x2 - x1, y2 - y1) ==
              doctest::Approx(p.L2).epsilon(1e-13));
      }
    }
  }
}

// ---------------------------------------------------------------------------
// Equilibria
// ---------------------------------------------------------------------------
TEST_CASE("the four equilibrium configurations are fixed points") {
  // With each rod pointing straight down or straight up and zero velocities
  // the net torque on every rod is zero (sin 0 = sin pi = 0), so the state
  // must not evolve.
  //  - (0, 0), stable: sin(0) = 0 exactly, the vector field is exactly zero
  //    and the state must remain bitwise unchanged.
  //  - The others are unstable: sin(pi) = 1.2e-16 in floating point is a tiny
  //    perturbation, which can grow at most like exp(lambda t) with lambda <
  //    ~5 /s for these parameters (about 100x in 1 s). The deviation after 1 s
  //    is therefore ~1e-14; the tolerance 1e-9 leaves a wide margin and is
  //    still ~9 orders of magnitude smaller than any real motion.
  auto const p = tu::parameter_sets[1];

  for (double th1 : {0., tu::pi}) {
    for (double th2 : {0., tu::pi}) {
      pnd::State const rest{th1, th2, 0., 0.};
      pnd::State state = rest;
      auto pend = tu::make_pendulum(p, state);

      for (int i = 0; i < 1000; ++i) {
        pend.evolution();
      }

      double const distance = tu::state_distance(pend.state(), rest);
      INFO("equilibrium (" << th1 << ", " << th2 << ")");
      if (th1 == 0. && th2 == 0.) {
        CHECK(distance == 0.);
      } else {
        CHECK_LT(distance, 1e-9);
      }
    }
  }
}

TEST_CASE("gravity is restoring near the stable equilibrium") {
  // Negative control for the test above (it must not pass because the
  // system never moves): released at rest at a small positive theta1, the first
  // rod is pulled back towards the vertical, i.e. towards negative theta1.
  auto const p = tu::parameter_sets[1];
  pnd::State state{0.1, 0., 0., 0.};
  auto pend = tu::make_pendulum(p, state);

  CHECK(pend.derive(pend.state()).omega1 < 0.);

  for (int i = 0; i < 100; ++i) { // 0.1 s
    pend.evolution();
  }
  CHECK(pend.state().theta1 < 0.1);
  CHECK(pend.state().omega1 < 0.);
}

// ---------------------------------------------------------------------------
// Small oscillations: exact normal modes
// ---------------------------------------------------------------------------
TEST_CASE("small oscillations reproduce the analytic normal modes") {
  // Released at rest along an eigenvector with a small amplitude, the system
  // must oscillate as a single harmonic: theta_i(t) = eps u_i cos(omega t).
  // This verifies the masses, the lengths and the coupling terms of the
  // equations at once (swapping L1/L2 or m1/m2 changes omega by O(1)).
  //
  // The only deviation from the harmonic solution is the leading nonlinear
  // correction, relative O(eps^2): measured max relative errors over 3 periods
  // are 3e-9 and 3e-8 for eps = 1e-4 (the RK4 error is negligible, ~1e-12).
  // Tolerance: 1e-6 of the amplitude (30x margin); a relative error of only
  // 1e-6 in omega would already give a phase error 2 pi * 3 * 1e-6 ~ 2e-5.
  constexpr double eps = 1e-4;
  constexpr double relative_tolerance = 1e-6;
  auto const p = tu::parameter_sets[1];

  auto const modes = normal_modes(p);

  // Sanity check of the test itself: closed-form frequency from the
  // literature, omega^2 = g / (2 m1 L1 L2) [ (m1 + m2)(L1 + L2) +- sqrt(...) ].
  double const disc =
      std::sqrt((p.m1 + p.m2) * (p.m1 + p.m2) * (p.L1 + p.L2) * (p.L1 + p.L2) -
                4. * p.m1 * (p.m1 + p.m2) * p.L1 * p.L2);
  double const prefactor = tu::g / (2. * p.m1 * p.L1 * p.L2);
  CHECK(modes[1].omega * modes[1].omega ==
        doctest::Approx(prefactor * ((p.m1 + p.m2) * (p.L1 + p.L2) + disc))
            .epsilon(1e-12));
  CHECK(modes[0].omega * modes[0].omega ==
        doctest::Approx(prefactor * ((p.m1 + p.m2) * (p.L1 + p.L2) - disc))
            .epsilon(1e-12));

  for (auto const &mode : modes) {
    // ... and the eigenvector must also satisfy the first linearised equation
    double const lambda = mode.omega * mode.omega;
    double const residual = (p.m1 + p.m2) * (tu::g - p.L1 * lambda) * mode.u1 -
                            p.m2 * p.L2 * lambda * mode.u2;
    CHECK(std::abs(residual) < 1e-12 * (p.m1 + p.m2) * tu::g);

    pnd::State state{eps * mode.u1, eps * mode.u2, 0., 0.};
    auto pend = tu::make_pendulum(p, state);

    double const period = 2. * tu::pi / mode.omega;
    int const steps = tu::steps_for(3. * period, pnd::constants::dt);
    double max_angle_error = 0.;
    double max_omega_error = 0.;
    for (int i = 1; i <= steps; ++i) {
      pend.evolution();
      double const t = i * pnd::constants::dt;
      double const c = std::cos(mode.omega * t);
      double const s = std::sin(mode.omega * t);
      auto const &now = pend.state();
      max_angle_error =
          std::max({max_angle_error, std::abs(now.theta1 - eps * mode.u1 * c),
                    std::abs(now.theta2 - eps * mode.u2 * c)});
      max_omega_error =
          std::max({max_omega_error,
                    std::abs(now.omega1 + eps * mode.u1 * mode.omega * s),
                    std::abs(now.omega2 + eps * mode.u2 * mode.omega * s)});
    }

    INFO("normal mode with omega = " << mode.omega << " rad/s");
    CHECK_LT(max_angle_error / eps, relative_tolerance);
    CHECK_LT(max_omega_error / (eps * mode.omega), relative_tolerance);
  }
}

// ---------------------------------------------------------------------------
// Limit m2 / m1 -> 0: simple pendulum
// ---------------------------------------------------------------------------
TEST_CASE("for m2/m1 -> 0 the first bob is an ideal simple pendulum") {
  // A vanishing second mass does not react back on the first rod, which then
  // swings as a simple pendulum of length L1 with the EXACT (non-harmonic)
  // period T = 4 K(k) / omega0. The release angles include 143 degrees, where
  // the period is ~60% longer than the small-angle value: this verifies that
  // sin(theta), not theta, is used in the equations.
  //
  // Deviations: the back-reaction of m2 is O(m2/m1) = 3e-10, RK4 error ~1e-11;
  // the crossing time is found by linear interpolation, whose error is
  // O(dt^3) because theta'' = 0 at the crossing. Measured: <= 2.4e-10 relative
  // error on the quarter period; tolerance 1e-8 (40x margin).
  constexpr double relative_tolerance = 1e-8;
  tu::Params const p{3.0, 1e-9, 0.8, 0.5};

  for (double theta0 : {0.3, tu::pi / 2., 2.5}) {
    pnd::State state{theta0, 0., 0., 0.};
    auto pend = tu::make_pendulum(p, state);

    double const expected = simple_pendulum_quarter_period(theta0, p.L1);
    int const max_steps = tu::steps_for(1.5 * expected, pnd::constants::dt);

    double crossing_time = -1.;
    double previous = pend.state().theta1;
    for (int i = 1; i <= max_steps; ++i) {
      pend.evolution();
      double const current = pend.state().theta1;
      if (previous > 0. && current <= 0.) {
        double const t_previous = (i - 1) * pnd::constants::dt;
        crossing_time =
            t_previous + pnd::constants::dt * previous / (previous - current);
        break;
      }
      previous = current;
    }

    INFO("release angle = " << theta0 << " rad");
    REQUIRE(crossing_time > 0.); // the pendulum did reach the vertical
    CHECK_LT(std::abs(crossing_time - expected) / expected, relative_tolerance);
  }
}

// ---------------------------------------------------------------------------
// Symmetries
// ---------------------------------------------------------------------------
TEST_CASE("mirror symmetry x -> -x is preserved by the evolution") {
  // The system is invariant under reflection through the vertical axis:
  // (theta1, theta2, omega1, omega2) -> -(theta1, theta2, omega1, omega2).
  // The vector field is odd in this change of variables (sin and the
  // velocity-quadratic terms flip sign consistently, cos is even), and
  // floating-point negation is exact, so the mirrored trajectory must be the
  // exact mirror image up to the last-bit behaviour of sin/cos. The tolerance
  // 1e-12 (measured: 0 on this run) is thus a pure rounding-level tolerance.
  auto const p = tu::parameter_sets[1];
  pnd::State a{1.1, -0.4, 0.3, 0.9};
  pnd::State b{-1.1, 0.4, -0.3, -0.9};
  auto pend_a = tu::make_pendulum(p, a);
  auto pend_b = tu::make_pendulum(p, b);

  double max_asymmetry = 0.;
  for (int i = 0; i < 5000; ++i) { // 5 s
    pend_a.evolution();
    pend_b.evolution();
    pnd::State const &x = pend_a.state();
    pnd::State const &y = pend_b.state();
    max_asymmetry = std::max(
        max_asymmetry,
        tu::state_distance(x, {-y.theta1, -y.theta2, -y.omega1, -y.omega2}));
  }
  CHECK_LT(max_asymmetry, 1e-12);
}

// ---------------------------------------------------------------------------
// Scaling laws (dimensional analysis)
// ---------------------------------------------------------------------------
TEST_CASE("only the mass ratio matters: scaling both masses scales E only") {
  // The equations of motion are homogeneous of degree 0 in (m1, m2): the
  // masses cancel except through m2/m1. Multiplying both by c must leave the
  // trajectory unchanged and multiply every energy by c. Only rounding
  // differences are allowed (measured 2e-14 on the state, 4e-15 on E);
  // tolerance 1e-12.
  constexpr double c = 7.3;
  auto const p = tu::parameter_sets[1];
  tu::Params const scaled{c * p.m1, c * p.m2, p.L1, p.L2};

  pnd::State a{1.1, -0.4, 0.3, 0.9};
  pnd::State b = a;
  auto pend_p = tu::make_pendulum(p, a);
  auto pend_scaled = tu::make_pendulum(scaled, b);

  for (int i = 0; i < 3000; ++i) { // 3 s
    pend_p.evolution();
    pend_scaled.evolution();
  }

  CHECK_LT(tu::state_distance(pend_p.state(), pend_scaled.state()), 1e-12);
  double const e_p = pend_p.kinetic() + pend_p.potential();
  double const e_scaled = pend_scaled.kinetic() + pend_scaled.potential();
  CHECK(e_scaled / e_p == doctest::Approx(c).epsilon(1e-12));
}

TEST_CASE(
    "scaling the lengths by 4 rescales time by 2 (dimensional analysis)") {
  // With g fixed, the only time scale is sqrt(L / g). If both lengths are
  // multiplied by lambda, the motion is the same function of t / sqrt(lambda):
  // theta'(t) = theta(t / sqrt(lambda)), omega'(t) = omega(t / sqrt(lambda)) /
  // sqrt(lambda). With lambda = 4 the factor sqrt(lambda) = 2 is a power of
  // two, so a run with dt' = 2 dt and the same number of steps is
  // algebraically self-similar: the results must coincide to rounding level
  // (measured: exactly 0). Tolerance 1e-12.
  tu::Params const p{1.0, 2.0, 0.5, 0.3};
  tu::Params const p4{1.0, 2.0, 4. * 0.5, 4. * 0.3};
  pnd::State a{1.1, -0.4, 0.3, 0.9};
  pnd::State b{1.1, -0.4, 0.3 / 2., 0.9 / 2.};
  auto pend = tu::make_pendulum(p, a);
  auto pend4 = tu::make_pendulum(p4, b);

  for (int i = 0; i < 3000; ++i) {
    pend.evolution(1e-3);
    pend4.evolution(2e-3);
  }

  pnd::State const &x = pend.state();
  pnd::State const &y = pend4.state();
  CHECK_LT(
      tu::state_distance(x, {y.theta1, y.theta2, 2. * y.omega1, 2. * y.omega2}),
      1e-12);
}

TEST_SUITE_END();
