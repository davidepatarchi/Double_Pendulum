#include "pendulum.hpp"

#include <cmath>
#include <stdexcept>

namespace pnd {
// STATE
State State::operator+(State const &state) {
  return {theta1 + state.theta1, theta2 + state.theta2, omega1 + state.omega1,
          omega2 + state.omega2};
}
State State::operator*(double x) {
  return {theta1 * x, theta2 * x, omega1 * x, omega2 * x};
}
// STATE

// DOUBLE PENDULUM
void DoublePendulum::check_mass() {
  if (mass1_ <= 0. || isfinite(mass1_) == false) {
    throw std::invalid_argument("Error: invalid mass1");
  }
  if (mass2_ <= 0. || isfinite(mass2_) == false) {
    throw std::invalid_argument("Error: invalid mass2");
  }
}
void DoublePendulum::check_lenght() {
  if (lenght1_ <= 0. || isfinite(lenght1_) == false) {
    throw std::invalid_argument("Error: invalid lenght1");
  }
  if (lenght2_ <= 0. || isfinite(lenght2_) == false) {
    throw std::invalid_argument("Error: invalid lenght2");
  }
}
void DoublePendulum::check_state() {
  if (isfinite(state_.theta1) == false || isfinite(state_.theta2) == false ||
      isfinite(state_.omega1) == false || isfinite(state_.omega2) == false) {
    throw std::invalid_argument("Error: invalid state");
  }
}

DoublePendulum::DoublePendulum(double lenght1, double mass1, double lenght2,
                               double mass2, State &state)
    : lenght1_{lenght1}, mass1_{mass1}, lenght2_{lenght2}, mass2_{mass2},
      state_{state} {
  check_lenght();
  check_mass();
  check_state();
}

double DoublePendulum::pos_x1(double origin_x) {
  return origin_x + lenght1_ * std::sin(state_.theta1);
}
double DoublePendulum::pos_y1(double origin_y) {
  return origin_y - lenght1_ * std::cos(state_.theta1);
}
double DoublePendulum::pos_x2(double origin_x) {
  return origin_x + lenght2_ * std::sin(state_.theta2);
}
double DoublePendulum::pos_y2(double origin_y) {
  return origin_y - lenght2_ * std::cos(state_.theta2);
}

double &DoublePendulum::lenght1() { return lenght1_; }
double &DoublePendulum::lenght2() { return lenght2_; }
double &DoublePendulum::mass1() { return mass1_; }
double &DoublePendulum::mass2() { return mass2_; }
State &DoublePendulum::state() { return state_; }

State DoublePendulum::derive(State const &state) {
  double m1 = mass1_;
  double m2 = mass2_;
  double L1 = lenght1_;
  double L2 = lenght2_;
  double th1 = state.theta1;
  double th2 = state.theta2;
  double w1 = state.omega1;
  double w2 = state.omega2;
  double delta = th1 - th2;
  State derived;

  derived.theta1 = w1;
  derived.theta2 = w2;
  derived.omega1 = (-m2 * L1 * w1 * w1 * std::sin(delta) * std::cos(delta) +
                    m2 * constants::g * std::sin(th2) * std::cos(delta) -
                    m2 * L2 * w2 * w2 * std::sin(delta) -
                    (m1 + m2) * constants::g * std::sin(th1)) /
                   (L1 * (m1 + m2 * std::sin(delta) * std::sin(delta)));
  derived.omega2 =
      ((m1 + m2) * (L1 * w1 * w1 * std::sin(delta) +
                    constants::g * std::sin(th1) * std::cos(delta) -
                    constants::g * std::sin(th2)) +
       m2 * L2 * w2 * w2 * std::sin(delta) * std::cos(delta)) /
      (L2 * (m1 + m2 * std::sin(delta) * std::sin(delta)));

  return derived;
}

void DoublePendulum::evolution() {
  // RK4 evolution
  State k1 = derive(state_);
  State k2 = derive(state_ + k1 * (constants::dt / 2));
  State k3 = derive(state_ + k2 * (constants::dt / 2));
  State k4 = derive(state_ + k3 * constants::dt);

  state_ = state_ + (k1 + k2 * 2. + k3 * 2. + k4) * (constants::dt / 6);

  check_state();
  if (state_.theta1 < -std::numbers::pi) {
    state_.theta1 += 2 * std::numbers::pi;
  }
  if (state_.theta1 > std::numbers::pi) {
    state_.theta1 -= 2 * std::numbers::pi;
  }
  if (state_.theta2 < -std::numbers::pi) {
    state_.theta2 += 2 * std::numbers::pi;
  }
  if (state_.theta2 > std::numbers::pi) {
    state_.theta2 -= 2 * std::numbers::pi;
  }
}

double DoublePendulum::kinetic() {
  double kinetic1 = 0.5 * (mass1_ + mass2_) * lenght1_ * lenght1_ *
                    state_.omega1 * state_.omega1;
  double kinetic2 =
      0.5 * mass2_ * lenght2_ * lenght2_ * state_.omega2 * state_.omega2;
  double kinetic12 = mass2_ * lenght1_ * lenght2_ * state_.omega1 *
                     state_.omega2 * std::cos(state_.theta1 - state_.theta2);

  return kinetic1 + kinetic2 + kinetic12;
}
double DoublePendulum::potential() {
  double potential1 = (mass1_ + mass2_) * constants::g * lenght1_ *
                      (1 - std::cos(state_.theta1));
  double potential2 =
      mass2_ * constants::g * lenght2_ * (1 - std::cos(state_.theta2));

  return potential1 + potential2;
}
// PENDULUM

} // namespace pnd