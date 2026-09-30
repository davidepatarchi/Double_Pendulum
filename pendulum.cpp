#include "pendulum.hpp"

#include <cmath>
#include <stdexcept>

namespace pnd {
// PENDULUM
void Pendulum::check_mass() {
  if (mass_ <= 0. || isfinite(mass_) == false) {
    throw std::invalid_argument("Error: invalid mass");
  }
}
void Pendulum::check_lenght() {
  if (lenght_ <= 0. || isfinite(lenght_) == false) {
    throw std::invalid_argument("Error: invalid lenght");
  }
}
void Pendulum::check_state() {
  if (isfinite(theta_) == false || isfinite(omega_) == false ||
      isfinite(alpha_) == false) {
    throw std::invalid_argument("Error: invalid state");
  }
}

Pendulum::Pendulum(double lenght, double mass, double theta)
    : lenght_{lenght}, mass_{mass}, theta_{theta} {

  omega_ = 0.;
  alpha_ = 0.;

  check_lenght();
  check_mass();
  check_state();
}

double Pendulum::pos_x(double origin_x) {
  return origin_x + lenght_ * std::sin(theta_);
}
double Pendulum::pos_y(double origin_y) {
  return origin_y - lenght_ * std::cos(theta_);
}

double &Pendulum::lenght() { return lenght_; }
double &Pendulum::mass() { return mass_; }
double &Pendulum::theta() { return theta_; }
double &Pendulum::omega() { return omega_; }
double &Pendulum::alpha() { return alpha_; }

void Pendulum::evolution(Pendulum &pendulum) {
  double alpha1 = alpha_;
  double alpha2 = pendulum.alpha();

  // position
  theta_ +=
      omega_ * constants::dt + 0.5 * alpha_ * constants::dt * constants::dt;
  pendulum.theta() += pendulum.omega() * constants::dt +
                      0.5 * pendulum.alpha() * constants::dt * constants::dt;

  // velocity prediction
  omega_ += alpha_ * constants::dt;
  pendulum.omega() += pendulum.alpha() * constants::dt;

  // acceleration
  double m1 = mass_;
  double m2 = pendulum.mass();
  double L1 = lenght_;
  double L2 = pendulum.lenght();
  double th1 = theta_;
  double th2 = pendulum.theta();
  double w1 = omega_;
  double w2 = pendulum.omega();
  double delta = theta_ - pendulum.theta();

  alpha_ = (-m2 * L1 * w1 * w1 * std::sin(delta) * std::cos(delta) +
            m2 * constants::g * std::sin(th2) * std::cos(delta) -
            m2 * L2 * w2 * w2 * std::sin(delta) -
            (m1 + m2) * constants::g * std::sin(th1)) /
           (L1 * (m1 + m2 * std::sin(delta) * std::sin(delta)));
  pendulum.alpha() =
      ((m1 + m2) * (L1 * w1 * w1 * std::sin(delta) +
                    constants::g * std::sin(th1) * std::cos(delta) -
                    constants::g * std::sin(th2)) +
       m2 * L2 * w2 * w2 * std::sin(delta) * std::cos(delta)) /
      (L2 * (m1 + m2 * std::sin(delta) * std::sin(delta)));

  // velocity
  omega_ += (alpha_ + alpha1) * 0.5 * constants::dt;
  pendulum.omega() += (pendulum.alpha() + alpha2) * 0.5 * constants::dt;

  check_state();
  if (theta_ < -std::numbers::pi) {
    theta_ += 2 * std::numbers::pi;
  }
  if (theta_ > std::numbers::pi) {
    theta_ -= 2 * std::numbers::pi;
  }
  if (pendulum.theta() < -std::numbers::pi) {
    pendulum.theta() += 2 * std::numbers::pi;
  }
  if (pendulum.theta() > std::numbers::pi) {
    pendulum.theta() -= 2 * std::numbers::pi;
  }
}
// PENDULUM

} // namespace pnd