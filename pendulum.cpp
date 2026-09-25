#include "pendulum.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace pnd {

// VEC2
vec2 vec2::operator+(vec2 const &v) { return {x + v.x, y + v.y}; }
vec2 vec2::operator-(vec2 const &v) { return {x - v.x, y - v.y}; }
vec2 vec2::operator*(double k) { return {x * k, y * k}; }
vec2 operator*(double k, vec2 const &v) { return {v.x * k, v.y * k}; }
vec2 vec2::operator/(double k) {
  if (k == 0.) {
    throw std::invalid_argument("Error: division by 0");
  }
  double val_x = x / k;
  double val_y = y / k;
  if (std::isfinite(val_x) == false || std::isfinite(val_y) == false) {
    throw std::invalid_argument("Error: impossible division");
  }
  return {val_x, val_y};
}
double vec2::norm() { return std::sqrt(x * x + y * y); }
double vec2::norm2() { return x * x + y * y; }
// VEC2

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
void Pendulum::check_angle() {
  if (isfinite(angle_) == false) {
    throw std::invalid_argument("Error: invalid angle");
  }
}
void Pendulum::check_position() {
  double distance = (position_ - origin_).norm();
  if (std::abs(distance - lenght_) > 0.5) {
    throw std::invalid_argument("Error: point out of trajectory");
  }
}
void Pendulum::check_velocity() {
  if (isfinite(velocity_.x) == false || isfinite(velocity_.y) == false) {
    throw std::invalid_argument("Error: invalid velocity");
  }
}
void Pendulum::check_acceleration() {
  if (isfinite(acceleration_.x) == false ||
      isfinite(acceleration_.y) == false) {
    throw std::invalid_argument("Error: invalid acceleration");
  }
}

Pendulum::Pendulum(vec2 const &origin, double lenght, double angle, double mass)
    : origin_{origin}, lenght_{lenght}, angle_{angle}, mass_{mass} {
  position_ =
      origin_ + vec2{vec2{lenght * std::cos(angle), lenght * std::sin(angle)}};
  velocity_ = {vec2{0., 0.}};
  acceleration_ = {vec2{0., 0.}};
  check_lenght();
  check_angle();
  check_mass();
  check_position();
  check_velocity();
  check_acceleration();
}

vec2 Pendulum::origin() { return origin_; }
double Pendulum::lenght() { return lenght_; }
double Pendulum::angle() { return angle_; }
vec2 Pendulum::position() { return position_; }

void Pendulum::evolution() {
  vec2 direction = origin_ - position_;

  // position
  position_ = position_ + velocity_ * constants::dt +
              0.5 * acceleration_ * constants::dt * constants::dt;

  // angle
  angle_ = std::atan(direction.y / direction.x);

  // acceleration
  vec2 init_acc = acceleration_;
  vec2 g = {0, -constants::g};
  acceleration_ =
      (velocity_.norm2() * direction / (lenght_ * lenght_) + g) / mass_;

  // velocity
  velocity_ = (acceleration_ + init_acc) * 0.5 * constants::dt;
}
// PENDULUM

} // namespace pnd