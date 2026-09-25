#ifndef PENDULUM_HPP
#define PENDULUM_HPP

namespace pnd {
namespace constants {
inline constexpr double dt = 0.001;
inline constexpr double g = 9.81;
} // namespace constants

struct vec2 {
  double x;
  double y;

  vec2 operator+(vec2 const &v);
  vec2 operator-(vec2 const &v);
  vec2 operator*(double k); // factor on the right
  vec2 operator/(double k);
  double norm();
  double norm2();
};
vec2 operator*(double k, vec2 const &v); // factor on the left

class Pendulum {
private:
  vec2 origin_;
  double lenght_;
  double angle_;
  double mass_;
  vec2 position_;
  vec2 velocity_;
  vec2 acceleration_;

  void check_lenght();
  void check_angle();
  void check_mass();
  void check_position();
  void check_velocity();
  void check_acceleration();

public:
  Pendulum(vec2 const &origin, double lenght, double angle, double mass);

  vec2 origin();
  double lenght();
  double angle();
  vec2 position();

  void evolution();
};

} // namespace pnd
#endif
