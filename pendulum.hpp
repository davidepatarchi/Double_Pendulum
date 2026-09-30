#ifndef PENDULUM_HPP
#define PENDULUM_HPP

#include <numbers>

namespace pnd {
namespace constants {
inline constexpr double dt = 1e-4;
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
  double lenght_;
  double mass_;
  double theta_;
  double omega_;
  double alpha_;

  void check_lenght();
  void check_mass();
  void check_state();

public:
  Pendulum(double lenght, double mass, double theta);

  double pos_x(double origin_x);
  double pos_y(double origin_y);
  double &lenght();
  double &mass();
  double &theta();
  double &omega();
  double &alpha();

  void evolution(Pendulum &pendulum);
};

} // namespace pnd
#endif
