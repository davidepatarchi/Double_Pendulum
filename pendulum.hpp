#ifndef PENDULUM_HPP
#define PENDULUM_HPP

#include <numbers>

namespace pnd {
namespace constants {
inline constexpr double dt = 0.001;
inline constexpr double g = 9.81;
} // namespace constants

struct State {
  double theta1;
  double theta2;
  double omega1;
  double omega2;

  State operator+(State const &state);
  State operator*(double x);
};

class DoublePendulum {
private:
  double lenght1_;
  double mass1_;
  double lenght2_;
  double mass2_;
  State state_;

  void check_lenght();
  void check_mass();
  void check_state();

public:
  DoublePendulum(double lenght1, double mass1, double lenght2, double mass2,
                 State &state);

  double pos_x1(double origin_x1);
  double pos_y1(double origin_y1);
  double pos_x2(double origin_x);
  double pos_y2(double origin_x);

  double &lenght1();
  double &lenght2();
  double &mass1();
  double &mass2();
  State &state();

  State derive(State const &state);
  void evolution();
  double energy();
};

} // namespace pnd
#endif
