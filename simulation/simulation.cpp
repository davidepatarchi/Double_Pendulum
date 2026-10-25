#include "simulation.hpp"
#include "pendulum.hpp"

#include <fstream>
#include <iomanip>
#include <stdexcept>

void simulation(double theta1, double theta2, const std::string &filename,
                long n_steps, double length1, double length2, double mass1,
                double mass2) {
  theta1 *= std::numbers::pi / 180.;
  theta2 *= std::numbers::pi / 180.;

  std::ofstream file(filename);
  if (!file.is_open())
    throw std::runtime_error("Could not open " + filename);

  file << std::setprecision(17);
  file << "Time,Kinetic,Potential,Total\n";

  pnd::State state{theta1, theta2, 0., 0.};
  pnd::DoublePendulum pend{length1, mass1, length2, mass2, state};

  auto write_row = [&](double t) {
    double K = pend.kinetic();
    double U = pend.potential();
    file << t << ',' << K << ',' << U << ',' << K + U << '\n';
  };

  write_row(0.); // stato iniziale
  for (long i = 0; i < n_steps; ++i) {
    pend.evolution();
    write_row(static_cast<double>(i + 1) * pnd::constants::dt);
  }
}
