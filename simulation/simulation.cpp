#include "pendulum.hpp"
#include <fstream>
#include <iostream>
#include <string>

std::ofstream simulation(double theta1, double theta2, std::string filename,
                         float lenght1 = 0.5, float lenght2 = 0.5,
                         double mass1 = 10., double mass2 = 10.,
                         double dt = pnd::constants::dt) {
  // setting initial variables
  double T{0.};
  double K{0.};
  double U{0.};
  double E{0.};

  // creating csv file w columns names: Time, Kinetic, Potential, Total
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Could not open energy.csv\n";
  }
  file << "Time,Kinetic,Potential,Total\n";

  // creating the double pendulum
  pnd::State state{theta1, theta2, 0., 0.};
  pnd::DoublePendulum pend{lenght1, mass1, lenght2, mass2, state};

  // evolving
  for (long int i{0}; i < 1e9; ++i) {
    pend.evolution();

    T += dt;
    K = pend.kinetic();
    U = pend.potential();
    E = K + U;
    file << T << ',' << K << ',' << U << ',' << E << '\n';
  }

  file.close();
  return file;
}
