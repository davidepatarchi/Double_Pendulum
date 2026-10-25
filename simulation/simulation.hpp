#pragma once
#include <string>

void simulation(double theta1, double theta2, const std::string &filename,
                long n_steps, double length1 = 0.5, double length2 = 0.5,
                double mass1 = 10., double mass2 = 10.);