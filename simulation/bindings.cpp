#include "simulation.hpp"
#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(pendulum_sim, m) {
  m.doc() = "Double pendulum simulation";

  m.def("simulation", &simulation, py::arg("theta1"), py::arg("theta2"),
        py::arg("filename"), py::arg("n_steps"), py::arg("length1") = 0.5,
        py::arg("length2") = 0.5, py::arg("mass1") = 10.,
        py::arg("mass2") = 10., py::call_guard<py::gil_scoped_release>());
}