# Double Pendulum

A 2D simulation of a double pendulum implemented in C++ using SFML.

The project numerically solves the equations of motion of a double pendulum and provides a graphical representation of its evolution. The initial configuration can be set through command-line arguments or interactively by dragging the masses with the mouse.

---

## 1. Compilation, Execution and Data Analysis

### Requirements

The project requires:

* C++20
* CMake
* Ninja
* SFML 3
* A C++ compiler supporting C++20

### Compilation

Clone the repository and move into the project directory:

```bash
git clone https://github.com/davidepatarchi/Double_Pendulum.git
cd Double_Pendulum
```

Configure the project using CMake:

```bash
cmake -S . -B build -G"Ninja Multi-Config"
```

Compile the project in Debug mode:

```bash
cmake --build build --config Debug
```

or in Release mode:

```bash
cmake --build build --config Release
```

### Execution

To start the simulation:

```bash
./build/Debug/progetto
```

or:

```bash
./build/Release/progetto
```

The simulation opens an `800 × 800` window containing the double pendulum.

The two masses can be dragged with the left mouse button to set the initial configuration. When the mouse button is released, the numerical evolution starts.

Initial angles can also be specified from the command line:

```bash
./build/Release/progetto 45 60
```

where:

* `45` is the initial angle of the first pendulum, $\theta_1$;
* `60` is the initial angle of the second pendulum, $\theta_2$.

The angles are measured with respect to the downward vertical direction.

Pressing `Space` resets the system to its equilibrium configuration.

### Data analysis

During the simulation, the total mechanical energy is calculated at every time step.

The collected data can be used to study:

* total mechanical energy;
* energy conservation;
* numerical errors introduced by the integration algorithm;
* dependence of the motion on the initial conditions;
* chaotic behaviour of the double pendulum.

A typical analysis consists of plotting the total energy as a function of time and studying its deviation from the initial value.

---

## 2. Double Pendulum Evolution

### Physical model

The system consists of two point masses connected by two rigid, massless rods.

The parameters used by the simulation are:

* $m_1$, $m_2$: masses of the two pendulums;
* $L_1$, $L_2$: lengths of the two rods;
* $\theta_1$, $\theta_2$: angular positions;
* $\omega_1 = \dot{\theta}_1$, $\omega_2 = \dot{\theta}_2$: angular velocities;
* $\alpha_1 = \ddot{\theta}_1$, $\alpha_2 = \ddot{\theta}_2$: angular accelerations;
* $g$: gravitational acceleration.

The simulation uses:

```text
g = 9.81 m/s²
dt = 10⁻⁴ s
```

The state of the system is therefore described by the four variables

$$
(\theta_1,\theta_2,\omega_1,\omega_2).
$$

### Position of the masses

The position of the first mass, measured from the fixed pivot, is

$$
x_1 = L_1\sin\theta_1,
$$

$$
y_1 = -L_1\cos\theta_1.
$$

The position of the second mass is

$$
x_2 = x_1 + L_2\sin\theta_2,
$$

$$
y_2 = y_1 - L_2\cos\theta_2.
$$

These relations convert the angular coordinates of the pendulums into Cartesian coordinates for rendering.

### Equations of motion

The motion of the two pendulums is coupled: the acceleration of one mass depends on the position and velocity of the other.

Defining

$$
\Delta = \theta_1-\theta_2,
$$

the angular accelerations are

$$
\alpha_1 =
\frac{
-m_2L_1\omega_1^2\sin\Delta\cos\Delta
+m_2g\sin\theta_2\cos\Delta
-m_2L_2\omega_2^2\sin\Delta
-(m_1+m_2)g\sin\theta_1
}{
L_1\left(m_1+m_2\sin^2\Delta\right)
},
$$

and

$$
\alpha_2 =
\frac{
(m_1+m_2)
\left[
L_1\omega_1^2\sin\Delta
+g\sin\theta_1\cos\Delta
-g\sin\theta_2
\right]
+m_2L_2\omega_2^2\sin\Delta\cos\Delta
}{
L_2\left(m_1+m_2\sin^2\Delta\right)
}.
$$

These equations describe the coupled dynamics of the ideal double pendulum, taking into account gravity and the interaction between the two pendulums.

### Numerical integration

The equations of motion are solved numerically because, for arbitrary initial conditions, the double pendulum does not generally have a simple analytical solution.

The simulation uses the **Velocity Verlet algorithm**.

The angular position is updated according to

$$
\theta(t+\Delta t)
=
\theta(t)
+
\omega(t)\Delta t
+
\frac{1}{2}\alpha(t)\Delta t^2.
$$

The angular velocity is first predicted using

$$
\omega^*
=
\omega(t)+\alpha(t)\Delta t.
$$

The new acceleration is then calculated using the updated state.

Finally, the angular velocity is corrected using

$$
\omega(t+\Delta t)
=
\omega(t)
+
\frac{\alpha(t)+\alpha(t+\Delta t)}{2}\Delta t.
$$

This procedure is applied to both pendulums at every time step.

The Velocity Verlet algorithm is particularly suitable for mechanical systems because it generally provides good long-term numerical behaviour and is less prone to energy drift than simple explicit Euler integration.

### Mechanical energy

The total mechanical energy is

$$
E=T+V.
$$

The kinetic energy is

$$
T =
\frac12m_1L_1^2\omega_1^2
+
\frac12m_2L_2^2\omega_2^2
+
m_2L_1L_2\omega_1\omega_2
\cos(\theta_1-\theta_2).
$$

The first two terms correspond to the individual kinetic contributions, while the third term is the coupling term between the two pendulums.

The potential energy is defined so that

$$
V=0
$$

when both pendulums are in their equilibrium configuration.

Therefore,

$$
V =
(m_1+m_2)gL_1(1-\cos\theta_1)
+
m_2gL_2(1-\cos\theta_2).
$$

The total energy is consequently

$$
E =
\frac12m_1L_1^2\omega_1^2
+
\frac12m_2L_2^2\omega_2^2
+
m_2L_1L_2\omega_1\omega_2
\cos(\theta_1-\theta_2)
+
(m_1+m_2)gL_1(1-\cos\theta_1)
+
m_2gL_2(1-\cos\theta_2).
$$

For an ideal double pendulum without dissipative forces, the exact mechanical energy is conserved:

$$
E(t)=E(0).
$$

In the numerical simulation, small variations can occur due to numerical errors associated with the finite time step and the integration algorithm.

---

## 3. Data Analysis

The purpose of the data analysis is to study the numerical behaviour of the simulation and verify the physical properties of the system.

### Energy conservation

The main quantity analysed is the total mechanical energy:

$$
E(t)=T(t)+V(t).
$$

For an ideal double pendulum,

$$
E(t)=E(0).
$$

The numerical result can therefore be studied by plotting the total energy as a function of time.

The relative energy error can be defined as

$$
\varepsilon_E(t)
=
\frac{|E(t)-E(0)|}{|E(0)|}.
$$

This quantity provides a measure of the numerical error introduced by the integration.

A small relative energy error indicates that the numerical solution remains close to the expected conservation law.

### Numerical stability

The data analysis can also be used to investigate the dependence of the simulation on the time step.

By repeating the simulation with different values of $\Delta t$, it is possible to compare the corresponding energy errors.

For example:

```text
dt = 10⁻²
dt = 10⁻³
dt = 10⁻⁴
dt = 10⁻⁵
```

A smaller time step generally reduces the local numerical error, at the cost of increasing the computational time.

### Chaotic behaviour

The double pendulum is a classical example of a chaotic dynamical system.

Small changes in the initial conditions can produce significantly different trajectories after a sufficiently long time.

This can be investigated by running two simulations with slightly different initial conditions:

$$
(\theta_1,\theta_2)
$$

and

$$
(\theta_1+\delta\theta_1,\theta_2+\delta\theta_2).
$$

The trajectories can then be compared to study the sensitivity of the system to its initial conditions.

This behaviour is one of the main reasons why the double pendulum is interesting from both a physical and numerical point of view.

---

## Project Structure

```text
Double_Pendulum/
│
├── CMakeLists.txt
├── main.cpp
├── pendulum.cpp
├── pendulum.hpp
├── render.cpp
├── render.hpp
└── README.md
```

### Main components

* `main.cpp` handles the simulation loop, user input and program execution.
* `pendulum.hpp` contains the definition of the physical model and simulation constants.
* `pendulum.cpp` implements the equations of motion, numerical integration and energy calculation.
* `render.hpp` and `render.cpp` contain the graphical rendering functions.
* `CMakeLists.txt` contains the CMake build configuration.

---

## References

The physical model is based on the classical ideal double pendulum with point masses and massless rods.

The equations of motion are integrated numerically using the **Velocity Verlet** algorithm.

The project combines concepts from:

* classical mechanics;
* numerical methods;
* computational physics;
* C++ programming;
* data analysis.
