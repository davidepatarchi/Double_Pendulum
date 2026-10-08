#!/usr/bin/env python3
"""
Independent reference solutions for the double-pendulum regression tests.

NOTHING here is taken from the C++ code under test.  The equations of motion
are derived symbolically from the Lagrangian (sympy), integrated with a
high-order adaptive explicit Runge-Kutta scheme (DOP853, tight tolerances) and
cross-checked with a different, implicit method (Radau).  The two solvers agree
to ~1e-12, which is far below every tolerance used in the C++ tests.  The
printed numbers are pasted into tests/test_dynamics.cpp.

Conventions (identical to the project, see README):
  theta_i : angle from the DOWNWARD vertical, absolute (not relative),
  omega_i = d theta_i / dt,
  g       : 9.81 m/s^2.

Usage:   python3 generate_reference.py            (needs numpy, scipy, sympy)
"""
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

G = 9.81

# ----------------------------------------------------------------------------
# 1. Equations of motion from the Lagrangian  L = T - V  (symbolic)
# ----------------------------------------------------------------------------
t = sp.symbols("t")
m1, m2, L1, L2, g = sp.symbols("m1 m2 L1 L2 g", positive=True)
th1, th2 = sp.Function("th1")(t), sp.Function("th2")(t)

x1, y1 = L1 * sp.sin(th1), -L1 * sp.cos(th1)
x2, y2 = x1 + L2 * sp.sin(th2), y1 - L2 * sp.cos(th2)

T = sp.Rational(1, 2) * m1 * (x1.diff(t) ** 2 + y1.diff(t) ** 2) \
  + sp.Rational(1, 2) * m2 * (x2.diff(t) ** 2 + y2.diff(t) ** 2)
V = g * (m1 * y1 + m2 * y2)
Lag = T - V

eqs = [sp.diff(Lag.diff(q.diff(t)), t) - Lag.diff(q) for q in (th1, th2)]
acc = sp.solve(eqs, [th1.diff(t, 2), th2.diff(t, 2)], dict=True)[0]

w1, w2, a, b = sp.symbols("w1 w2 a b")
subs = {th1.diff(t): w1, th2.diff(t): w2, th1: a, th2: b}
alpha1 = acc[th1.diff(t, 2)].subs(subs)
alpha2 = acc[th2.diff(t, 2)].subs(subs)
rhs = sp.lambdify((a, b, w1, w2, m1, m2, L1, L2, g),
                  [w1, w2, alpha1, alpha2], "numpy")

# Energy with V = 0 at the stable equilibrium (theta1 = theta2 = 0)
E_expr = (T + V + g * (m1 + m2) * L1 + g * m2 * L2) \
    .subs({th1.diff(t): w1, th2.diff(t): w2}).subs({th1: a, th2: b})
energy = sp.lambdify((a, b, w1, w2, m1, m2, L1, L2, g), E_expr, "numpy")


def f(_t, y, p):
    return rhs(*y, p["m1"], p["m2"], p["L1"], p["L2"], G)


def solve(p, y0, t_eval, method="DOP853", tol=1e-13):
    sol = solve_ivp(f, (0.0, t_eval[-1]), y0, args=(p,), method=method,
                    t_eval=t_eval, rtol=tol, atol=tol)
    assert sol.success
    return sol.y.T


# ----------------------------------------------------------------------------
# 2. Reference cases (same parameters as the C++ regression tests)
# ----------------------------------------------------------------------------
CASES = {
    # Parameters used by main.cpp (L = 150 px / 300 px/m = 0.5 m, m = 10 kg),
    # initial angles of the README example "./progetto 45 60".
    "main_params_45_60": dict(
        p=dict(m1=10.0, m2=10.0, L1=0.5, L2=0.5),
        y0=[np.radians(45.0), np.radians(60.0), 0.0, 0.0],
        times=[0.5, 1.0, 2.0, 3.0]),
    # Asymmetric masses/lengths, large amplitude, non-zero initial velocities:
    # exercises every term of the equations (no accidental symmetry).
    "asymmetric_large_amplitude": dict(
        p=dict(m1=2.0, m2=0.5, L1=1.3, L2=0.7),
        y0=[2.0, -1.0, 0.5, -0.3],
        times=[0.5, 1.0, 2.0]),
    # Both bobs launched with large angular velocities: the angles leave
    # [-pi, pi] several times, exercising the wrap-around in evolution().
    # (Angles printed here are UNWRAPPED; the C++ test compares modulo 2 pi.)
    "rotating_wraps_angles": dict(
        p=dict(m1=1.0, m2=1.0, L1=1.0, L2=1.0),
        y0=[0.0, 0.0, 10.0, 6.0],
        times=[0.5, 1.0, 2.0]),
}


def main():
    for name, c in CASES.items():
        p, y0, times = c["p"], c["y0"], np.array(c["times"])
        ref = solve(p, y0, times, "DOP853", 1e-13)
        alt = solve(p, y0, times, "Radau", 1e-13)
        e0 = energy(*y0, p["m1"], p["m2"], p["L1"], p["L2"], G)
        print(f"\n=== {name}")
        print(f"params = {p},  y0 = {y0}")
        print(f"E(0) = {e0:.17g} J")
        spread = np.max(np.abs(ref - alt))
        print(f"max |DOP853 - Radau|        = {spread:.2e}")
        for tt, row in zip(times, ref):
            e = energy(*row, p["m1"], p["m2"], p["L1"], p["L2"], G)
            print(f"t = {tt:4.1f}:  {{{row[0]:.17g}, {row[1]:.17g}, "
                  f"{row[2]:.17g}, {row[3]:.17g}}},   "
                  f"|dE/E0| = {abs(e - e0) / e0:.1e}")


if __name__ == "__main__":
    main()
