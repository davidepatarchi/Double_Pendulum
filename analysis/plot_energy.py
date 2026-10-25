import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "build" / "simulation" / "Release"))

import pandas as pd
import matplotlib.pyplot as plt
import pendulum_sim

N_STEPS = 1_005_000

theta1_values = [60.0, 60.1, 60.2, 60.3]
theta2 = 30.0

data = []

for i, theta1 in enumerate(theta1_values, start=1):
    file = str(ROOT / "data" / f"test{i}.csv")
    pendulum_sim.simulation(theta1, theta2, file, N_STEPS)
    data.append(pd.read_csv(file))


fig, ax = plt.subplots(2, 2, figsize=(10, 8))

axes = ax.flatten()

for i, (dataframe, theta1) in enumerate(zip(data, theta1_values)):
    axes[i].plot(dataframe["Time"], dataframe["Total"], label="Total")
    axes[i].plot(dataframe["Time"], dataframe["Kinetic"], label="Kinetic")
    axes[i].plot(dataframe["Time"], dataframe["Potential"], label="Potential")

    axes[i].set_title(f"$\\theta_1 = {theta1}°$, $\\theta_2 = {theta2}°$")
    axes[i].set_xlabel("Time [s]")
    axes[i].set_ylabel("Energy [J]")
    axes[i].legend()
    axes[i].grid()

plt.tight_layout()
plt.show()
