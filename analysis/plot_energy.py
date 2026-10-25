import sys

sys.path.insert(0, "../build/Release")

import pandas as pd
import matplotlib.pyplot as plt
import pendulum_sim

N_STEPS = 1e6

theta1 = 60
theta2 = 30
file1 = "../data/test1.csv"
pendulum_sim.simulation(theta1, theta2, file1)
data1 = pd.read_csv(file1)

theta1 = 60.1
theta2 = 30
file2 = "../data/test2.csv"
pendulum_sim.simulation(theta1, theta2, file2)
data2 = pd.read_csv(file2)

theta1 = 60.2
theta2 = 30
file3 = "../data/test3.csv"
pendulum_sim.simulation(theta1, theta2, file3)
data3 = pd.read_csv(file3)

theta1 = 60.3
theta2 = 30
file4 = "../data/test4.csv"
pendulum_sim.simulation(theta1, theta2, file4)
data4 = pd.read_csv(file4)

data = [data1, data2, data3, data4]
fig, ax = plt.subplots(2, 2, figsize=(10, 8))

for file in data:
    ax[0, 0].plot(file["Time"], file["Energy"])
    ax[0, 0].plot(file["Time"], file["Kinetic"])
    ax[0, 0].plot(file["Time"], file["Potential"])
    ax[0, 0].set_title(f"Theta = 60")

plt.show()
