import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("data/energy.csv")

plt.plot(data["Time"], data["Kinetic"], label="Kinetic")
plt.plot(data["Time"], data["Potential"], label="Potential")
plt.plot(data["Time"], data["Total"], label="Total")

plt.xlabel("Time [s]")
plt.ylabel("Energy [J]")
plt.legend()
plt.grid()

plt.show()
