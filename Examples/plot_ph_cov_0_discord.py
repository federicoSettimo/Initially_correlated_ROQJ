import numpy as np
import matplotlib.pyplot as plt

# Read data from avg_0_z.txt
with open("avg_0_z.txt", "r") as f:
    avg_data = [float(line.strip()) for line in f]

# Read data from exact_0_z.txt
with open("exact_0_z.txt", "r") as f:
    exact_data = [float(line.strip()) for line in f]

# Generate time values from 0 to 10 with the same number of points as avg_data
time = np.linspace(0, 5, len(avg_data))

# Plotting
markevery = int(len(avg_data)/30)
plt.plot(time, avg_data, label='Avgerage', marker = 'o', markersize = 3, color = "red", markevery = markevery, linewidth = 0)
plt.plot(time, exact_data, label='Exact', color = "black")

# Add labels and title
plt.xlabel(r'$t$')
plt.ylabel(r'$\operatorname{tr}[X \sigma_z]$')
plt.legend()

# Display the plot
plt.show()
