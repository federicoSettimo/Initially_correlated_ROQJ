import matplotlib.pyplot as plt
import numpy as np

with open('error.txt', 'r') as file:
    lines = file.readlines()

Nstates = []
t = []
err = []
for line in lines:
    data = line.split()
    Nstates.append(float(data[0]))
    t.append(float(data[1]))
    err.append(float(data[2]))

with open('error_partial_info.txt', 'r') as file:
    lines = file.readlines()

Nstates_part = []
t_part = []
err_part = []
for line in lines:
    data = line.split()
    Nstates_part.append(float(data[0]))
    t_part.append(float(data[1]))
    err_part.append(float(data[2]))

Nstates_filtered = []
err_filtered = []
Nstates_part_filtered = []
err_part_filtered = []
target_t = 1.3
for i in range(len(t)):
    if t[i] == target_t:
        Nstates_filtered.append(Nstates[i])
        err_filtered.append(err[i])
for i in range(len(t_part)):
    if t_part[i] == target_t:
        Nstates_part_filtered.append(Nstates_part[i])
        err_part_filtered.append(err_part[i])

filein = open('gamma_error.txt')
tf = float(filein.readline())
dt = float(filein.readline())
tt = np.arange(0,tf,dt)
Npoints = int(tf/dt)
g = np.zeros(Npoints)
for i in range(Npoints):
    g[i] = filein.readline()

fig, ax = plt.subplots(1,3, figsize=(16,4), sharex=False, sharey=False)
plt.subplots_adjust(wspace=.4)

plot0 = ax[0].scatter(Nstates, t, c=np.log10(err), cmap='viridis', alpha=0.75)
ax[0].set_xscale('log')
plt.colorbar(plot0, ax=ax[0], label=r'log$_{10}(TD(\rho_{ex}, \rho_{avg}))$')
ax[0].set_xlabel(r'Number of realizations')
ax[0].set_ylabel(r'$t$')
ax[0].set_title('Using full information')

plot1 = ax[1].scatter(Nstates_part, t_part, c=np.log10(err_part), cmap='viridis', alpha=0.75)
ax[1].set_xscale('log')
plt.colorbar(plot1, ax=ax[1], label=r'log$_{10}(TD(\rho_{ex}, \rho_{avg}))$')
ax[1].set_xlabel(r'Number of realizations')
ax[1].set_ylabel(r'$t$')
ax[1].set_title('Using partial information')

ax[2].loglog(Nstates_filtered, err_filtered, label = "Full information")
ax[2].loglog(Nstates_part_filtered, err_part_filtered, '--', label = "Partial information")
ax[2].loglog(Nstates_filtered, 1./np.sqrt(Nstates_filtered)/100, color = "red", linewidth=.5, label = r'$1/\sqrt{N}$')
ax[2].legend(loc = "lower right")
ax[2].set_xlabel(r'Number of realizations')
ax[2].set_ylabel(r'$TD(\rho_{ex}, \rho_{avg}$')
ax[2].set_title(fr'Error for $ t = ${target_t}')
ax[2].set_ylim([2*10**(-7),5*10**(-3)])

axx = ax[2].inset_axes([0.15, 0.1, 0.3, 0.2])
axx.plot(tt, g)
axx.set_title(r'$\gamma$')

plt.savefig('error_2_ME.png')

#plt.show()
