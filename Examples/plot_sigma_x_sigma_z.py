import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import sys

filein = open("params.txt")
Ncopies = int(filein.readline())
Nensemble = int(filein.readline())
ti = float(filein.readline())
tf = float (filein.readline())
dt = float(filein.readline())
print_traj = bool(filein.readline())
Ntraj = int(filein.readline())
dimH = int(filein.readline())
Npoints = int((tf-ti)/dt)+1

t = np.arange(ti,tf,dt)
markevery = int(Npoints/20)

fig, ax = plt.subplots(1,3, figsize=(19,4), sharex=True, sharey=True)

trajectories_p = np.zeros((Ntraj, Npoints))
exact_p = np.zeros(Npoints)
avg_obs_p = np.zeros(Npoints)
err_obs_p = np.zeros(Npoints)
gp = np.zeros(Npoints)
gm = np.zeros(Npoints)
non_P_div = np.zeros(Npoints)
if print_traj == True:
    filein = open("trajectories_p.txt")
f_exact = open("analytic_p.txt")
f_avg = open("average_p.txt")
f_err = open("error_p.txt")
f_g = open("gammas.txt")
for i in range(Npoints):
    exact_p[i] = f_exact.readline()
    avg_obs_p[i] = f_avg.readline()
    err_obs_p[i] = f_err.readline()
    if print_traj == True:
        j = 0
        line = filein.readline()
        for x in line.split():
            trajectories_p[j,i] = x
            j+=1
    g = f_g.readline()
    gm[i] = g.split()[0]
    gp[i] = g.split()[1]
if print_traj == True:
    for i in range(Ntraj):
        ax[0].plot(t, trajectories_p[i,:], alpha=.1)
ax[0].plot(t,exact_p,color='black', label="Exact")
ax[0].errorbar(t,avg_obs_p,err_obs_p, marker='o', markersize=3, color='red', label="Average", errorevery=markevery, markevery=markevery, linewidth=0, elinewidth=1)
axx = ax[2].inset_axes([.1,.1,.4,.3])
axx.plot(t, gm, label=r'$\gamma_-$')
axx.plot(t, gp, '--', label=r'$\gamma_+$')
axx.plot(t, np.zeros(Npoints), color="black", linewidth=.3)
axx.axhline(0, color="black", linewidth=.5)
#axx.legend(loc = "lower left")

trajectories_m = np.zeros((Ntraj, Npoints))
exact_m = np.zeros(Npoints)
avg_obs_m = np.zeros(Npoints)
err_obs_m = np.zeros(Npoints)
exact = np.zeros(Npoints)
avg_obs = np.zeros(Npoints)
err_obs = np.zeros(Npoints)
if print_traj == True:
    filein = open("trajectories_m.txt")
f_exact = open("analytic_m.txt")
f_avg = open("average_m.txt")
f_err = open("error_m.txt")
for i in range(Npoints):
    exact_m[i] = f_exact.readline()
    avg_obs_m[i] = f_avg.readline()
    err_obs_m[i] = f_err.readline()
    exact[i] = .5*(exact_p[i]-exact_m[i])
    avg_obs[i] = .5*(avg_obs_p[i]-avg_obs_m[i])
    err_obs[i] = np.sqrt(err_obs_p[i]**2 + err_obs_m[i]**2)
    if print_traj == True:
        j = 0
        line = filein.readline()
        for x in line.split():
            trajectories_m[j,i] = x
            j+=1
if print_traj == True:
    for i in range(Ntraj):
        ax[1].plot(t, trajectories_m[i,:], alpha=.1)
ax[1].plot(t,exact_m,color='black', label="Exact")
ax[1].errorbar(t,avg_obs_m,err_obs_m, marker='o', markersize=3, color='red', label="Average", errorevery=markevery, markevery=markevery, linewidth=0, elinewidth=1)
ax[2].plot(t,exact,color='black', label="Exact")
ax[2].errorbar(t,avg_obs,err_obs, marker='o', markersize=3, color='red', label="Average", errorevery=markevery, markevery=markevery, linewidth=0, elinewidth=1)

for i in range(3):
    ax[i].legend(loc="lower right")
    ax[i].set_xlabel(r'$t$')
    ax[i].set_ylim([-1.1,1.1])
ax[0].set_title(r"Positive part of $Q_\alpha$: $Q_\alpha^+$")
ax[1].set_title(r"Negative part of $Q_\alpha$: $Q_\alpha^-$")
ax[2].set_title(r"$Q_\alpha = Q_\alpha^+ - Q_\alpha^-$")
axx.set_title("Rates")
if sys.argv.__len__() > 1:
    plt.suptitle(sys.argv[1])

if sys.argv.__len__() > 2:
    ax[0].set_ylabel(sys.argv[2])

if sys.argv.__len__() > 3:
    plt.savefig("Examples/"+sys.argv[3])

plt.show()