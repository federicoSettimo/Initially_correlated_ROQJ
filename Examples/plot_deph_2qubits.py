import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(2,3, figsize=(13,8), sharex=True, sharey=True)
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
colors = colors+colors+colors

obs1_avg = []
obs2_avg = []
obs1_ex = []
obs2_ex = []

filein = open("params_0.txt")
tf = float (filein.readline())
dt = float(filein.readline())
Npoints = int((tf)/dt)
t = np.arange(0,tf,dt)
ex_1 = np.zeros(Npoints)
ex_2 = np.zeros(Npoints)
avg_1 = np.zeros(Npoints)
avg_2 = np.zeros(Npoints)
f_ex = open("exact_0.txt")
#f_avg = open("exact_0.txt") # Now I am cheating for Q0, but we already know it is ok (Markov) - juts sth wrong when writing eigsys in c++
f_avg = open("avg_0.txt")
for i in range(Npoints):
    ex = f_ex.readline()
    ex_1[i] = float(ex.split()[0])
    ex_2[i] = float(ex.split()[1])
    avg = f_avg.readline()
    avg_1[i] = float(avg.split()[0])
    avg_2[i] = float(avg.split()[1])
obs1_avg.append(avg_1)
obs2_avg.append(avg_2)
obs1_ex.append(ex_1)
obs2_ex.append(ex_2)
ax[0,0].plot(t, ex_1, color=colors[0], label=r'$Q_0$')
ax[0,0].plot(t, avg_1, marker = 'o', markersize = 3, color=colors[0], markevery = int(Npoints/20), linewidth = 0)
ax[1,0].plot(t, ex_2, color=colors[0], label=r'$Q_0$')
ax[1,0].plot(t, avg_2, marker = 'o', markersize = 3, color=colors[0], markevery = int(Npoints/20), linewidth = 0)

icol = 0
for i in range(4):
    for j in range(i):
        icol = icol+1
        ex_1 = np.zeros(Npoints)
        ex_2 = np.zeros(Npoints)
        avg_1 = np.zeros(Npoints)
        avg_2 = np.zeros(Npoints)
        f_ex = open("exact_"+str(i)+"_"+str(j)+"_x.txt")
        f_avg = open("avg_"+str(i)+"_"+str(j)+"_x.txt")
        for ii in range(Npoints):
            ex = f_ex.readline()
            ex_1[ii] = float(ex.split()[0])
            ex_2[ii] = float(ex.split()[1])
            avg = f_avg.readline()
            avg_1[ii] = float(avg.split()[0])
            avg_2[ii] = float(avg.split()[1])
        obs1_avg.append(avg_1)
        obs2_avg.append(avg_2)
        obs1_ex.append(ex_1)
        obs2_ex.append(ex_2)
        ax[0,0].plot(t, ex_1, color=colors[icol], label=r"$Q_{{{},{}}}^x$".format(i,j))
        ax[0,0].plot(t, avg_1, marker = 'o', markersize = 3, color=colors[icol], markevery = int(Npoints/20), linewidth = 0)
        ax[1,0].plot(t, ex_2, color=colors[icol], label=r"$Q_{{{},{}}}^x$".format(i,j))
        ax[1,0].plot(t, avg_2, marker = 'o', markersize = 3, color=colors[icol], markevery = int(Npoints/20), linewidth = 0)

        ex_1 = np.zeros(Npoints)
        ex_2 = np.zeros(Npoints)
        avg_1 = np.zeros(Npoints)
        avg_2 = np.zeros(Npoints)
        f_ex = open("exact_"+str(i)+"_"+str(j)+"_y.txt")
        f_avg = open("avg_"+str(i)+"_"+str(j)+"_y.txt")
        for ii in range(Npoints):
            ex = f_ex.readline()
            ex_1[ii] = float(ex.split()[0])
            ex_2[ii] = float(ex.split()[1])
            avg = f_avg.readline()
            avg_1[ii] = float(avg.split()[0])
            avg_2[ii] = float(avg.split()[1])
        obs1_avg.append(avg_1)
        obs2_avg.append(avg_2)
        obs1_ex.append(ex_1)
        obs2_ex.append(ex_2)
        ax[0,1].plot(t, ex_1, color=colors[icol-1], label=r"$Q_{{{},{}}}^y$".format(i,j))
        ax[0,1].plot(t, avg_1, marker = 'o', markersize = 3, color=colors[icol-1], markevery = int(Npoints/20), linewidth = 0)
        ax[1,1].plot(t, ex_2, color=colors[icol-1], label=r"$Q_{{{},{}}}^y$".format(i,j))
        ax[1,1].plot(t, avg_2, marker = 'o', markersize = 3, color=colors[icol-1], markevery = int(Npoints/20), linewidth = 0)


w = np.zeros(len(obs1_avg))+1
w[2] = 2
w[7] = 1.3
w[9] = 1.1
w[10] = 1.6
avg = np.zeros(len(t))
ex = np.zeros(len(t))
for i in range(len(obs1_avg)):
    avg += w[i]*obs1_avg[i]
    ex += w[i]*obs1_ex[i]
ax[0,2].plot(t, ex, color=colors[0])
ax[0,2].plot(t, avg, marker = 'o', markersize = 3, color=colors[0], markevery = int(Npoints/20), linewidth = 0)
avg = np.zeros(len(t))
ex = np.zeros(len(t))
for i in range(len(obs1_avg)):
    avg += w[i]*obs2_avg[i]
    ex += w[i]*obs2_ex[i]
ax[1,2].plot(t, ex, color=colors[0])
ax[1,2].plot(t, avg, marker = 'o', markersize = 3, color=colors[0], markevery = int(Npoints/20), linewidth = 0)

w = np.zeros(len(obs1_avg))+1
w[2] = 0
w[3] = 0
w[4] = 1.5
w[5] = 0
w[6] = 0
w[7] = 0
w[9] = 0
w[11] = 1.3
w[12] = 1.2
avg = np.zeros(len(t))
ex = np.zeros(len(t))
for i in range(len(obs1_avg)):
    avg += w[i]*obs1_avg[i]
    ex += w[i]*obs1_ex[i]
ax[0,2].plot(t, ex, color=colors[1])
ax[0,2].plot(t, avg, marker = 'o', markersize = 3, color=colors[1], markevery = int(Npoints/20), linewidth = 0)
avg = np.zeros(len(t))
ex = np.zeros(len(t))
for i in range(len(obs1_avg)):
    avg += w[i]*obs2_avg[i]
    ex += w[i]*obs2_ex[i]
ax[1,2].plot(t, ex, color=colors[1])
ax[1,2].plot(t, avg, marker = 'o', markersize = 3, color=colors[1], markevery = int(Npoints/20), linewidth = 0)


ax[0,0].set_title(r'$\Phi^{\alpha,x}_t(Q_\alpha^x)$')
ax[0,0].set_ylabel(r'$<2|Q_\alpha|0> + <1|Q_\alpha|0>$ (local coherences)')
ax[0,1].set_title(r'$\Phi^{\alpha,y}_t(Q_\alpha^y)$')
ax[0,2].set_title(r'$\sum_\alpha w_\alpha^{s}\Phi^{\alpha,x,y}_t(Q_\alpha^{x,y})$')
ax[1,0].set_ylabel(r'$<3|Q_\alpha|0>$ (non-local coherences)')
ax[1,0].set_xlabel(r'$t$')
ax[1,1].set_xlabel(r'$t$')
ax[1,2].set_xlabel(r'$t$')
for i in range(2):
    for j in range(2):
        ax[i,j].legend(loc="lower right")


plt.savefig("2qubits_B+_Markov.png", dpi=300)
plt.show()
