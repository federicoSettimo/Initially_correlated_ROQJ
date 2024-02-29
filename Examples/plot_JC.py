import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(1,4, figsize=(25,4), sharex=True, sharey=True)

# Reading Phi_x(sigma_x)
filein = open("params.txt")
tf = float (filein.readline())
dt = float(filein.readline())
Delta = float(filein.readline())
g = float(filein.readline())
n0 = int(filein.readline())
nx = int(filein.readline())
ny = int(filein.readline())
nz = int(filein.readline())
Npoints = int((tf)/dt)
t = np.arange(0,tf,dt)



exact_x_x = np.zeros(Npoints)
exact_y_x = np.zeros(Npoints)
exact_z_x = np.zeros(Npoints)
avg_x_x = np.zeros(Npoints)
avg_y_x = np.zeros(Npoints)
avg_z_x = np.zeros(Npoints)
gamma_x = np.zeros(Npoints)
b_x = np.zeros(Npoints)
f_ex = open("exact_x.txt")
f_avg = open("avg_x.txt")
f_g = open("gamma_x.txt")
for i in range(Npoints):
    ex = f_ex.readline()
    exact_x_x[i] = float(ex.split()[0])
    exact_y_x[i] = float(ex.split()[1])
    exact_y_x[i] = float(ex.split()[2])
    avg = f_avg.readline()
    avg_x_x[i] = float(avg.split()[0])
    avg_y_x[i] = float(avg.split()[1])
    avg_y_x[i] = float(avg.split()[2])
    gg = f_g.readline()
    gamma_x[i] = float(gg.split()[0])
    b_x[i] = float(gg.split()[1])
ax[0].plot(t, gamma_x, color = "green", label = r'$\gamma^x$')
ax[0].plot(t, b_x, '--', color = "green", label = r'$b^x$')
ax[1].plot(t, exact_x_x, color = "green", label = r'$\sigma_x$')
ax[1].plot(t, avg_x_x, marker = 'o', markersize = 3, color = "green", markevery = int(Npoints/20), linewidth = 0)
ax[2].plot(t, exact_y_x, color = "green", label = r'$\sigma_x$')
ax[2].plot(t, avg_y_x, marker = 'o', markersize = 3, color = "green", markevery = int(Npoints/20), linewidth = 0)
ax[3].plot(t, exact_z_x, color = "green", label = r'$\sigma_x$')
ax[3].plot(t, avg_z_x, marker = 'o', markersize = 3, color = "green", markevery = int(Npoints/20), linewidth = 0)



 # Reading Phi_y(sigma_y)
exact_x_y = np.zeros(Npoints)
exact_y_y = np.zeros(Npoints)
exact_z_y = np.zeros(Npoints)
avg_x_y = np.zeros(Npoints)
avg_y_y = np.zeros(Npoints)
avg_z_y = np.zeros(Npoints)
gamma_y = np.zeros(Npoints)
b_y = np.zeros(Npoints)
f_ex = open("exact_y.txt")
f_avg = open("avg_y.txt")
f_g = open("gamma_y.txt")
for i in range(Npoints):
    ex = f_ex.readline()
    exact_x_y[i] = float(ex.split()[0])
    exact_y_y[i] = float(ex.split()[1])
    exact_z_y[i] = float(ex.split()[2])
    avg = f_avg.readline()
    avg_x_y[i] = float(avg.split()[0])
    avg_y_y[i] = float(avg.split()[1])
    avg_z_y[i] = float(avg.split()[2])
    gg = f_g.readline()
    gamma_y[i] = float(gg.split()[0])
    b_y[i] = float(gg.split()[1])
ax[0].plot(t, gamma_y, color = "blue", label = r'$\gamma^y$')
ax[0].plot(t, b_y, '--', color = "blue", label = r'$b^y$')
ax[1].plot(t, exact_x_y, color = "blue", label = r'$\sigma_y$')
ax[1].plot(t, avg_x_y, marker = 'o', markersize = 3, color = "blue", markevery = int(Npoints/20), linewidth = 0)
ax[2].plot(t, exact_y_y, color = "blue", label = r'$\sigma_y$')
ax[2].plot(t, avg_y_y, marker = 'o', markersize = 3, color = "blue", markevery = int(Npoints/20), linewidth = 0)
ax[3].plot(t, exact_z_y, color = "blue", label = r'$\sigma_y$')
ax[3].plot(t, avg_z_y, marker = 'o', markersize = 3, color = "blue", markevery = int(Npoints/20), linewidth = 0)



 # Reading Phi_z(sigma_z)
exact_x_z = np.zeros(Npoints)
exact_y_z = np.zeros(Npoints)
exact_z_z = np.zeros(Npoints)
avg_x_z = np.zeros(Npoints)
avg_y_z = np.zeros(Npoints)
avg_z_z = np.zeros(Npoints)
gamma_z = np.zeros(Npoints)
f_ex = open("exact_z.txt")
f_avg = open("avg_z.txt")
f_g = open("gamma_z.txt")
for i in range(Npoints):
    ex = f_ex.readline()
    exact_x_z[i] = float(ex.split()[0])
    exact_y_z[i] = float(ex.split()[1])
    exact_z_z[i] = float(ex.split()[2])
    avg = f_avg.readline()
    avg_x_z[i] = float(avg.split()[0])
    avg_y_z[i] = float(avg.split()[1])
    avg_z_z[i] = float(avg.split()[2])
    gamma_z[i] = float(f_g.readline())
ax[0].plot(t, gamma_z, color = "red", label = r'$\gamma_-^z$')
ax[1].plot(t, exact_x_z, color = "red", label = r'$\sigma_z$')
ax[1].plot(t, avg_x_z, marker = 'o', markersize = 3, color = "red", markevery = int(Npoints/20), linewidth = 0)
ax[2].plot(t, exact_y_z, color = "red", label = r'$\sigma_z$')
ax[2].plot(t, avg_y_z, marker = 'o', markersize = 3, color = "red", markevery = int(Npoints/20), linewidth = 0)
ax[3].plot(t, exact_z_z, color = "red", label = r'$\sigma_z$')
ax[3].plot(t, avg_z_z, marker = 'o', markersize = 3, color = "red", markevery = int(Npoints/20), linewidth = 0)




# Reading Phi_0(sigma_0)
exact_x_0 = np.zeros(Npoints)
exact_y_0 = np.zeros(Npoints)
exact_z_0 = np.zeros(Npoints)
avg_x_0 = np.zeros(Npoints)
avg_y_0 = np.zeros(Npoints)
avg_z_0 = np.zeros(Npoints)
gamma_p_0 = np.zeros(Npoints)
gamma_m_0 = np.zeros(Npoints)
beta = np.zeros(Npoints)
gamma_z = np.zeros(Npoints)
f_ex = open("exact_0.txt")
f_avg = open("avg_0.txt")
f_g = open("gamma_0.txt")
for i in range(Npoints):
    ex = f_ex.readline()
    exact_x_0[i] = float(ex.split()[0])
    exact_y_0[i] = float(ex.split()[1])
    exact_z_0[i] = float(ex.split()[2])
    avg = f_avg.readline()
    avg_x_0[i] = float(avg.split()[0])
    avg_y_0[i] = float(avg.split()[1])
    avg_z_0[i] = float(avg.split()[2])
    gg = f_g.readline()
    gamma_p_0[i] = float(gg.split()[0])
    gamma_m_0[i] = float(gg.split()[1])
    beta[i] = float(gg.split()[2])
    gamma_z[i] = float(gg.split()[3]) + .5*float(gg.split()[4])
ax[0].plot(t, gamma_p_0, color = "orange", label = r'$\gamma_+^0$')
ax[0].plot(t, gamma_m_0, '--', color = "orange", label = r'$\gamma_-^0$')
ax[0].plot(t, gamma_z, '-.', color = "orange", label = r'$\gamma_z^\pm + \frac{1}{2}\gamma_+^\pm$')
ax[0].plot(t, beta, color = "orange", label = r'$b^0$', linewidth = .5)
ax[1].plot(t, .5*exact_x_0, color = "orange", label = r'$Q_0$')
ax[1].plot(t, .5*avg_x_0, marker = 'o', markersize = 3, color = "orange", markevery = int(Npoints/20), linewidth = 0)
ax[2].plot(t, .5*exact_y_0, color = "orange", label = r'$Q_0$')
ax[2].plot(t, .5*avg_y_0, marker = 'o', markersize = 3, color = "orange", markevery = int(Npoints/20), linewidth = 0)
ax[3].plot(t, .5*exact_z_0, color = "orange", label = r'$Q_0$')
ax[3].plot(t, .5*avg_z_0, marker = 'o', markersize = 3, color = "orange", markevery = int(Npoints/20), linewidth = 0)


for i in range(4):
    ax[i].set_xlabel(r'$t$')
    ax[i].legend(loc = "lower left")
ax[0].set_title(r'Functions of the ME')
ax[1].set_title(r'$<1|Q_\alpha(t)|1>$')
ax[2].set_title(r'Re$<0|Q_\alpha(t)|1>$')
ax[3].set_title(r'Im$<0|Q_\alpha(t)|1>$')

plt.suptitle(r"JC, $\Delta = {{{}}}$".format(Delta)+r", $g = {{{}}}$".format(g)+r". Env $\rho_\alpha=|n_\alpha><n_\alpha|$, $n_0 = {{{}}}$".format(n0)+r", $n_x={{{}}}$".format(nx)+r"$, n_y={{{}}}$".format(ny)+r", $n_z={{{}}}$".format(nz))

plt.savefig("JC_n_env.png", dpi = 300)
plt.show()