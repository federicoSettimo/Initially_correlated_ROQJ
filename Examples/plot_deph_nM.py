import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(1,3, figsize=(16,4), sharex=True, sharey=True)

# Reading Phi_0(sigma_0)
filein = open("params_0.txt")
tf = float (filein.readline())
dt = float(filein.readline())
Npoints = int((tf)/dt)
t = np.arange(0,tf,dt)
exact_x_0 = np.zeros(Npoints)
exact_y_0 = np.zeros(Npoints)
avg_x_0 = np.zeros(Npoints)
avg_y_0 = np.zeros(Npoints)
gamma_0 = np.zeros(Npoints)
f_ex = open("exact_0.txt")
f_avg = open("avg_0.txt")
f_g = open("gamma_0.txt")
for i in range(Npoints):
    ex = f_ex.readline()
    exact_x_0[i] = float(ex.split()[0])
    exact_y_0[i] = float(ex.split()[1])
    avg = f_avg.readline()
    avg_x_0[i] = float(avg.split()[0])
    avg_y_0[i] = float(avg.split()[1])
    gamma_0[i] = float(f_g.readline())
ax[0].plot(t, gamma_0, color = "orange", label = r'$0$')
ax[1].plot(t, exact_x_0, color = "orange", label = r'$Q_0=(1-\sum_i\sigma_i)/2$')
ax[1].plot(t, exact_y_0, '--', color = "orange")
ax[1].plot(t, avg_x_0, marker = 'o', markersize = 3, color = "orange", markevery = int(Npoints/20), linewidth = 0)
ax[1].plot(t, avg_y_0, marker = 'x', markersize = 3, color = "orange", markevery = int(Npoints/20), linewidth = 0)

# Reading Phi_x(sigma_x)
filein = open("params_x.txt")
tf = float (filein.readline())
dt = float(filein.readline())
Npoints = int((tf)/dt)
t = np.arange(0,tf,dt)
exact_x_x = np.zeros(Npoints)
exact_y_x = np.zeros(Npoints)
avg_x_x = np.zeros(Npoints)
avg_y_x = np.zeros(Npoints)
gamma_x = np.zeros(Npoints)
f_ex = open("exact_x.txt")
f_avg = open("avg_x.txt")
f_g = open("gamma_x.txt")
for i in range(Npoints):
    ex = f_ex.readline()
    exact_x_x[i] = float(ex.split()[0])
    exact_y_x[i] = float(ex.split()[1])
    avg = f_avg.readline()
    avg_x_x[i] = float(avg.split()[0])
    avg_y_x[i] = float(avg.split()[1])
    gamma_x[i] = float(f_g.readline())
ax[0].plot(t, gamma_x, color = "green", label = r'$x$')
ax[1].plot(t, exact_x_x, color = "green", label = r'$Q_x=\sigma_x/2$')
ax[1].plot(t, exact_y_x, '--', color = "green")
ax[1].plot(t, avg_x_x, marker = 'o', markersize = 3, color = "green", markevery = int(Npoints/20), linewidth = 0)
ax[1].plot(t, avg_y_x, marker = 'x', markersize = 3, color = "green", markevery = int(Npoints/20), linewidth = 0)


# Reading Phi_y(sigma_y)
filein = open("params_y.txt")
tf = float (filein.readline())
dt = float(filein.readline())
Npoints = int((tf)/dt)
t = np.arange(0,tf,dt)
exact_x_y = np.zeros(Npoints)
exact_y_y = np.zeros(Npoints)
avg_x_y = np.zeros(Npoints)
avg_y_y = np.zeros(Npoints)
gamma_y = np.zeros(Npoints)
f_ex = open("exact_y.txt")
f_avg = open("avg_y.txt")
f_g = open("gamma_y.txt")
for i in range(Npoints):
    ex = f_ex.readline()
    exact_x_y[i] = float(ex.split()[0])
    exact_y_y[i] = float(ex.split()[1])
    avg = f_avg.readline()
    avg_x_y[i] = float(avg.split()[0])
    avg_y_y[i] = float(avg.split()[1])
    gamma_y[i] = float(f_g.readline())
ax[0].plot(t, gamma_y, color = "blue", label = r'$y$')
ax[1].plot(t, exact_x_y, color = "blue", label = r'$Q_y=\sigma_y/2$')
ax[1].plot(t, exact_y_y, '--', color = "blue")
ax[1].plot(t, avg_x_y, marker = 'o', markersize = 3, color = "blue", markevery = int(Npoints/20), linewidth = 0)
ax[1].plot(t, avg_y_y, marker = 'x', markersize = 3, color = "blue", markevery = int(Npoints/20), linewidth = 0)


# Reading Phi_z(sigma_z)
filein = open("params_z.txt")
tf = float (filein.readline())
dt = float(filein.readline())
Npoints = int((tf)/dt)
t = np.arange(0,tf,dt)
exact_x_z = np.zeros(Npoints)
exact_y_z = np.zeros(Npoints)
avg_x_z = np.zeros(Npoints)
avg_y_z = np.zeros(Npoints)
gamma_z = np.zeros(Npoints)
f_ex = open("exact_z.txt")
f_avg = open("avg_z.txt")
f_g = open("gamma_z.txt")
for i in range(Npoints):
    ex = f_ex.readline()
    exact_x_z[i] = float(ex.split()[0])
    exact_y_z[i] = float(ex.split()[1])
    avg = f_avg.readline()
    avg_x_z[i] = float(avg.split()[0])
    avg_y_z[i] = float(avg.split()[1])
    gamma_z[i] = float(f_g.readline())
ax[0].plot(t, gamma_z, color = "red", label = r'$z$')
ax[1].plot(t, exact_x_z, color = "red", label = r'$Q_z=\sigma_z/2$')
ax[1].plot(t, exact_y_z, '--', color = "red")
ax[1].plot(t, avg_x_z, marker = 'o', markersize = 3, color = "red", markevery = int(Npoints/20), linewidth = 0)
ax[1].plot(t, avg_y_z, marker = 'x', markersize = 3, color = "red", markevery = int(Npoints/20), linewidth = 0)

# Plotting some states
w0 = 1
wx = 1
wy = 1
wz = 1
exact_x = w0*exact_x_0 + wx*exact_x_x + wy*exact_x_y + wz*exact_x_z;
exact_y = w0*exact_y_0 + wx*exact_y_x + wy*exact_y_y + wz*exact_y_z;
avg_x = w0*avg_x_0 + wx*avg_x_x + wy*avg_x_y + wz*avg_x_z;
avg_y = w0*avg_y_0 + wx*avg_y_x + wy*avg_y_y + wz*avg_y_z;
ax[2].plot(t, exact_x, color = "red", label = r'|11>+|00>')
ax[2].plot(t, exact_y, '--', color = "red")
ax[2].plot(t, avg_x, marker = 'o', markersize = 3, color = "red", markevery = int(Npoints/20), linewidth = 0)
ax[2].plot(t, avg_y, marker = 'x', markersize = 3, color = "red", markevery = int(Npoints/20), linewidth = 0)

wx = 5/3
wy = 1
wz = 4/3
exact_x = w0*exact_x_0 + wx*exact_x_x + wy*exact_x_y + wz*exact_x_z;
exact_y = w0*exact_y_0 + wx*exact_y_x + wy*exact_y_y + wz*exact_y_z;
avg_x = w0*avg_x_0 + wx*avg_x_x + wy*avg_x_y + wz*avg_x_z;
avg_y = w0*avg_y_0 + wx*avg_y_x + wy*avg_y_y + wz*avg_y_z;
ax[2].plot(t, exact_x, color = "green", label = r'|11>+|10>+|01>')
ax[2].plot(t, exact_y, '--', color = "green")
ax[2].plot(t, avg_x, marker = 'o', markersize = 3, color = "green", markevery = int(Npoints/20), linewidth = 0)
ax[2].plot(t, avg_y, marker = 'x', markersize = 3, color = "green", markevery = int(Npoints/20), linewidth = 0)

wx = 1
wy = 5/3
wz = 4/3
exact_x = w0*exact_x_0 + wx*exact_x_x + wy*exact_x_y + wz*exact_x_z;
exact_y = w0*exact_y_0 + wx*exact_y_x + wy*exact_y_y + wz*exact_y_z;
avg_x = w0*avg_x_0 + wx*avg_x_x + wy*avg_x_y + wz*avg_x_z;
avg_y = w0*avg_y_0 + wx*avg_y_x + wy*avg_y_y + wz*avg_y_z;
ax[2].plot(t, exact_x, color = "blue", label = r'|11>+i|10>+i|01>')
ax[2].plot(t, exact_y, '--', color = "blue")
ax[2].plot(t, avg_x, marker = 'o', markersize = 3, color = "blue", markevery = int(Npoints/20), linewidth = 0)
ax[2].plot(t, avg_y, marker = 'x', markersize = 3, color = "blue", markevery = int(Npoints/20), linewidth = 0)

wx = 1
wy = 17/11
wz = 20/11
exact_x = w0*exact_x_0 + wx*exact_x_x + wy*exact_x_y + wz*exact_x_z;
exact_y = w0*exact_y_0 + wx*exact_y_x + wy*exact_y_y + wz*exact_y_z;
avg_x = w0*avg_x_0 + wx*avg_x_x + wy*avg_x_y + wz*avg_x_z;
avg_y = w0*avg_y_0 + wx*avg_y_x + wy*avg_y_y + wz*avg_y_z;
ax[2].plot(t, exact_x, color = "orange", label = r'3|11>+i|10>+i|01>')
ax[2].plot(t, exact_y, '--', color = "orange")
ax[2].plot(t, avg_x, marker = 'o', markersize = 3, color = "orange", markevery = int(Npoints/20), linewidth = 0)
ax[2].plot(t, avg_y, marker = 'x', markersize = 3, color = "orange", markevery = int(Npoints/20), linewidth = 0)


for i in range(3):
    ax[i].set_xlabel(r'$t$')
    ax[i].legend(loc = "lower right")
ax[0].set_title(r'$\gamma^\alpha(t)$')
ax[1].set_title(r'tr$[\Phi^\alpha_t(Q_\alpha)\sigma_{x,(y)}]$ solid (dashed)')
ax[2].set_title(r'tr$[\sum_\alpha w_\alpha\Phi^\alpha_t(Q_\alpha)\sigma_{x,(y)}]$ solid (dashed)')

plt.savefig("Initially_correlated_nM_partial_info.png", dpi = 300)
plt.show()