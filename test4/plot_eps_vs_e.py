import matplotlib.pyplot as plt
from matplotlib import cm
import csv
import numpy as np
import pandas as pd

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 20
plt.rcParams['axes.linewidth'] = 2


# Generate 2 colors from the 'tab10' colormap
colors = cm.get_cmap('tab10', 2)
ax = plt.subplot(111)
ax.xaxis.set_tick_params(which='major', size=7, width=2, direction='in')
ax.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in')
ax.yaxis.set_tick_params(which='major', size=7, width=2, direction='in')
ax.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in')


def errorPerEpsilon(dt):
	dt = str(dt)

	c_data_8 = pd.read_csv('test4/data/' + dt + 'eps8.csv').values
	pos_c_8 = c_data_8[:,:3]
	vel_c_8 = c_data_8[:,3:6]
	time_c_8 = c_data_8[:,6]

	c_data_7 = pd.read_csv('test4/data/' + dt + 'eps7.csv').values
	pos_c_7 = c_data_7[:,:3]
	vel_c_7 = c_data_7[:,3:6]
	time_c_7 = c_data_7[:,6]

	c_data_6 = pd.read_csv('test4/data/' + dt + 'eps6.csv').values
	pos_c_6 = c_data_6[:,:3]
	vel_c_6 = c_data_6[:,3:6]
	time_c_6 = c_data_6[:,6]

	c_data_5 = pd.read_csv('test4/data/' + dt + 'eps5.csv').values
	pos_c_5 = c_data_5[:,:3]
	vel_c_5 = c_data_5[:,3:6]
	time_c_5 = c_data_5[:,6]
	
	c_data_4 = pd.read_csv('test4/data/' + dt + 'eps4.csv').values
	pos_c_4 = c_data_4[:,:3]
	vel_c_4 = c_data_4[:,3:6]
	time_c_4 = c_data_4[:,6]

	c_data_3 = pd.read_csv('test4/data/' + dt + 'eps3.csv').values
	pos_c_3 = c_data_3[:,:3]
	vel_c_3 = c_data_3[:,3:6]
	time_c_3 = c_data_3[:,6]

	c_data_2 = pd.read_csv('test4/data/' + dt + 'eps2.csv').values
	pos_c_2 = c_data_2[:,:3]
	vel_c_2 = c_data_2[:,3:6]
	time_c_2 = c_data_2[:,6]

	c_data_1 = pd.read_csv('test4/data/' + dt + 'eps1.csv').values
	pos_c_1 = c_data_1[:,:3]
	vel_c_1 = c_data_1[:,3:6]
	time_c_1 = c_data_1[:,6]

	def energy(pos, vel, eps):
		return (0.5*(vel[:,0]**2 + vel[:,1]**2 + vel[:,2]**2)) - pos[:,0] - (eps*(pos[:,0]**2) / 2)

	return np.array([1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]), np.array([abs(max(energy(pos_c_8, vel_c_8, 1e-8))), abs(max(energy(pos_c_7, vel_c_7, 1e-7))), abs(max(energy(pos_c_6, vel_c_6, 1e-6))),abs(max(energy(pos_c_5, vel_c_5, 1e-5))), abs(max(energy(pos_c_4, vel_c_4, 1e-4))),abs(max(energy(pos_c_3, vel_c_3, 1e-3))), abs(max(energy(pos_c_2, vel_c_2, 1e-2))), abs(max(energy(pos_c_1, vel_c_1, 1e-1)))])

plt.rcParams['axes.titlesize'] = 16
plt.title("Energy Error vs ε for \nVarying Elecric and Magnetic Field")

eps1, err1 = errorPerEpsilon(1)
eps2, err2 = errorPerEpsilon(2)
eps5, err5 = errorPerEpsilon(5)

ax.plot(eps1, err1, label = "∆t = T")
ax.plot(eps2, err2, label = "∆t = T/2")
ax.plot(eps5, err5, label = "∆t = T/5")



ax.set_xscale("log")
ax.set_yscale("log")
ax.set_ylim(1e-16, 1e2)
ax.set_xlim(1e-7,1e-1)

ax.set_xlabel("ε")
ax.set_ylabel("Energy Error")
ax.legend(loc = "best", frameon=False, fontsize=12)
plt.tight_layout()
plt.savefig('test4/e_vs_eps.eps', format='eps')
plt.show()

