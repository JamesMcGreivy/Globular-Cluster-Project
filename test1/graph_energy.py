import matplotlib.pyplot as plt
from matplotlib import cm
import csv
import numpy as np
import pandas as pd

plt.gcf().subplots_adjust(left=0.15)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 14
plt.rcParams['axes.linewidth'] = 2
# Generate 2 colors from the 'tab10' colormap
colors = cm.get_cmap('tab10', 2)
ax = plt.subplot(111)
ax.xaxis.set_tick_params(which='major', size=7, width=2, direction='in')
ax.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in')
ax.yaxis.set_tick_params(which='major', size=7, width=2, direction='in')
ax.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in')


#Get all of the data for boris 0.01
boris_data_01 = pd.read_csv('test1/data/boris100.csv').values[::10]
pos_boris_01 = boris_data_01[:,:3]
vel_boris_01 = boris_data_01[:,3:6]
time_boris_01 = boris_data_01[:,6]

#Get all of the data for boris 0.001
boris_data_001 = pd.read_csv('test1/data/boris1000.csv').values[::100]
pos_boris_001 = boris_data_001[:,:3]
vel_boris_001 = boris_data_001[:,3:6]
time_boris_001 = boris_data_001[:,6]

#Get all of the data for collocation 0.5
c_data_2 = pd.read_csv('test1/data/collocation2.csv').values
pos_c_2 = c_data_2[:,:3]
vel_c_2 = c_data_2[:,3:6]
time_c_2 = c_data_2[:,6]


def total_energy(position, velocity):

	E0 = 0.5*(velocity[0,0]**2 + velocity[0,1]**2 + velocity[0,2]**2)

	vel_sq = velocity[:,0]**2 + velocity[:,1]**2 + velocity[:,2]**2

	return np.abs(0.5*(vel_sq) - E0)


plt.plot(time_c_2, total_energy(pos_c_2, vel_c_2), label = "Mixed-Basis Collocation, ∆t = T/2",linewidth=1.0)
plt.plot(time_boris_001, total_energy(pos_boris_001, vel_boris_001),label = "Boris, ∆t = T/1000",linewidth=1.0)
plt.plot(time_boris_01, total_energy(pos_boris_01, vel_boris_01),label = "Boris, ∆t = T/10",linewidth=1.0)

#plt.plot(time_boris_2_001,delH_boris_2_001,label = "Boris (tan), Ωc∆t = 0.001",linewidth=0.75)
#plt.plot(time_boris_2_1,delH_boris_2_1,label = "Boris (tan), Ωc∆t = 0.1",linewidth=0.75)

#plt.ylim([1e-15,10e-14])
plt.xlabel("orbits")
plt.ylabel("Relative Change in Energy")

ax.legend(loc = "upper left", frameon=False, fontsize=12)
plt.savefig('test1/energy_errors.eps', format='eps')




