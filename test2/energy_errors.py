import matplotlib.pyplot as plt
from matplotlib import cm
import csv
import numpy as np
import pandas as pd

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

plt.gcf().subplots_adjust(left=0.15)

#Get all of the data for boris 0.01
boris_data_01 = pd.read_csv('test2/data/boris100.csv').values[::10]
pos_boris_01 = boris_data_01[:,:3]
vel_boris_01 = boris_data_01[:,3:6]
time_boris_01 = boris_data_01[:,6]

#Get all of the data for boris 0.001
boris_data_001 = pd.read_csv('test2/data/boris1000.csv').values[::100]
pos_boris_001 = boris_data_001[:,:3]
vel_boris_001 = boris_data_001[:,3:6]
time_boris_001 = boris_data_001[:,6]


#Get all of the data for collocation 0.2
c_data_5 = pd.read_csv('test2/data/collocation5.csv').values
pos_c_5 = c_data_5[:,:3]
vel_c_5 = c_data_5[:,3:6]
time_c_5 = c_data_5[:,6]

#Get all of the data for collocation 0.1
c_data_1 = pd.read_csv('test2/data/collocation10.csv').values
pos_c_10 = c_data_10[:,:3]
vel_c_10 = c_data_10[:,3:6]
time_c_10 = c_data_10[:,6]

# Functions for determining change in energy

def total_energy(position,velocity):
	vel_sq = velocity[:,0]**2 + velocity[:,1]**2 + velocity[:,2]**2
	return (0.5*(vel_sq)) - 0.1*np.cos(position[:,0])

def get_energy_error(positions, velocities):

	E0 = 0.5*(velocities[0][0]**2 + velocities[0][1]**2 + velocities[0][2]**2) - 0.1*np.cos(positions[0][0])
	
	delEnergy = total_energy(positions, velocities) - E0	
	return delEnergy


plt.xlabel("t")
plt.ylabel("Change in Energy")
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

ax.plot(time_boris_001, get_energy_error(pos_boris_001, vel_boris_001), label = "Boris, ∆t = T/1000",linewidth=1.0, color="indianred")
ax.plot(time_boris_01, get_energy_error(pos_boris_01, vel_boris_01), label = "Boris, ∆t = T/100",linewidth=1.0, color="orange")
ax.plot(time_c_10, get_energy_error(pos_c_10, vel_c_10), label = "Mixed Basis Collocation, ∆t = T/10",linewidth=1.0, color="mediumblue")
ax.plot(time_c_5, get_energy_error(pos_c_5, vel_c_5), label = "Mixed Basis Collocation, ∆t = T/5",linewidth=1.0)

ax.legend(loc = "upper left", frameon=False, fontsize=12)
plt.yscale("log")
ax.set_ylim([1e-15,1e2])
plt.savefig('test2/energy_errors.eps', format='eps')
plt.show()

