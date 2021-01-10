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


#Get all of the data for boris 50
boris_data_10 = pd.read_csv('test3/data/boris10.csv').values
pos_boris_10 = boris_data_10[:,:3]
vel_boris_10 = boris_data_10[:,3:6]
time_boris_10 = boris_data_10[:,6]

#Get all of the data for boris 500
boris_data_100 = pd.read_csv('test3/data/boris100.csv').values
pos_boris_100 = boris_data_100[:,:3]
vel_boris_100 = boris_data_100[:,3:6]
time_boris_100 = boris_data_100[:,6]


#Get all of the data for collocation 2
c_data_2 = pd.read_csv('test3/data/collocation2.csv').values
pos_c_2 = c_data_2[:,:3]
vel_c_2 = c_data_2[:,3:6]
time_c_2 = c_data_2[:,6]

#Get all of the data for collocation 5
c_data_5 = pd.read_csv('test3/data/collocation5.csv').values
pos_c_5 = c_data_5[:,:3]
vel_c_5 = c_data_5[:,3:6]
time_c_5 = c_data_5[:,6]

#Get all of the data for collocation 10
c_data_10 = pd.read_csv('test3/data/collocation10.csv').values
pos_c_10 = c_data_10[:,:3]
vel_c_10 = c_data_10[:,3:6]
time_c_10 = c_data_10[:,6]


# Functions for determining change in energy

def total_energy(position,velocity):

	vel_sq = velocity[0]**2 + velocity[1]**2 + velocity[2]**2
	return (0.5*(vel_sq)) - position[1]

def get_energy_error(time, pos, vel):

	E0 = total_energy(pos[0],vel[0])
	delEnergy = []

	for i in range( len( time) ):
		energy = total_energy(pos[i], vel[i])
		delEnergy.append(abs(energy - E0)/E0)

	return delEnergy



plt.xlabel("t")
plt.ylabel("Relative Change in Energy")
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

ax.plot(time_boris_10, get_energy_error(time_boris_10, pos_boris_10, vel_boris_10), label = "Boris, ∆t = T/10",linewidth=1.0,color="darkgreen")
ax.plot(time_boris_100, get_energy_error(time_boris_100, pos_boris_100, vel_boris_100), label = "Boris, ∆t = T/100",linewidth=1.0,color="orange")
ax.plot(time_c_10, get_energy_error(time_c_10, pos_c_10, vel_c_10), label = "15th Order Collocation, ∆t = T/10",linewidth=1.2)
ax.plot(time_c_5, get_energy_error(time_c_5, pos_c_5, vel_c_5), label = "15th Order Collocation, ∆t = T/5",linewidth=1.2)
ax.plot(time_c_2, get_energy_error(time_c_2, pos_c_2, vel_c_2), label = "15th Order Collocation, ∆t = T/2",linewidth=1.2,color="lightblue")

plt.yscale("log")
plt.ylim(1e-15, 1e4)

ax.legend(loc = "upper left", frameon=False, fontsize=12)
plt.tight_layout()
plt.savefig('test3/energy_errors.eps', format='eps')




